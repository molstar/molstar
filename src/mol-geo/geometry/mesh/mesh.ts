/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from '../../../mol-util';
import { Vec3, Mat4, Mat3, Vec4 } from '../../../mol-math/linear-algebra';
import { Sphere3D } from '../../../mol-math/geometry';
import { transformPositionArray, transformDirectionArray, computeIndexedVertexNormals, GroupMapping, createGroupMapping } from '../../util';
import { GeometryUtils } from '../geometry';
import { createMarkers } from '../marker-data';
import { TransformData } from '../transform-data';
import { LocationIterator, PositionLocation } from '../../util/location-iterator';
import { createColors } from '../color-data';
import { ChunkedArray, hashFnv32a, invertCantorPairing, sortedCantorPairing } from '../../../mol-data/util';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { calculateInvariantBoundingSphere, calculateTransformBoundingSphere } from '../../../mol-gl/renderable/util';
import { Theme } from '../../../mol-theme/theme';
import { MeshValues } from '../../../mol-gl/renderable/mesh';
import { Color } from '../../../mol-util/color';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';
import { createEmptyClipping } from '../clipping-data';
import { RenderableState } from '../../../mol-gl/renderable';
import { arraySetAdd } from '../../../mol-util/array';
import { degToRad } from '../../../mol-math/misc';
import { createEmptySubstance } from '../substance-data';
import { createEmptyEmissive } from '../emissive-data';

export interface Mesh {
    readonly kind: 'mesh',

    /** Number of vertices in the mesh */
    vertexCount: number,
    /** Number of triangles in the mesh */
    triangleCount: number,

    /** Vertex buffer as array of xyz values wrapped in a value cell */
    readonly vertexBuffer: ValueCell<Float32Array>,
    /** Index buffer as array of vertex index triplets wrapped in a value cell */
    readonly indexBuffer: ValueCell<Uint32Array>,
    /** Normal buffer as array of xyz values for each vertex wrapped in a value cell */
    readonly normalBuffer: ValueCell<Float32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,
    /** Indicates that group may vary within a triangle, wrapped in a value cell */
    readonly varyingGroup: ValueCell<boolean>,

    /** Bounding sphere of the mesh */
    readonly boundingSphere: Sphere3D
    /** Maps group ids to vertex indices */
    readonly groupMapping: GroupMapping

    setBoundingSphere(boundingSphere: Sphere3D): void

    readonly meta: { [k: string]: unknown }
}

export namespace Mesh {
    export function create(vertices: Float32Array, indices: Uint32Array, normals: Float32Array, groups: Float32Array, vertexCount: number, triangleCount: number, mesh?: Mesh): Mesh {
        return mesh ?
            update(vertices, indices, normals, groups, vertexCount, triangleCount, mesh) :
            fromArrays(vertices, indices, normals, groups, vertexCount, triangleCount);
    }

    export function createEmpty(mesh?: Mesh): Mesh {
        const vb = mesh ? mesh.vertexBuffer.ref.value : new Float32Array(0);
        const ib = mesh ? mesh.indexBuffer.ref.value : new Uint32Array(0);
        const nb = mesh ? mesh.normalBuffer.ref.value : new Float32Array(0);
        const gb = mesh ? mesh.groupBuffer.ref.value : new Float32Array(0);
        return create(vb, ib, nb, gb, 0, 0, mesh);
    }

    function hashCode(mesh: Mesh) {
        return hashFnv32a([
            mesh.vertexCount, mesh.triangleCount,
            mesh.vertexBuffer.ref.version, mesh.indexBuffer.ref.version,
            mesh.normalBuffer.ref.version, mesh.groupBuffer.ref.version
        ]);
    }

    function fromArrays(vertices: Float32Array, indices: Uint32Array, normals: Float32Array, groups: Float32Array, vertexCount: number, triangleCount: number): Mesh {

        const boundingSphere = Sphere3D();
        let groupMapping: GroupMapping;

        let currentHash = -1;
        let currentGroup = -1;

        const mesh = {
            kind: 'mesh' as const,
            vertexCount,
            triangleCount,
            vertexBuffer: ValueCell.create(vertices),
            indexBuffer: ValueCell.create(indices),
            normalBuffer: ValueCell.create(normals),
            groupBuffer: ValueCell.create(groups),
            varyingGroup: ValueCell.create(false),
            get boundingSphere() {
                const newHash = hashCode(mesh);
                if (newHash !== currentHash) {
                    const b = calculateInvariantBoundingSphere(mesh.vertexBuffer.ref.value, mesh.vertexCount, 1);
                    Sphere3D.copy(boundingSphere, b);
                    currentHash = newHash;
                }
                return boundingSphere;
            },
            get groupMapping() {
                if (mesh.groupBuffer.ref.version !== currentGroup) {
                    groupMapping = createGroupMapping(mesh.groupBuffer.ref.value, mesh.vertexCount);
                    currentGroup = mesh.groupBuffer.ref.version;
                }
                return groupMapping;
            },
            setBoundingSphere(sphere: Sphere3D) {
                Sphere3D.copy(boundingSphere, sphere);
                currentHash = hashCode(mesh);
            },
            meta: {}
        };
        return mesh;
    }

    function update(vertices: Float32Array, indices: Uint32Array, normals: Float32Array, groups: Float32Array, vertexCount: number, triangleCount: number, mesh: Mesh) {
        mesh.vertexCount = vertexCount;
        mesh.triangleCount = triangleCount;
        ValueCell.update(mesh.vertexBuffer, vertices);
        ValueCell.update(mesh.indexBuffer, indices);
        ValueCell.update(mesh.normalBuffer, normals);
        ValueCell.update(mesh.groupBuffer, groups);
        return mesh;
    }

    export function computeNormals(mesh: Mesh) {
        const { vertexCount, triangleCount } = mesh;
        const vertices = mesh.vertexBuffer.ref.value;
        const indices = mesh.indexBuffer.ref.value;

        const normals = mesh.normalBuffer.ref.value.length >= vertexCount * 3
            ? mesh.normalBuffer.ref.value
            : new Float32Array(vertexCount * 3);

        if (normals === mesh.normalBuffer.ref.value) {
            normals.fill(0, 0, vertexCount * 3);
        }

        computeIndexedVertexNormals(vertices, indices, normals, vertexCount, triangleCount);
        ValueCell.update(mesh.normalBuffer, normals);
    }

    export function checkForDuplicateVertices(mesh: Mesh, fractionDigits = 3) {
        const v = mesh.vertexBuffer.ref.value;

        const map = new Map<string, number>();
        const hash = (v: Vec3, d: number) => `${v[0].toFixed(d)}|${v[1].toFixed(d)}|${v[2].toFixed(d)}`;
        let duplicates = 0;

        const a = Vec3();
        for (let i = 0, il = mesh.vertexCount; i < il; ++i) {
            Vec3.fromArray(a, v, i * 3);
            const k = hash(a, fractionDigits);
            const count = map.get(k);
            if (count !== undefined) {
                duplicates += 1;
                map.set(k, count + 1);
            } else {
                map.set(k, 1);
            }
        }
        return duplicates;
    }

    const tmpMat3 = Mat3();
    export function transform(mesh: Mesh, t: Mat4) {
        const v = mesh.vertexBuffer.ref.value;
        transformPositionArray(t, v, 0, mesh.vertexCount);
        if (!Mat4.isTranslationAndUniformScaling(t)) {
            const n = Mat3.directionTransform(tmpMat3, t);
            transformDirectionArray(n, mesh.normalBuffer.ref.value, 0, mesh.vertexCount);
        }
        ValueCell.update(mesh.vertexBuffer, v);
    }

    export type OriginalData = {
        indexBuffer: Uint32Array
        vertexCount: number
        triangleCount: number
    }

    /** Meshes may contain some original data in case any processing was done. */
    export function getOriginalData(x: Mesh | MeshValues) {
        const { originalData } = 'kind' in x ? x.meta : x.meta.ref.value as Mesh['meta'];
        return originalData as OriginalData | undefined;
    }

    /**
     * Ensure that each vertices of each triangle have the same group id.
     * Note that normals are copied over and can't be re-created from the new mesh.
     */
    export function uniformTriangleGroup(mesh: Mesh, splitTriangles = true) {
        const { indexBuffer, vertexBuffer, groupBuffer, normalBuffer, triangleCount, vertexCount } = mesh;
        const ib = indexBuffer.ref.value;
        const vb = vertexBuffer.ref.value;
        const gb = groupBuffer.ref.value;
        const nb = normalBuffer.ref.value;

        // new
        const index = ChunkedArray.create(Uint32Array, 3, 1024, triangleCount);

        // re-use
        const vertex = ChunkedArray.create(Float32Array, 3, 1024, vb);
        vertex.currentIndex = vertexCount * 3;
        vertex.elementCount = vertexCount;
        const normal = ChunkedArray.create(Float32Array, 3, 1024, nb);
        normal.currentIndex = vertexCount * 3;
        normal.elementCount = vertexCount;
        const group = ChunkedArray.create(Float32Array, 1, 1024, gb);
        group.currentIndex = vertexCount;
        group.elementCount = vertexCount;

        const vi = Vec3();
        const vj = Vec3();
        const vk = Vec3();
        const ni = Vec3();
        const nj = Vec3();
        const nk = Vec3();

        function add(i: number) {
            Vec3.fromArray(vi, vb, i * 3);
            Vec3.fromArray(ni, nb, i * 3);
            ChunkedArray.add3(vertex, vi[0], vi[1], vi[2]);
            ChunkedArray.add3(normal, ni[0], ni[1], ni[2]);
        }

        function addMid(i: number, j: number) {
            Vec3.fromArray(vi, vb, i * 3);
            Vec3.fromArray(vj, vb, j * 3);
            Vec3.scale(vi, Vec3.add(vi, vi, vj), 0.5);
            Vec3.fromArray(ni, nb, i * 3);
            Vec3.fromArray(nj, nb, j * 3);
            Vec3.scale(ni, Vec3.add(ni, ni, nj), 0.5);
            ChunkedArray.add3(vertex, vi[0], vi[1], vi[2]);
            ChunkedArray.add3(normal, ni[0], ni[1], ni[2]);
        }

        function addCenter(i: number, j: number, k: number) {
            Vec3.fromArray(vi, vb, i * 3);
            Vec3.fromArray(vj, vb, j * 3);
            Vec3.fromArray(vk, vb, k * 3);
            Vec3.scale(vi, Vec3.add(vi, Vec3.add(vi, vi, vj), vk), 1 / 3);
            Vec3.fromArray(ni, nb, i * 3);
            Vec3.fromArray(nj, nb, j * 3);
            Vec3.fromArray(nk, nb, k * 3);
            Vec3.scale(ni, Vec3.add(ni, Vec3.add(ni, ni, nj), nk), 1 / 3);
            ChunkedArray.add3(vertex, vi[0], vi[1], vi[2]);
            ChunkedArray.add3(normal, ni[0], ni[1], ni[2]);
        }

        function split2(i0: number, i1: number, i2: number, g0: number, g1: number) {
            ++newTriangleCount;
            add(i0); addMid(i0, i1); addMid(i0, i2);
            ChunkedArray.add3(index, newVertexCount, newVertexCount + 1, newVertexCount + 2);
            for (let j = 0; j < 3; ++j) ChunkedArray.add(group, g0);
            newVertexCount += 3;

            newTriangleCount += 2;
            add(i1); add(i2); addMid(i0, i1); addMid(i0, i2);
            ChunkedArray.add3(index, newVertexCount, newVertexCount + 1, newVertexCount + 3);
            ChunkedArray.add3(index, newVertexCount, newVertexCount + 3, newVertexCount + 2);
            for (let j = 0; j < 4; ++j) ChunkedArray.add(group, g1);
            newVertexCount += 4;
        }

        let newVertexCount = vertexCount;
        let newTriangleCount = 0;

        if (splitTriangles) {
            for (let i = 0, il = triangleCount; i < il; ++i) {
                const i0 = ib[i * 3], i1 = ib[i * 3 + 1], i2 = ib[i * 3 + 2];
                const g0 = gb[i0], g1 = gb[i1], g2 = gb[i2];
                if (g0 === g1 && g0 === g2) {
                    ++newTriangleCount;
                    ChunkedArray.add3(index, i0, i1, i2);
                } else if (g0 === g1) {
                    split2(i2, i0, i1, g2, g0);
                } else if (g0 === g2) {
                    split2(i1, i2, i0, g1, g2);
                } else if (g1 === g2) {
                    split2(i0, i1, i2, g0, g1);
                } else {
                    newTriangleCount += 2;
                    add(i0); addMid(i0, i1); addMid(i0, i2); addCenter(i0, i1, i2);
                    ChunkedArray.add3(index, newVertexCount, newVertexCount + 1, newVertexCount + 3);
                    ChunkedArray.add3(index, newVertexCount, newVertexCount + 3, newVertexCount + 2);
                    for (let j = 0; j < 4; ++j) ChunkedArray.add(group, g0);
                    newVertexCount += 4;

                    newTriangleCount += 2;
                    add(i1); addMid(i1, i2); addMid(i1, i0); addCenter(i0, i1, i2);
                    ChunkedArray.add3(index, newVertexCount, newVertexCount + 1, newVertexCount + 3);
                    ChunkedArray.add3(index, newVertexCount, newVertexCount + 3, newVertexCount + 2);
                    for (let j = 0; j < 4; ++j) ChunkedArray.add(group, g1);
                    newVertexCount += 4;

                    newTriangleCount += 2;
                    add(i2); addMid(i2, i1); addMid(i2, i0); addCenter(i0, i1, i2);
                    ChunkedArray.add3(index, newVertexCount + 3, newVertexCount + 1, newVertexCount);
                    ChunkedArray.add3(index, newVertexCount + 2, newVertexCount + 3, newVertexCount);
                    for (let j = 0; j < 4; ++j) ChunkedArray.add(group, g2);
                    newVertexCount += 4;
                }
            }
        } else {
            for (let i = 0, il = triangleCount; i < il; ++i) {
                const i0 = ib[i * 3], i1 = ib[i * 3 + 1], i2 = ib[i * 3 + 2];
                const g0 = gb[i0], g1 = gb[i1], g2 = gb[i2];
                if (g0 !== g1 || g0 !== g2) {
                    ++newTriangleCount;
                    add(i0); add(i1); add(i2);
                    ChunkedArray.add3(index, newVertexCount, newVertexCount + 1, newVertexCount + 2);
                    const g = g1 === g2 ? g1 : g0;
                    for (let j = 0; j < 3; ++j) ChunkedArray.add(group, g);
                    newVertexCount += 3;
                } else {
                    ++newTriangleCount;
                    ChunkedArray.add3(index, i0, i1, i2);
                }
            }
        }

        const newIb = ChunkedArray.compact(index);
        const newVb = ChunkedArray.compact(vertex);
        const newNb = ChunkedArray.compact(normal);
        const newGb = ChunkedArray.compact(group);

        mesh.vertexCount = newVertexCount;
        mesh.triangleCount = newTriangleCount;

        ValueCell.update(vertexBuffer, newVb);
        ValueCell.update(groupBuffer, newGb);
        ValueCell.update(indexBuffer, newIb);
        ValueCell.update(normalBuffer, newNb);

        // keep some original data, e.g., for geometry export
        (mesh.meta.originalData as OriginalData) = { indexBuffer: ib, vertexCount, triangleCount };

        return mesh;
    }

    //

    function getNeighboursMap(mesh: Mesh) {
        const { vertexCount, triangleCount } = mesh;
        const elements = mesh.indexBuffer.ref.value;

        const neighboursMap: number[][] = [];
        for (let i = 0; i < vertexCount; ++i) {
            neighboursMap[i] = [];
        }

        for (let i = 0; i < triangleCount; ++i) {
            const v1 = elements[i * 3];
            const v2 = elements[i * 3 + 1];
            const v3 = elements[i * 3 + 2];
            arraySetAdd(neighboursMap[v1], v2);
            arraySetAdd(neighboursMap[v1], v3);
            arraySetAdd(neighboursMap[v2], v1);
            arraySetAdd(neighboursMap[v2], v3);
            arraySetAdd(neighboursMap[v3], v1);
            arraySetAdd(neighboursMap[v3], v2);
        }
        return neighboursMap;
    }

    function getEdgeCounts(mesh: Mesh) {
        const { triangleCount } = mesh;
        const elements = mesh.indexBuffer.ref.value;

        const edgeCounts = new Map<number, number>();
        const add = (a: number, b: number) => {
            const z = sortedCantorPairing(a, b);
            const c = edgeCounts.get(z) || 0;
            edgeCounts.set(z, c + 1);
        };

        for (let i = 0; i < triangleCount; ++i) {
            const a = elements[i * 3];
            const b = elements[i * 3 + 1];
            const c = elements[i * 3 + 2];
            add(a, b); add(a, c); add(b, c);
        }
        return edgeCounts;
    }

    function getBorderVertices(edgeCounts: Map<number, number>) {
        const borderVertices = new Set<number>();
        const pair: [number, number] = [0, 0];
        edgeCounts.forEach((c, z) => {
            if (c === 1) {
                invertCantorPairing(pair, z);
                borderVertices.add(pair[0]);
                borderVertices.add(pair[1]);
            }
        });

        return borderVertices;
    }

    function getBorderNeighboursMap(neighboursMap: number[][], borderVertices: Set<number>, edgeCounts: Map<number, number>) {
        const borderNeighboursMap = new Map<number, number[]>();
        const add = (v: number, nb: number) => {
            if (borderNeighboursMap.has(v)) arraySetAdd(borderNeighboursMap.get(v)!, nb);
            else borderNeighboursMap.set(v, [nb]);
        };

        borderVertices.forEach(v => {
            const neighbours = neighboursMap[v];
            for (const nb of neighbours) {
                if (borderVertices.has(nb) && edgeCounts.get(sortedCantorPairing(v, nb)) === 1) {
                    add(v, nb);
                }
            }
        });

        return borderNeighboursMap;
    }

    function trimEdges(mesh: Mesh, neighboursMap: number[][]) {
        const { indexBuffer, triangleCount } = mesh;
        const ib = indexBuffer.ref.value;

        // new
        const index = ChunkedArray.create(Uint32Array, 3, 1024, triangleCount);

        let newTriangleCount = 0;
        for (let i = 0; i < triangleCount; ++i) {
            const a = ib[i * 3];
            const b = ib[i * 3 + 1];
            const c = ib[i * 3 + 2];
            if (neighboursMap[a].length === 2 ||
                neighboursMap[b].length === 2 ||
                neighboursMap[c].length === 2) continue;

            ChunkedArray.add3(index, a, b, c);
            newTriangleCount += 1;
        }

        const newIb = ChunkedArray.compact(index);
        mesh.triangleCount = newTriangleCount;
        ValueCell.update(indexBuffer, newIb);

        return mesh;
    }

    function fillEdges(mesh: Mesh, neighboursMap: number[][], borderNeighboursMap: Map<number, number[]>, maxLengthSquared: number) {
        const { vertexBuffer, indexBuffer, normalBuffer, triangleCount } = mesh;
        const vb = vertexBuffer.ref.value;
        const ib = indexBuffer.ref.value;
        const nb = normalBuffer.ref.value;

        // new
        const index = ChunkedArray.create(Uint32Array, 3, 1024, triangleCount);

        let newTriangleCount = 0;
        for (let i = 0; i < triangleCount; ++i) {
            ChunkedArray.add3(index, ib[i * 3], ib[i * 3 + 1], ib[i * 3 + 2]);
            newTriangleCount += 1;
        }

        const vA = Vec3();
        const vB = Vec3();
        const vC = Vec3();
        const vD = Vec3();
        const vAB = Vec3();
        const vAC = Vec3();
        const vAD = Vec3();
        const vABC = Vec3();

        const vAN = Vec3();
        const vN = Vec3();

        const AngleThreshold = degToRad(120);
        const added = new Set<number>();

        const indices = Array.from(borderNeighboursMap.keys())
            .filter(v => borderNeighboursMap.get(v)!.length < 2)
            .map(v => {
                const bnd = borderNeighboursMap.get(v)!;

                Vec3.fromArray(vA, vb, v * 3);
                Vec3.fromArray(vB, vb, bnd[0] * 3);
                Vec3.fromArray(vC, vb, bnd[1] * 3);
                Vec3.sub(vAB, vB, vA);
                Vec3.sub(vAC, vC, vA);

                return [v, Vec3.angle(vAB, vAC)];
            });

        // start with the smallest angle
        indices.sort(([, a], [, b]) => a - b);

        for (const [v, angle] of indices) {
            if (added.has(v) || angle > AngleThreshold) continue;

            const nbs = borderNeighboursMap.get(v)!;
            if (neighboursMap[nbs[0]].includes(nbs[1]) &&
                !borderNeighboursMap.get(nbs[0])?.includes(nbs[1])
            ) continue;

            Vec3.fromArray(vA, vb, v * 3);
            Vec3.fromArray(vB, vb, nbs[0] * 3);
            Vec3.fromArray(vC, vb, nbs[1] * 3);
            Vec3.sub(vAB, vB, vA);
            Vec3.sub(vAC, vC, vA);
            Vec3.add(vABC, vAB, vAC);

            if (Vec3.squaredDistance(vA, vB) >= maxLengthSquared) continue;

            let add = false;
            for (const nb of neighboursMap[v]) {
                if (nbs.includes(nb)) continue;

                Vec3.fromArray(vD, vb, nb * 3);
                Vec3.sub(vAD, vD, vA);
                if (Vec3.dot(vABC, vAD) < 0) {
                    add = true;
                    break;
                }
            }
            if (!add) continue;

            Vec3.fromArray(vAN, nb, v * 3);
            Vec3.triangleNormal(vN, vA, vB, vC);
            if (Vec3.dot(vN, vAN) > 0) {
                ChunkedArray.add3(index, v, nbs[0], nbs[1]);
            } else {
                ChunkedArray.add3(index, nbs[1], nbs[0], v);
            }
            added.add(v); added.add(nbs[0]); added.add(nbs[1]);
            newTriangleCount += 1;
        }

        const newIb = ChunkedArray.compact(index);
        mesh.triangleCount = newTriangleCount;
        ValueCell.update(indexBuffer, newIb);

        return mesh;
    }

    function laplacianEdgeSmoothing(mesh: Mesh, borderNeighboursMap: Map<number, number[]>, options: { iterations: number, lambda: number }) {
        const { iterations, lambda } = options;

        const a = Vec3();
        const b = Vec3();
        const c = Vec3();
        const t = Vec3();

        const mu = -lambda;

        let dst = new Float32Array(mesh.vertexBuffer.ref.value.length);

        const step = (f: number) => {
            const pos = mesh.vertexBuffer.ref.value;
            dst.set(pos);

            borderNeighboursMap.forEach((nbs, v) => {
                if (nbs.length !== 2) return;

                Vec3.fromArray(a, pos, v * 3);
                Vec3.fromArray(b, pos, nbs[0] * 3);
                Vec3.fromArray(c, pos, nbs[1] * 3);

                const wab = 1 / Vec3.distance(a, b);
                const wac = 1 / Vec3.distance(a, c);
                Vec3.scale(b, b, wab);
                Vec3.scale(c, c, wac);

                Vec3.add(t, b, c);
                Vec3.scale(t, t, 1 / (wab + wac));
                Vec3.sub(t, t, a);

                Vec3.scale(t, t, f);
                Vec3.add(t, a, t);

                Vec3.toArray(t, dst, v * 3);
            });

            const tmp = mesh.vertexBuffer.ref.value;
            ValueCell.update(mesh.vertexBuffer, dst);
            dst = tmp;
        };

        for (let k = 0; k < iterations; ++k) {
            step(lambda);
            step(mu);
        }
    }

    export function smoothEdges(mesh: Mesh, options: { iterations: number, maxNewEdgeLength: number }) {
        trimEdges(mesh, getNeighboursMap(mesh));

        for (let k = 0; k < 10; ++k) {
            const oldTriangleCount = mesh.triangleCount;
            const edgeCounts = getEdgeCounts(mesh);
            const neighboursMap = getNeighboursMap(mesh);
            const borderVertices = getBorderVertices(edgeCounts);
            const borderNeighboursMap = getBorderNeighboursMap(neighboursMap, borderVertices, edgeCounts);
            fillEdges(mesh, neighboursMap, borderNeighboursMap, options.maxNewEdgeLength * options.maxNewEdgeLength);
            if (mesh.triangleCount === oldTriangleCount) break;
        }

        const edgeCounts = getEdgeCounts(mesh);
        const neighboursMap = getNeighboursMap(mesh);
        const borderVertices = getBorderVertices(edgeCounts);
        const borderNeighboursMap = getBorderNeighboursMap(neighboursMap, borderVertices, edgeCounts);
        laplacianEdgeSmoothing(mesh, borderNeighboursMap, { iterations: options.iterations, lambda: 0.5 });
        return mesh;
    }

    //

    export const Params = {
        ...BaseGeometry.Params,
        doubleSided: PD.Boolean(false, BaseGeometry.CustomQualityParamInfo),
        flipSided: PD.Boolean(false, BaseGeometry.ShadingCategory),
        flatShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
        ignoreLight: PD.Boolean(false, BaseGeometry.ShadingCategory),
        celShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
        xrayShaded: PD.Select<boolean | 'inverted'>(false, [[false, 'Off'], [true, 'On'], ['inverted', 'Inverted']], BaseGeometry.ShadingCategory),
        transparentBackfaces: PD.Select('off', PD.arrayToOptions(['off', 'on', 'opaque'] as const), BaseGeometry.ShadingCategory),
        bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
        bumpAmplitude: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }, BaseGeometry.ShadingCategory),
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<Mesh, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState,
        updateRenderableState,
        createPositionIterator
    };

    function createPositionIterator(mesh: Mesh, transform: TransformData): LocationIterator {
        const groupCount = mesh.vertexCount;
        const instanceCount = transform.instanceCount.ref.value;
        const location = PositionLocation();
        const p = location.position;
        const n = location.normal;
        const vs = mesh.vertexBuffer.ref.value;
        const ns = mesh.normalBuffer.ref.value;
        const m = transform.aTransform.ref.value;
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            if (instanceIndex < 0) {
                Vec3.fromArray(p, vs, groupIndex * 3);
                Vec3.fromArray(n, ns, groupIndex * 3);
            } else {
                Vec3.transformMat4Offset(p, vs, m, 0, groupIndex * 3, instanceIndex * 16);
                Vec3.transformDirectionOffset(n, ns, m, 0, groupIndex * 3, instanceIndex * 16);
            }
            return location;
        };
        return LocationIterator(groupCount, instanceCount, 1, getLocation);
    }

    function createValues(mesh: Mesh, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): MeshValues {
        const { instanceCount, groupCount } = locationIt;
        const positionIt = createPositionIterator(mesh, transform);

        const color = createColors(locationIt, positionIt, theme.color);
        const marker = props.instanceGranularity
            ? createMarkers(instanceCount, 'instance')
            : createMarkers(instanceCount * groupCount, 'groupInstance');
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const emissive = createEmptyEmissive();
        const material = createEmptySubstance();
        const clipping = createEmptyClipping();

        const counts = { drawCount: mesh.triangleCount * 3, vertexCount: mesh.vertexCount, groupCount, instanceCount };

        const invariantBoundingSphere = Sphere3D.clone(mesh.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount, 0);

        return {
            dGeometryType: ValueCell.create('mesh'),

            aPosition: mesh.vertexBuffer,
            aNormal: mesh.normalBuffer,
            aGroup: mesh.groupBuffer,
            elements: mesh.indexBuffer,
            dVaryingGroup: mesh.varyingGroup,
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),
            ...color,
            ...marker,
            ...overpaint,
            ...transparency,
            ...emissive,
            ...material,
            ...clipping,
            ...transform,

            ...BaseGeometry.createValues(props, counts),
            uDoubleSided: ValueCell.create(props.doubleSided),
            dFlatShaded: ValueCell.create(props.flatShaded),
            dFlipSided: ValueCell.create(props.flipSided),
            dIgnoreLight: ValueCell.create(props.ignoreLight),
            dCelShaded: ValueCell.create(props.celShaded),
            dXrayShaded: ValueCell.create(props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off'),
            dTransparentBackfaces: ValueCell.create(props.transparentBackfaces),
            uBumpFrequency: ValueCell.create(props.bumpFrequency),
            uBumpAmplitude: ValueCell.create(props.bumpAmplitude),

            meta: ValueCell.create(mesh.meta),
        };
    }

    function createValuesSimple(mesh: Mesh, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(mesh, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: MeshValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props);
        ValueCell.updateIfChanged(values.uDoubleSided, props.doubleSided);
        ValueCell.updateIfChanged(values.dFlatShaded, props.flatShaded);
        ValueCell.updateIfChanged(values.dFlipSided, props.flipSided);
        ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
        ValueCell.updateIfChanged(values.dCelShaded, props.celShaded);
        ValueCell.updateIfChanged(values.dXrayShaded, props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off');
        ValueCell.updateIfChanged(values.dTransparentBackfaces, props.transparentBackfaces);
        ValueCell.updateIfChanged(values.uBumpFrequency, props.bumpFrequency);
        ValueCell.updateIfChanged(values.uBumpAmplitude, props.bumpAmplitude);
    }

    function updateBoundingSphere(values: MeshValues, mesh: Mesh) {
        const invariantBoundingSphere = Sphere3D.clone(mesh.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value, 0);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
    }

    function createRenderableState(props: PD.Values<Params>): RenderableState {
        const state = BaseGeometry.createRenderableState(props);
        updateRenderableState(state, props);
        return state;
    }

    function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        BaseGeometry.updateRenderableState(state, props);
        state.opaque = state.opaque && !props.xrayShaded;
        state.writeDepth = state.opaque;
    }
}
