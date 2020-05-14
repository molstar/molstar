/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from '../../../mol-util';
import { Vec3, Mat4, Mat3, Vec4 } from '../../../mol-math/linear-algebra';
import { Sphere3D } from '../../../mol-math/geometry';
import { transformPositionArray, transformDirectionArray, computeIndexedVertexNormals, GroupMapping, createGroupMapping} from '../../util';
import { GeometryUtils } from '../geometry';
import { createMarkers } from '../marker-data';
import { TransformData } from '../transform-data';
import { LocationIterator } from '../../util/location-iterator';
import { createColors } from '../color-data';
import { ChunkedArray, hashFnv32a } from '../../../mol-data/util';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { calculateInvariantBoundingSphere, calculateTransformBoundingSphere } from '../../../mol-gl/renderable/util';
import { Theme } from '../../../mol-theme/theme';
import { MeshValues } from '../../../mol-gl/renderable/mesh';
import { Color } from '../../../mol-util/color';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';
import { createEmptyClipping } from '../clipping-data';

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

    /** Bounding sphere of the mesh */
    readonly boundingSphere: Sphere3D
    /** Maps group ids to vertex indices */
    readonly groupMapping: GroupMapping

    setBoundingSphere(boundingSphere: Sphere3D): void
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
            }
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

        ValueCell.update(vertexBuffer, newVb) as ValueCell<Float32Array>;
        ValueCell.update(groupBuffer, newGb) as ValueCell<Float32Array>;
        ValueCell.update(indexBuffer, newIb) as ValueCell<Uint32Array>;
        ValueCell.update(normalBuffer, newNb) as ValueCell<Float32Array>;

        return mesh;
    }

    //

    export const Params = {
        ...BaseGeometry.Params,
        doubleSided: PD.Boolean(false, BaseGeometry.CustomQualityParamInfo),
        flipSided: PD.Boolean(false, BaseGeometry.ShadingCategory),
        flatShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
        ignoreLight: PD.Boolean(false, BaseGeometry.ShadingCategory),
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<Mesh, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState: BaseGeometry.createRenderableState,
        updateRenderableState: BaseGeometry.updateRenderableState
    };

    function createValues(mesh: Mesh, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): MeshValues {
        const { instanceCount, groupCount } = locationIt;
        if (instanceCount !== transform.instanceCount.ref.value) {
            throw new Error('instanceCount values in TransformData and LocationIterator differ');
        }

        const color = createColors(locationIt, theme.color);
        const marker = createMarkers(instanceCount * groupCount);
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const clipping = createEmptyClipping();

        const counts = { drawCount: mesh.triangleCount * 3, groupCount, instanceCount };

        const invariantBoundingSphere = Sphere3D.clone(mesh.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount);

        return {
            aPosition: mesh.vertexBuffer,
            aNormal: mesh.normalBuffer,
            aGroup: mesh.groupBuffer,
            elements: mesh.indexBuffer,
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),
            ...color,
            ...marker,
            ...overpaint,
            ...transparency,
            ...clipping,
            ...transform,

            ...BaseGeometry.createValues(props, counts),
            dDoubleSided: ValueCell.create(props.doubleSided),
            dFlatShaded: ValueCell.create(props.flatShaded),
            dFlipSided: ValueCell.create(props.flipSided),
            dIgnoreLight: ValueCell.create(props.ignoreLight),
        };
    }

    function createValuesSimple(mesh: Mesh, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(mesh, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: MeshValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props);
        ValueCell.updateIfChanged(values.dDoubleSided, props.doubleSided);
        ValueCell.updateIfChanged(values.dFlatShaded, props.flatShaded);
        ValueCell.updateIfChanged(values.dFlipSided, props.flipSided);
        ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
    }

    function updateBoundingSphere(values: MeshValues, mesh: Mesh) {
        const invariantBoundingSphere = Sphere3D.clone(mesh.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
    }
}
