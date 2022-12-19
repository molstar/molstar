/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

/** Helper functions for manipulation with mesh data. */

import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { CIF, CifFile } from '../../mol-io/reader/cif';
import { Box3D } from '../../mol-math/geometry';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { volumeFromDensityServerData } from '../../mol-model-formats/volume/density-server';
import { Grid } from '../../mol-model/volume';
import { ColorNames } from '../../mol-util/color/names';
import { TypedArray } from '../../mol-util/type-helpers';

import { CIF_schema_mesh } from './mesh-cif-schema';


type MeshModificationParams = {
    scale?: [number, number, number],
    shift?: [number, number, number],
    matrix?: Mat4,
    group?: number,
    invertSides?: boolean
};

/** Modify mesh in-place */
export function modify(m: Mesh, params: MeshModificationParams) {
    if (params.scale !== undefined) {
        const [qx, qy, qz] = params.scale;
        const vertices = m.vertexBuffer.ref.value;
        for (let i = 0; i < vertices.length; i += 3) {
            vertices[i] *= qx;
            vertices[i + 1] *= qy;
            vertices[i + 2] *= qz;
        }
    }
    if (params.shift !== undefined) {
        const [dx, dy, dz] = params.shift;
        const vertices = m.vertexBuffer.ref.value;
        for (let i = 0; i < vertices.length; i += 3) {
            vertices[i] += dx;
            vertices[i + 1] += dy;
            vertices[i + 2] += dz;
        }
    }
    if (params.matrix !== undefined) {
        const r = m.vertexBuffer.ref.value;
        const matrix = params.matrix;
        const size = 3 * m.vertexCount;
        for (let i = 0; i < size; i += 3) {
            Vec3.transformMat4Offset(r, r, matrix, i, i, 0);
        }
    }
    if (params.group !== undefined) {
        const groups = m.groupBuffer.ref.value;
        for (let i = 0; i < groups.length; i++) {
            groups[i] = params.group;
        }
    }
    if (params.invertSides) {
        const indices = m.indexBuffer.ref.value;
        let tmp;
        for (let i = 0; i < indices.length; i += 3) {
            tmp = indices[i];
            indices[i] = indices[i + 1];
            indices[i + 1] = tmp;
        }
        const normals = m.normalBuffer.ref.value;
        for (let i = 0; i < normals.length; i++) {
            normals[i] *= -1;
        }
    }
}

/** Create a copy a mesh, possibly modified */
export function copy(m: Mesh, modification?: MeshModificationParams): Mesh {
    const nVertices = m.vertexCount;
    const nTriangles = m.triangleCount;
    const vertices = new Float32Array(m.vertexBuffer.ref.value);
    const indices = new Uint32Array(m.indexBuffer.ref.value);
    const normals = new Float32Array(m.normalBuffer.ref.value);
    const groups = new Float32Array(m.groupBuffer.ref.value);
    const result = Mesh.create(vertices, indices, normals, groups, nVertices, nTriangles);
    if (modification) {
        modify(result, modification);
    }
    return result;
}

/** Join more meshes into one */
export function concat(...meshes: Mesh[]): Mesh {
    const nVertices = sum(meshes.map(m => m.vertexCount));
    const nTriangles = sum(meshes.map(m => m.triangleCount));
    const vertices = concatArrays(Float32Array, meshes.map(m => m.vertexBuffer.ref.value));
    const normals = concatArrays(Float32Array, meshes.map(m => m.normalBuffer.ref.value));
    const groups = concatArrays(Float32Array, meshes.map(m => m.groupBuffer.ref.value));
    const newIndices = [];
    let offset = 0;
    for (const m of meshes) {
        newIndices.push(m.indexBuffer.ref.value.map(i => i + offset));
        offset += m.vertexCount;
    }
    const indices = concatArrays(Uint32Array, newIndices);
    return Mesh.create(vertices, indices, normals, groups, nVertices, nTriangles);
}

/** Return Mesh from CIF data and mesh IDs (group IDs).
 * Assume the CIF contains coords in grid space,
 * transform the output mesh to `space` */
export async function meshFromCif(data: CifFile, invertSides: boolean | undefined = undefined, outSpace: 'grid' | 'fractional' | 'cartesian' = 'cartesian'): Promise<{ mesh: Mesh, meshIds: number[] }> {
    const volumeInfoBlock = data.blocks.find(b => b.header === 'VOLUME_INFO');
    const meshesBlock = data.blocks.find(b => b.header === 'MESHES');
    if (!volumeInfoBlock || !meshesBlock) throw new Error('Missing VOLUME_INFO or MESHES block in mesh CIF file');
    const volumeInfoCif = CIF.schema.densityServer(volumeInfoBlock);
    const meshCif = CIF_schema_mesh(meshesBlock);

    const nVertices = meshCif.mesh_vertex._rowCount;
    const nTriangles = Math.floor(meshCif.mesh_triangle._rowCount / 3);

    const mesh_id = meshCif.mesh.id.toArray();
    const vertex_meshId = meshCif.mesh_vertex.mesh_id.toArray();
    const x = meshCif.mesh_vertex.x.toArray();
    const y = meshCif.mesh_vertex.y.toArray();
    const z = meshCif.mesh_vertex.z.toArray();
    const triangle_meshId = meshCif.mesh_triangle.mesh_id.toArray();
    const triangle_vertexId = meshCif.mesh_triangle.vertex_id.toArray();

    // Shift indices from within-mesh indices to overall indices
    const indices = new Uint32Array(3 * nTriangles);
    const offsets = offsetMap(vertex_meshId);
    for (let i = 0; i < 3 * nTriangles; i++) {
        const offset = offsets.get(triangle_meshId[i])!;
        indices[i] = offset + triangle_vertexId[i];
    }
    const vertices = flattenCoords(x, y, z);
    const normals = new Float32Array(3 * nVertices);
    const groups = new Float32Array(vertex_meshId);
    const mesh = Mesh.create(vertices, indices, normals, groups, nVertices, nTriangles);

    invertSides ??= isInverted(mesh);
    if (invertSides) {
        modify(mesh, { invertSides: true }); // Vertex orientation convention is opposite in Volseg API and in MolStar
    }

    if (outSpace === 'cartesian') {
        const volume = await volumeFromDensityServerData(volumeInfoCif).run();
        const gridToCartesian = Grid.getGridToCartesianTransform(volume.grid);
        modify(mesh, { matrix: gridToCartesian });
    } else if (outSpace === 'fractional') {
        const gridSize = volumeInfoCif.volume_data_3d_info.sample_count.value(0);
        const originFract = volumeInfoCif.volume_data_3d_info.origin.value(0);
        const dimensionFract = volumeInfoCif.volume_data_3d_info.dimensions.value(0);
        if (dimensionFract[0] !== 1 || dimensionFract[1] !== 1 || dimensionFract[2] !== 1) throw new Error(`Asserted the fractional dimensions are [1,1,1], but are actually [${dimensionFract}]`);
        const scale: [number, number, number] = [1 / gridSize[0], 1 / gridSize[1], 1 / gridSize[2]];
        modify(mesh, { scale: scale, shift: Array.from(originFract) as any });
    }

    Mesh.computeNormals(mesh); // normals only necessary if flatShaded==false

    // const boxMesh = makeMeshFromBox([[0,0,0], [1,1,1]], 1);
    // const gridSize = volumeInfoCif.volume_data_3d_info.sample_count.value(0); const boxMesh = makeMeshFromBox([[0,0,0], Array.from(gridSize)] as any, 1);
    // const cellSize = volumeInfoCif.volume_data_3d_info.spacegroup_cell_size.value(0); const boxMesh = makeMeshFromBox([[0, 0, 0], Array.from(cellSize)] as any, 1);
    // mesh = concat(mesh, boxMesh);  // debug
    return { mesh: mesh, meshIds: Array.from(mesh_id) };
}

function isInverted(mesh: Mesh): boolean {
    const vertices = mesh.vertexBuffer.ref.value;
    const indices = mesh.indexBuffer.ref.value;
    const center = meshCenter(mesh);
    const center3 = Vec3.create(3 * center[0], 3 * center[1], 3 * center[2]);

    let dirMetric = 0.0;
    const [a, b, c, u, v, normal, radius] = [Vec3(), Vec3(), Vec3(), Vec3(), Vec3(), Vec3(), Vec3()];
    for (let i = 0; i < indices.length; i += 3) {
        Vec3.fromArray(a, vertices, 3 * indices[i]);
        Vec3.fromArray(b, vertices, 3 * indices[i + 1]);
        Vec3.fromArray(c, vertices, 3 * indices[i + 2]);
        Vec3.sub(u, b, a);
        Vec3.sub(v, c, b);
        Vec3.cross(normal, u, v); // direction of the surface
        Vec3.add(radius, a, b);
        Vec3.add(radius, radius, c);
        Vec3.sub(radius, radius, center3); // direction center -> this triangle
        dirMetric += Vec3.dot(radius, normal);
    }
    return dirMetric < 0;
}

function meshCenter(mesh: Mesh) {
    const vertices = mesh.vertexBuffer.ref.value;
    const n = vertices.length;
    let x = 0.0;
    let y = 0.0;
    let z = 0.0;
    for (let i = 0; i < vertices.length; i += 3) {
        x += vertices[i];
        y += vertices[i + 1];
        z += vertices[i + 2];
    }
    return Vec3.create(x / n, y / n, z / n);
}

function flattenCoords(x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number>): Float32Array {
    const n = x.length;
    const out = new Float32Array(3 * n);
    for (let i = 0; i < n; i++) {
        out[3 * i] = x[i];
        out[3 * i + 1] = y[i];
        out[3 * i + 2] = z[i];
    }
    return out;
}

/** Get mapping of unique values to the position of their first occurrence */
function offsetMap(values: ArrayLike<number>) {
    const result = new Map<number, number>();
    for (let i = 0; i < values.length; i++) {
        if (!result.has(values[i])) {
            result.set(values[i], i);
        }
    }
    return result;
}

/** Return bounding box */
export function bbox(mesh: Mesh): Box3D | null { // Is there no function for this?
    const nVertices = mesh.vertexCount;
    const coords = mesh.vertexBuffer.ref.value;
    if (nVertices === 0) {
        return null;
    }
    let minX = coords[0], minY = coords[1], minZ = coords[2];
    let maxX = minX, maxY = minY, maxZ = minZ;
    for (let i = 0; i < 3 * nVertices; i += 3) {
        const x = coords[i], y = coords[i + 1], z = coords[i + 2];
        if (x < minX) minX = x;
        if (y < minY) minY = y;
        if (z < minZ) minZ = z;
        if (x > maxX) maxX = x;
        if (y > maxY) maxY = y;
        if (z > maxZ) maxZ = z;
    }
    return Box3D.create(Vec3.create(minX, minY, minZ), Vec3.create(maxX, maxY, maxZ));
}

/** Example mesh - 1 triangle */
export function fakeFakeMesh1(): Mesh {
    const nVertices = 3;
    const nTriangles = 1;
    const vertices = new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]);
    const indices = new Uint32Array([0, 1, 2]);
    const normals = new Float32Array([0, 0, 1]);
    const groups = new Float32Array([0]);
    return Mesh.create(vertices, indices, normals, groups, nVertices, nTriangles);
}

/** Example mesh - irregular tetrahedron */
export function fakeMesh4(): Mesh {
    const nVertices = 4;
    const nTriangles = 4;
    const vertices = new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1]);
    const indices = new Uint32Array([0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3]);
    const normals = new Float32Array([-1, -1, -1, 1, 0, 0, 0, 1, 0, 0, 0, 1]);
    const groups = new Float32Array([0, 1, 2, 3]);
    return Mesh.create(vertices, indices, normals, groups, nVertices, nTriangles);
}

/** Return a box-shaped mesh */
export function meshFromBox(box: [[number, number, number], [number, number, number]], group: number = 0) {
    const [[x0, y0, z0], [x1, y1, z1]] = box;
    const vertices = new Float32Array([
        x0, y0, z0,
        x1, y0, z0,
        x0, y1, z0,
        x1, y1, z0,
        x0, y0, z1,
        x1, y0, z1,
        x0, y1, z1,
        x1, y1, z1,
    ]);
    const indices = new Uint32Array([
        2, 1, 0, 1, 2, 3,
        1, 4, 0, 4, 1, 5,
        3, 5, 1, 5, 3, 7,
        2, 7, 3, 7, 2, 6,
        0, 6, 2, 6, 0, 4,
        4, 7, 6, 7, 4, 5,
    ]);
    const groups = new Float32Array([group, group, group, group, group, group, group, group]);
    const normals = new Float32Array(8);
    const mesh = Mesh.create(vertices, indices, normals, groups, 8, 12);
    Mesh.computeNormals(mesh); // normals only necessary if flatShaded==false
    return mesh;
}

function sum(array: number[]): number {
    return array.reduce((a, b) => a + b, 0);
}

function concatArrays<T extends TypedArray>(t: new (len: number) => T, arrays: T[]): T {
    const totalLength = arrays.map(a => a.length).reduce((a, b) => a + b, 0);
    const result: T = new t(totalLength);
    let offset = 0;
    for (const array of arrays) {
        result.set(array, offset);
        offset += array.length;
    }
    return result;
}

/** Generate random colors (in a cycle) */
export const ColorGenerator = function* () {
    const colors = shuffleArray(Object.values(ColorNames));
    let i = 0;
    while (true) {
        yield colors[i];
        i++;
        if (i >= colors.length) i = 0;
    }
}();
function shuffleArray<T>(array: T[]): T[] {
    // Stealed from https://www.w3docs.com/snippets/javascript/how-to-randomize-shuffle-a-javascript-array.html
    let curId = array.length;
    // There remain elements to shuffle
    while (0 !== curId) {
        // Pick a remaining element
        const randId = Math.floor(Math.random() * curId);
        curId -= 1;
        // Swap it with the current element.
        const tmp = array[curId];
        array[curId] = array[randId];
        array[randId] = tmp;
    }
    return array;
}

