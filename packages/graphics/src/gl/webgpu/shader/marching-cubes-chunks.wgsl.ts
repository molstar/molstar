/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { IsosurfaceChunkSize, IsosurfaceTableHeader } from '../../compute/marching-cubes/chunk-size.js';
import { mc_common_wgsl, mc_vertex_wgsl } from './marching-cubes.wgsl.js';

/**
 * Marching cubes into chunks of a vertex arena (see `computeChunkedIsosurface`), indexed: a vertex per crossed edge of
 * the grid, owned by the edge's lower point (its +x, +y and +z edges), shared by the triangles of the cells around
 * it; a triangle is the 3 indices of its vertices. A page of the arena holds `ChunkSize` vertices or as many
 * triangles (3 words each). The page table of a surface: its vertex count (`TableCount`), vertex page count
 * (`TablePages`), overflow status (`TableStatus`), drawn count (`TableDrawn`: 3 per triangle drawn), the box its
 * positions are quantized in (`TableBox`: min and step, f32 bits), its triangle count (`TableTriangles`), triangle
 * page count (`TableTrianglePages`) and the vertex pages it holds (`TableVertexCapacity`), then the ids of its vertex
 * pages and, after the capacity, of its triangle pages (`iso_vertex_wgsl`). The free pages are a stack: its size, then
 * the page ids.
 */
const mc_chunk_common_wgsl = /* wgsl */ `
const ChunkSize = ${IsosurfaceChunkSize}u;
const TableCount = 0u;
const TablePages = 1u;
const TableStatus = 2u;
const TableDrawn = 3u;
const TableBox = 4u;
const TableTriangles = 10u;
const TableTrianglePages = 11u;
const TableVertexCapacity = 12u;
const TableHeader = ${IsosurfaceTableHeader}u;

fn workgroupIndex(wid: vec3<u32>, count: vec3<u32>) -> u32 {
    return wid.x + wid.y * count.x + wid.z * count.x * count.y;
}
`;

/** Which own edges of a grid point the surface crosses, of the points and cells of the grid */
const mc_chunk_edges_wgsl = /* wgsl */ `
const Axes = array<vec3<i32>, 3>(vec3<i32>(1, 0, 0), vec3<i32>(0, 1, 0), vec3<i32>(0, 0, 1));

fn isInside(p: vec3<i32>) -> bool {
    return value(p) <= uParams.uIsoValue;
}

/** whether the edge from \`p\` along \`axis\` is in the grid and its ends on either side (as \`cubeCase\` takes them) */
fn isCrossed(p: vec3<i32>, inside: bool, axis: u32) -> bool {
    let q = p + Axes[axis];
    return q[axis] < uParams.uGridDim[axis] && isInside(q) != inside;
}

fn pointIndex(p: vec3<i32>) -> u32 {
    let d = vec3<u32>(uParams.uGridDim);
    let u = vec3<u32>(p);
    return u.x + u.y * d.x + u.z * d.x * d.y;
}

/**
 * \`pointBases\` of a point: its first vertex, its crossed edges (a bit per axis) in the top 3 bits; a surface has
 * fewer than 2^29 vertices (a storage buffer of 4 GiB holds 358 M)
 */
const EdgeShift = 29u;
const BaseMask = (1u << EdgeShift) - 1u;
`;

/**
 * Per grid point the vertices of its crossed edges and per cell its triangles, reserved in the surface's counts, once
 * per workgroup: the point's first vertex and its crossed edges into \`pointBases\`, the workgroup's first triangle
 * into \`bases\`
 */
export const mc_chunk_count_wgsl = /* wgsl */ `
${mc_common_wgsl}
${mc_chunk_common_wgsl}
${mc_chunk_edges_wgsl}

var<workgroup> wgVertices: atomic<u32>;
var<workgroup> wgTriangles: atomic<u32>;
var<workgroup> wgVertexBase: u32;

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(
    @builtin(global_invocation_id) id: vec3<u32>,
    @builtin(local_invocation_index) local: u32,
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(num_workgroups) workgroups: vec3<u32>,
) {
    let dim = vec3<u32>(uParams.uGridDim);
    let isPoint = all(id < dim);
    let p = vec3<i32>(id);
    var offset = 0u;
    var edges = 0u;
    if (isPoint) {
        let inside = isInside(p);
        for (var axis = 0u; axis < 3u; axis++) {
            if (isCrossed(p, inside, axis)) {
                edges |= 1u << axis;
            }
        }
        let count = countOneBits(edges);
        if (count > 0u) {
            offset = atomicAdd(&wgVertices, count);
        }
        if (all(id < dim - 1u)) {
            let triangles = triCounts[cubeCase(p)];
            if (triangles > 0u) {
                atomicAdd(&wgTriangles, triangles);
            }
        }
    }
    workgroupBarrier();
    if (local == 0u) {
        let vertices = atomicLoad(&wgVertices);
        let triangles = atomicLoad(&wgTriangles);
        var vertexBase = 0u;
        if (vertices > 0u) {
            vertexBase = atomicAdd(&table[TableCount], vertices);
        }
        var triangleBase = 0u;
        if (triangles > 0u) {
            triangleBase = atomicAdd(&table[TableTriangles], triangles);
        }
        wgVertexBase = vertexBase;
        bases[workgroupIndex(wid, workgroups)] = triangleBase;
    }
    workgroupBarrier();
    if (isPoint) {
        pointBases[pointIndex(p)] = (wgVertexBase + offset) | (edges << EdgeShift);
    }
}
`;

/**
 * The pages the surface's counts need, the vertices' first: popped from the free stack, or the surplus pushed back;
 * without enough free pages (or room in the table) the status is set. The drawn count (3 per triangle; none when
 * vertices are missing, which its triangles would index) into the draw arguments (\`drawIndirect\`, the instance count
 * is the render item's); the box of the positions and the vertex capacity into the table
 */
export const mc_chunk_allocate_wgsl = /* wgsl */ `
${mc_chunk_common_wgsl}

@compute @workgroup_size(1)
fn main() {
    let vertexCapacity = u32(uParams.uVertexCapacity);
    let triangleCapacity = arrayLength(&table) - TableHeader - vertexCapacity;
    let vertexCount = table[TableCount];
    let triangleCount = table[TableTriangles];
    let wantedVertices = (vertexCount + ChunkSize - 1u) / ChunkSize;
    let wantedTriangles = (triangleCount + ChunkSize - 1u) / ChunkSize;
    let needVertices = min(wantedVertices, vertexCapacity);
    let needTriangles = min(wantedTriangles, triangleCapacity);
    var status = select(0u, 1u, wantedVertices > vertexCapacity || wantedTriangles > triangleCapacity);
    var haveVertices = table[TablePages];
    var haveTriangles = table[TableTrianglePages];
    let triangleIds = TableHeader + vertexCapacity;
    var freeCount = free[0];
    // the surplus first, for the pages either needs more of
    for (var i = needVertices; i < haveVertices; i++) {
        freeCount += 1u;
        free[freeCount] = table[TableHeader + i];
    }
    haveVertices = min(haveVertices, needVertices);
    for (var i = needTriangles; i < haveTriangles; i++) {
        freeCount += 1u;
        free[freeCount] = table[triangleIds + i];
    }
    haveTriangles = min(haveTriangles, needTriangles);
    let takeVertices = min(needVertices - haveVertices, freeCount);
    for (var i = 0u; i < takeVertices; i++) {
        table[TableHeader + haveVertices + i] = free[freeCount - i];
    }
    freeCount -= takeVertices;
    haveVertices += takeVertices;
    let takeTriangles = min(needTriangles - haveTriangles, freeCount);
    for (var i = 0u; i < takeTriangles; i++) {
        table[triangleIds + haveTriangles + i] = free[freeCount - i];
    }
    freeCount -= takeTriangles;
    haveTriangles += takeTriangles;
    if (haveVertices < needVertices || haveTriangles < needTriangles) {
        status = 1u;
    }
    free[0] = freeCount;
    var drawn = min(triangleCount, haveTriangles * ChunkSize) * 3u;
    if (haveVertices < wantedVertices) {
        drawn = 0u;
    }
    table[TablePages] = haveVertices;
    table[TableTrianglePages] = haveTriangles;
    table[TableVertexCapacity] = vertexCapacity;
    table[TableStatus] = status;
    table[TableDrawn] = drawn;
    drawArgs[0] = drawn;
    drawArgs[2] = 0u;
    drawArgs[3] = 0u;
    for (var i = 0u; i < 3u; i++) {
        table[TableBox + i] = bitcast<u32>(uParams.uBoxMin[i]);
        table[TableBox + 3u + i] = bitcast<u32>(uParams.uBoxStep[i]);
    }
}
`;

/**
 * Per grid point the vertices of its crossed edges (\`pointBases\`), from its first vertex in the order of the axes;
 * per cell its triangles, from its workgroup's base and its offset within the workgroup (in the order the cells
 * reserve them), each the indices of the vertices of its edges; both written through the page table into
 * \`vertices\` (see \`iso_vertex_wgsl\`)
 */
export const mc_chunk_emit_wgsl = /* wgsl */ `
${mc_common_wgsl}
${mc_vertex_wgsl}
${mc_chunk_common_wgsl}
${mc_chunk_edges_wgsl}

var<workgroup> wgOffset: atomic<u32>;

/** rgb (as \`packGroup\`) as the integer of its bytes */
fn rgbBits(rgb: vec3<f32>) -> u32 {
    let v = vec3<u32>(floor(rgb * 255.0 + 0.5));
    return (v.r << 16u) | (v.g << 8u) | v.b;
}

/** the unit vector \`n\` octahedral, 12 bits per axis: x, then y in the upper 12 bits; +z for no direction */
fn octEncode(n: vec3<f32>) -> u32 {
    let l1 = abs(n.x) + abs(n.y) + abs(n.z);
    var p = vec2<f32>(0.0);
    if (l1 > 0.0) {
        p = n.xy / l1;
        if (n.z < 0.0) {
            p = (1.0 - abs(p.yx)) * select(vec2<f32>(-1.0), vec2<f32>(1.0), p >= vec2<f32>(0.0));
        }
    }
    let o = vec2<u32>(round((p * 0.5 + 0.5) * 4095.0));
    return o.x | (o.y << 12u);
}

/** the packed group of grid point \`p\`, 0 for none */
fn pointGroup(p: vec3<i32>) -> u32 {
    #if defined(dInput_buffer)
        let d = uParams.uDataDim;
        return rgbBits(packGroup(u32(dot(select(p, p - d, p >= d), uParams.uStrides))));
    #elif defined(dPackedGroup)
        return rgbBits(textureLoad(tVolume, p, 0).rgb);
    #else
        return rgbBits(packGroup(pointIndex(p)));
    #endif
}

/** the slot of item \`i\` of the pages whose ids start at \`ids\` of the table */
fn pageSlot(ids: u32, i: u32) -> u32 {
    return table[ids + i / ChunkSize] * ChunkSize + i % ChunkSize;
}

/** the vertex on the edge of \`p\` along \`axis\`, at \`index\` of the surface */
fn emitVertex(p: vec3<i32>, axis: u32, index: u32) {
    let b0 = p;
    let b1 = p + Axes[axis];
    let v0 = value(b0);
    let v1 = value(b1);
    let t = (uParams.uIsoValue - v0) / (v0 - v1);
    let p0 = vec3<f32>(b0);
    let p1 = vec3<f32>(b1);
    let position = (uParams.uTransform * vec4<f32>(p0 + t * (p0 - p1), 1.0)).xyz;

    let n0 = gradientNormal(b0);
    let n1 = gradientNormal(b1);
    var normal = -(n0 + t * (n0 - n1));
    if (uParams.uInvert != 0u) {
        normal = -normal;
    }
    normal = uParams.uTransformAdjoint * normal;
    let normalLength = length(normal);
    if (normalLength > 0.0) {
        normal /= normalLength;
    }

    // the group of the nearer end (at \`-t\` from \`b0\`), else of the other, as the CPU marching cubes; a triangle's
    // is of its first vertex (flat interpolation)
    var near = b0;
    var far = b1;
    if (-t >= 0.5) {
        near = b1;
        far = b0;
    }
    var group = pointGroup(near);
    if (group == 0u) {
        group = pointGroup(far);
    }

    let slot = pageSlot(TableHeader, index);
    let box = clamp((position - uParams.uBoxMin) * uParams.uBoxScale, vec3<f32>(0.0), vec3<f32>(65535.0));
    let q = vec3<u32>(round(box));
    let oct = octEncode(normal);
    vertices[slot * 3u] = q.x | (q.y << 16u);
    vertices[slot * 3u + 1u] = q.z | ((oct & 0xffffu) << 16u);
    vertices[slot * 3u + 2u] = group | ((oct >> 16u) << 24u);
}

/** the index of the vertex on edge \`e\` (of \`EdgeCorners\`) of the cell at \`c\`: of the edge's lower point */
fn edgeVertex(c: vec3<i32>, e: u32) -> u32 {
    let edge = EdgeCorners[e];
    let a = c + Corners[edge.x];
    let b = c + Corners[edge.y];
    let d = abs(b - a);
    let axis = select(select(2u, 1u, d.y != 0), 0u, d.x != 0);
    let point = pointBases[pointIndex(min(a, b))];
    // after the crossed edges of the axes before
    return (point & BaseMask) + countOneBits((point >> EdgeShift) & ((1u << axis) - 1u));
}

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(
    @builtin(global_invocation_id) id: vec3<u32>,
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(num_workgroups) workgroups: vec3<u32>,
) {
    let dim = vec3<u32>(uParams.uGridDim);
    if (any(id >= dim)) {
        return;
    }
    let p = vec3<i32>(id);
    // beyond the pages the surface got, see \`mc-chunk-allocate\`
    let vertexLimit = table[TablePages] * ChunkSize;
    let triangleLimit = table[TableTrianglePages] * ChunkSize;

    let point = pointBases[pointIndex(p)];
    // the set bits in turn: a loop of a varying count, not unrolled into a copy of the vertex per axis
    var edges = point >> EdgeShift;
    var index = point & BaseMask;
    while (edges != 0u && index < vertexLimit) {
        emitVertex(p, firstTrailingBit(edges), index);
        edges &= edges - 1u;
        index += 1u;
    }

    if (any(id >= dim - 1u)) {
        return;
    }
    let k = cubeCase(p);
    let count = triCounts[k];
    if (count == 0u) {
        return;
    }
    let base = bases[workgroupIndex(wid, workgroups)] + atomicAdd(&wgOffset, count);
    let triangleIds = TableHeader + table[TableVertexCapacity];
    for (var t = 0u; t < count; t++) {
        let triangle = base + t;
        if (triangle >= triangleLimit) {
            return;
        }
        let slot = pageSlot(triangleIds, triangle);
        for (var corner = 0u; corner < 3u; corner++) {
            // the same winding for negative and positive iso values
            let source = select(corner, 2u - corner, uParams.uInvert != 0u);
            vertices[slot * 3u + corner] = edgeVertex(p, triTable[k * 16u + t * 3u + source]);
        }
    }
}
`;

/** The surface's pages pushed back onto the free stack, its counts reset */
export const mc_chunk_free_wgsl = /* wgsl */ `
${mc_chunk_common_wgsl}

@compute @workgroup_size(1)
fn main() {
    var freeCount = free[0];
    let vertexPages = table[TablePages];
    for (var i = 0u; i < vertexPages; i++) {
        freeCount += 1u;
        free[freeCount] = table[TableHeader + i];
    }
    let triangleIds = TableHeader + table[TableVertexCapacity];
    let trianglePages = table[TableTrianglePages];
    for (var i = 0u; i < trianglePages; i++) {
        freeCount += 1u;
        free[freeCount] = table[triangleIds + i];
    }
    free[0] = freeCount;
    table[TableCount] = 0u;
    table[TablePages] = 0u;
    table[TableStatus] = 0u;
    table[TableDrawn] = 0u;
    table[TableTriangles] = 0u;
    table[TableTrianglePages] = 0u;
}
`;

/** The pages \`uFirst\` to \`uFirst + uCount\` (new pages of a grown arena) pushed onto the free stack */
export const mc_chunk_push_wgsl = /* wgsl */ `
@compute @workgroup_size(1)
fn main() {
    var freeCount = free[0];
    for (var i = 0u; i < u32(uParams.uCount); i++) {
        freeCount += 1u;
        free[freeCount] = u32(uParams.uFirst) + i;
    }
    free[0] = freeCount;
}
`;

/**
 * The vertices of the arena (see `computeChunkedIsosurface`) through the page table `table`, from the `u32` words
 * of `vertices`, 3 a vertex: `x | y << 16`, `z | n << 16` and `group | (n >> 16) << 24`, the position quantized
 * to 16 bits per axis in the surface's box (`TableBox` of the table), `n` the normal octahedral in 24 bits
 * (`octEncode`), the group packed as rgb (24 bits, as `packGroup`). `isoSlot` takes the `i`th vertex drawn (3 per
 * triangle, unindexed) through its triangle's indices; `isoDrawn` is the drawn count.
 */
export function iso_vertex_wgsl(table: string, vertices: string) {
  return /* wgsl */ `
fn isoSlot(i: u32) -> u32 {
    let triangle = i / 3u;
    let triangleIds = ${IsosurfaceTableHeader}u + ${table}[12u];
    let triangleSlot = ${table}[triangleIds + triangle / ${IsosurfaceChunkSize}u] * ${IsosurfaceChunkSize}u + triangle % ${IsosurfaceChunkSize}u;
    let v = ${vertices}[triangleSlot * 3u + i % 3u];
    return ${table}[${IsosurfaceTableHeader}u + v / ${IsosurfaceChunkSize}u] * ${IsosurfaceChunkSize}u + v % ${IsosurfaceChunkSize}u;
}

fn isoDrawn() -> u32 {
    return ${table}[3u];
}

fn isoPosition(slot: u32) -> vec3<f32> {
    let xy = ${vertices}[slot * 3u];
    let z = ${vertices}[slot * 3u + 1u] & 0xffffu;
    let q = vec3<f32>(f32(xy & 0xffffu), f32(xy >> 16u), f32(z));
    let boxMin = bitcast<vec3<f32>>(vec3<u32>(${table}[4u], ${table}[5u], ${table}[6u]));
    let boxStep = bitcast<vec3<f32>>(vec3<u32>(${table}[7u], ${table}[8u], ${table}[9u]));
    return boxMin + q * boxStep;
}

/** of unit length */
fn isoNormal(slot: u32) -> vec3<f32> {
    let n = (${vertices}[slot * 3u + 1u] >> 16u) | ((${vertices}[slot * 3u + 2u] >> 24u) << 16u);
    let e = vec2<f32>(f32(n & 0xfffu), f32(n >> 12u)) / 4095.0 * 2.0 - 1.0;
    var v = vec3<f32>(e, 1.0 - abs(e.x) - abs(e.y));
    let t = max(-v.z, 0.0);
    v.x += select(t, -t, v.x >= 0.0);
    v.y += select(t, -t, v.y >= 0.0);
    return normalize(v);
}

/** the group packed as rgb, the integer of its bytes: the group plus 1 (as \`unpackRGBToInt\` reads it) */
fn isoPackedGroup(slot: u32) -> u32 {
    return ${vertices}[slot * 3u + 2u] & 0xffffffu;
}
`;
}
