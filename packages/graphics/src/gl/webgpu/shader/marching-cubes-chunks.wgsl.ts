/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { IsosurfaceChunkSize } from '../../compute/marching-cubes/chunk-size.js';
import { mc_common_wgsl, mc_vertex_wgsl } from './marching-cubes.wgsl.js';

/**
 * Marching cubes into chunks of a vertex arena (see `computeChunkedIsosurface`). The page table of a surface: its
 * vertex count (`TableCount`), page count (`TablePages`), overflow status (`TableStatus`) and drawn vertex count
 * (`TableDrawn`), then its page ids; vertex `i` is at `page[i / ChunkSize] * ChunkSize + i % ChunkSize` of the arena.
 * The free pages are a stack: its size, then the page ids.
 */
const mc_chunk_common_wgsl = /* wgsl */ `
const ChunkSize = ${IsosurfaceChunkSize}u;
const TableCount = 0u;
const TablePages = 1u;
const TableStatus = 2u;
const TableDrawn = 3u;
const TableHeader = 4u;

fn workgroupIndex(wid: vec3<u32>, count: vec3<u32>) -> u32 {
    return wid.x + wid.y * count.x + wid.z * count.x * count.y;
}
`;

/** Per workgroup the vertex count of its cells, reserved in the surface's count; the workgroup's first vertex in `bases` */
export const mc_chunk_count_wgsl = /* wgsl */ `
${mc_common_wgsl}
${mc_chunk_common_wgsl}

var<workgroup> wgCount: atomic<u32>;

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(
    @builtin(global_invocation_id) id: vec3<u32>,
    @builtin(local_invocation_index) local: u32,
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(num_workgroups) workgroups: vec3<u32>,
) {
    let cells = vec3<u32>(uParams.uGridDim) - 1u;
    if (all(id < cells)) {
        let count = triCounts[cubeCase(vec3<i32>(id))] * 3u;
        if (count > 0u) {
            atomicAdd(&wgCount, count);
        }
    }
    workgroupBarrier();
    if (local == 0u) {
        let total = atomicLoad(&wgCount);
        var base = 0u;
        if (total > 0u) {
            base = atomicAdd(&table[TableCount], total);
        }
        bases[workgroupIndex(wid, workgroups)] = base;
    }
}
`;

/**
 * The pages the surface's count needs: popped from the free stack, or the surplus pushed back; without enough free
 * pages the status is set and the drawn count clamped to the pages it got. The drawn count (whole triangles) into the
 * draw arguments (\`drawIndirect\`, the instance count is the render item's)
 */
export const mc_chunk_allocate_wgsl = /* wgsl */ `
${mc_chunk_common_wgsl}

@compute @workgroup_size(1)
fn main() {
    let count = table[TableCount];
    let maxPages = arrayLength(&table) - TableHeader;
    let wanted = (count + ChunkSize - 1u) / ChunkSize;
    let need = min(wanted, maxPages);
    var status = select(0u, 1u, wanted > maxPages);
    var have = table[TablePages];
    var freeCount = free[0];
    if (need > have) {
        let take = min(need - have, freeCount);
        for (var i = 0u; i < take; i++) {
            table[TableHeader + have + i] = free[freeCount - i];
        }
        freeCount -= take;
        have += take;
        if (have < need) {
            status = 1u;
        }
    } else {
        for (var i = need; i < have; i++) {
            freeCount += 1u;
            free[freeCount] = table[TableHeader + i];
        }
        have = need;
    }
    free[0] = freeCount;
    var drawn = min(count, have * ChunkSize);
    drawn -= drawn % 3u;
    table[TablePages] = have;
    table[TableStatus] = status;
    table[TableDrawn] = drawn;
    drawArgs[0] = drawn;
    drawArgs[2] = 0u;
    drawArgs[3] = 0u;
}
`;

/**
 * The vertices of each cell, from its workgroup's base and its offset within the workgroup (in the order the cells
 * reserve them), written through the page table: position and group (packed rgba bits) in \`positions\`, the normal
 * as snorm16 in \`normals\`
 */
export const mc_chunk_emit_wgsl = /* wgsl */ `
${mc_common_wgsl}
${mc_vertex_wgsl}
${mc_chunk_common_wgsl}

var<workgroup> wgOffset: atomic<u32>;

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(
    @builtin(global_invocation_id) id: vec3<u32>,
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(num_workgroups) workgroups: vec3<u32>,
) {
    let cells = vec3<u32>(uParams.uGridDim) - 1u;
    if (any(id >= cells)) {
        return;
    }
    let c = vec3<i32>(id);
    let k = cubeCase(c);
    let count = triCounts[k] * 3u;
    if (count == 0u) {
        return;
    }
    let base = bases[workgroupIndex(wid, workgroups)] + atomicAdd(&wgOffset, count);
    let drawn = table[TableDrawn];

    // a group per cell, so triangles have a constant group
    #if defined(dInput_buffer)
        let group = vec4<f32>(packGroup(u32(dot(c, uParams.uStrides))), 1.0);
    #elif defined(dPackedGroup)
        let group = vec4<f32>(textureLoad(tVolume, c, 0).rgb, 1.0);
    #else
        let group = vec4<f32>(packGroup(cellIndex(id, cells)), 1.0);
    #endif
    let groupBits = bitcast<f32>(pack4x8unorm(group));

    for (var v = 0u; v < count; v++) {
        let i = base + v;
        // beyond the pages the surface got, see \`mc-chunk-allocate\`
        if (i >= drawn) {
            return;
        }
        // the same winding for negative and positive iso values
        var cv = v;
        if (uParams.uInvert != 0u) {
            let m = (v + 1u) % 3u;
            if (m == 1u) {
                cv = v + 2u;
            } else if (m == 0u) {
                cv = v - 2u;
            }
        }
        let edge = EdgeCorners[triTable[k * 16u + cv]];
        let b0 = c + Corners[edge.x];
        let b1 = c + Corners[edge.y];
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

        let slot = table[TableHeader + i / ChunkSize] * ChunkSize + i % ChunkSize;
        positions[slot] = vec4<f32>(position, groupBits);
        normals[slot] = vec2<u32>(pack2x16snorm(normal.xy), pack2x16snorm(vec2<f32>(normal.z, 0.0)));
    }
}
`;

/** The surface's pages pushed back onto the free stack, its counts reset */
export const mc_chunk_free_wgsl = /* wgsl */ `
${mc_chunk_common_wgsl}

@compute @workgroup_size(1)
fn main() {
    var freeCount = free[0];
    let have = table[TablePages];
    for (var i = 0u; i < have; i++) {
        freeCount += 1u;
        free[freeCount] = table[TableHeader + i];
    }
    free[0] = freeCount;
    table[TableCount] = 0u;
    table[TablePages] = 0u;
    table[TableStatus] = 0u;
    table[TableDrawn] = 0u;
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
