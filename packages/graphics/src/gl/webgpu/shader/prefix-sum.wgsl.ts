/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Exclusive scan of one block of `WorkgroupSizeX` elements per workgroup (Hillis-Steele in workgroup memory) */
export const scan_block_wgsl = /* wgsl */ `
var<workgroup> temp: array<u32, WorkgroupSizeX>;

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    // workgroups beyond the per-dimension limit are dispatched in 2d
    let block = wid.x + wid.y * u32(uParams.uGroupsX);
    let i = block * WorkgroupSizeX + lid.x;
    let count = u32(uParams.uCount);

    var value = 0u;
    if (i < count) {
        value = input[i];
    }
    temp[lid.x] = value;
    workgroupBarrier();

    for (var offset = 1u; offset < WorkgroupSizeX; offset *= 2u) {
        var t = 0u;
        if (lid.x >= offset) {
            t = temp[lid.x - offset];
        }
        workgroupBarrier();
        temp[lid.x] += t;
        workgroupBarrier();
    }

    if (i < count) {
        output[i] = temp[lid.x] - value;
    }
    if (lid.x == WorkgroupSizeX - 1u && block < u32(uParams.uBlockCount)) {
        blockSums[block] = temp[lid.x];
    }
}
`;

/** Adds the scanned total of the preceding blocks to each element */
export const scan_add_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    let block = wid.x + wid.y * u32(uParams.uGroupsX);
    let i = block * WorkgroupSizeX + lid.x;
    if (i < u32(uParams.uCount) && block < u32(uParams.uBlockCount)) {
        output[i] += blockSums[block];
    }
}
`;
