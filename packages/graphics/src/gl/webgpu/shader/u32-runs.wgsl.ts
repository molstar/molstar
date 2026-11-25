/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { U32RunRecordVec4s } from '../../compute/u32-runs.js';

/** Per word (see `u32-runs.ts`): its run by a binary search of the runs' first threads, `(start + step i) | bits` */
export const u32_runs_wgsl = /* wgsl */ `
fn record(r: u32, k: u32) -> vec4<u32> {
    return runs[r * ${U32RunRecordVec4s}u + k];
}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
    // rows of workgroups past the per-dimension limit (\`getDispatchSize\`)
    let t = id.x + id.y * groups.x * WorkgroupSizeX;
    if (t >= u32(uParams.uThreadCount)) {
        return;
    }
    // the last run whose first thread is at most t
    var lo = 0u;
    var hi = u32(uParams.uRunCount) - 1u;
    while (lo < hi) {
        let mid = (lo + hi + 1u) / 2u;
        if (record(mid, 1u).y <= t) {
            lo = mid;
        } else {
            hi = mid - 1u;
        }
    }
    let a = record(lo, 0u);
    let b = record(lo, 1u);
    let i = t - b.y;
    dst[a.x + i] = (a.z + a.w * i) | b.x;
}
`;
