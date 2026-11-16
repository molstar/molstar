/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CulledPoseRecordVec4s } from '../../compute/culled-poses.js';

/**
 * Per slot (see `culled-poses.ts`): its part by a binary search of the parts' first threads, its element's position,
 * rotation and scale read as stored, the pose written as `writePoses` (`culled-data.ts`) writes it: `t = p - s R c`
 * and `s` (1 not above 0), the unit quaternion (none: the identity); the similarity folded in (`M t`, `k s`, `qm q`;
 * the identity without one, exact)
 */
export const culled_poses_wgsl = /* wgsl */ `
fn record(p: u32, k: u32) -> vec4<u32> {
    return parts[p * ${CulledPoseRecordVec4s}u + k];
}
fn recordF(p: u32, k: u32) -> vec4<f32> {
    return bitcast<vec4<f32>>(record(p, k));
}
fn source(o: u32) -> f32 {
    return bitcast<f32>(sources[o]);
}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
    // rows of workgroups past the per-dimension limit (\`getDispatchSize\`)
    let t = id.x + id.y * groups.x * WorkgroupSizeX;
    if (t >= u32(uParams.uThreadCount)) {
        return;
    }
    // the last part whose first thread is at most t
    var lo = 0u;
    var hi = u32(uParams.uPartCount) - 1u;
    while (lo < hi) {
        let mid = (lo + hi + 1u) / 2u;
        if (record(mid, 8u).y <= t) {
            lo = mid;
        } else {
            hi = mid - 1u;
        }
    }
    let p = lo;
    let r0 = record(p, 0u);
    let r1 = record(p, 1u);
    let r8 = record(p, 8u);
    let j = t - r8.y;

    var position: vec3<f32>;
    if (r0.y == 0u) {
        let o = r0.x + 3u * j;
        position = vec3<f32>(source(o), source(o + 1u), source(o + 2u));
    } else {
        // on a grid (\`GridPositions\`): 21 bits an axis in two words
        let o = r0.x + 2u * j;
        let a = sources[o];
        let b = sources[o + 1u];
        let g = recordF(p, 2u);
        let steps = vec3<u32>(a & 0x1fffffu, (a >> 21u) | ((b & 0x3ffu) << 11u), (b >> 10u) & 0x1fffffu);
        position = g.xyz + g.w * vec3<f32>(steps);
    }
    var q = vec4<f32>(0.0, 0.0, 0.0, 1.0);
    if (r0.w != 0u) {
        if (r0.w == 1u) {
            let o = r0.z + 4u * j;
            q = vec4<f32>(source(o), source(o + 1u), source(o + 2u), source(o + 3u));
        } else {
            let o = r0.z + 2u * j;
            q = vec4<f32>(unpack2x16snorm(sources[o]), unpack2x16snorm(sources[o + 1u]));
        }
        let n = length(q);
        q = select(vec4<f32>(0.0, 0.0, 0.0, 1.0), q / n, n > 0.0);
    }
    var scale = bitcast<f32>(r1.z);
    if (r1.y != 0u) {
        scale = source(r1.x + j);
    }
    if (!(scale > 0.0)) {
        scale = 1.0;
    }

    // R c: c + 2 w (q × c) + 2 q × (q × c)
    let c = recordF(p, 3u);
    let tc = 2.0 * cross(q.xyz, c.xyz);
    let rc = c.xyz + q.w * tc + cross(q.xyz, tc);
    let base = position - scale * rc;
    let m0 = recordF(p, 5u);
    let m1 = recordF(p, 6u);
    let m2 = recordF(p, 7u);
    let moved = vec3<f32>(dot(m0.xyz, base) + m0.w, dot(m1.xyz, base) + m1.w, dot(m2.xyz, base) + m2.w);
    // qm q (Hamilton product: R(qm q) = Rm R)
    let a = recordF(p, 4u);
    let rotation = vec4<f32>(a.w * q.xyz + q.w * a.xyz + cross(a.xyz, q.xyz), a.w * q.w - dot(a.xyz, q.xyz));

    let o = r8.x + 2u * j;
    data[o] = vec4<f32>(moved, c.w * scale);
    data[o + 1u] = rotation;
}
`;
