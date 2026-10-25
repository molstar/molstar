/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { InstanceDataVec4s } from '../../renderable/merged.js';

/**
 * Per cell: its sphere tested against the frustum and the object's distance range as `cullSegment` tests cells; the
 * levels whose distance range it overlaps as bits (see `cell-cull.ts`)
 */
export const cell_cull_wgsl = /* wgsl */ `
fn outside(plane: vec4<f32>, sphere: vec4<f32>) -> bool {
    return dot(plane.xyz, sphere.xyz) + plane.w < -sphere.w;
}

// the occlusion test of \`HiZPass.isOccluded\` against the pyramid of the previous frame (depth reversed: near at 1)
fn hiZViewZ(depth: f32) -> f32 {
    let near = uParams.uHiZDepthRange.x;
    let far = uParams.uHiZDepthRange.y;
    if (uParams.uHiZProjection[2][3] == -1.0) {
        return -(near * far) / (near + (far - near) * depth);
    }
    return depth * (far - near) - far;
}

/** the screen rectangle (min x, min y, max x, max y in [0, 1]) of a sphere at \`p\` (view space) */
fn hiZProjectSphere(p: vec3<f32>, r: f32) -> vec4<f32> {
    let P = uParams.uHiZProjection;
    if (P[2][3] == -1.0) {
        // Mara and McGuire 2013, as specialized by Kapoulkine (niagara)
        let pr = p * r;
        let pzr2 = p.z * p.z - r * r;
        let vx = sqrt(p.x * p.x + pzr2);
        let minx = ((vx * p.x - pr.z) / (vx * p.z + pr.x)) * P[0][0];
        let maxx = ((vx * p.x + pr.z) / (vx * p.z - pr.x)) * P[0][0];
        let vy = sqrt(p.y * p.y + pzr2);
        let miny = ((vy * p.y - pr.z) / (vy * p.z + pr.y)) * P[1][1];
        let maxy = ((vy * p.y + pr.z) / (vy * p.z - pr.y)) * P[1][1];
        return vec4<f32>(maxx * -0.5 + 0.5, miny * 0.5 + 0.5, minx * -0.5 + 0.5, maxy * 0.5 + 0.5);
    }
    let minx = (p.x + r) * P[0][0];
    let maxx = (p.x - r) * P[0][0];
    let miny = (p.y + r) * P[1][1];
    let maxy = (p.y - r) * P[1][1];
    return vec4<f32>(maxx * 0.5 + 0.5, miny * -0.5 + 0.5, minx * 0.5 + 0.5, maxy * -0.5 + 0.5);
}

fn hiZDepth(x: i32, y: i32) -> f32 {
    return hiZViewZ(textureLoad(tHiZ, vec2<i32>(x, y), 0).r);
}

fn occluded(sphere: vec4<f32>) -> bool {
    if (uParams.uHiZ == 0u) {
        return false;
    }
    let vp = (uParams.uHiZView * vec4<f32>(sphere.xyz, 1.0)).xyz;
    let r = sphere.w * 1.2 + 1.52;
    let z = vp.z + r;
    if (-z < uParams.uHiZDepthRange.x) {
        return false;
    }
    let aabb = hiZProjectSphere(vp, r);
    let w = aabb.z - aabb.x;
    let h = aabb.w - aabb.y;
    let pr = max(w * uParams.uHiZViewport.x, h * uParams.uHiZViewport.y);
    let lod = i32(ceil(log2(pr / 2.0)));
    let levels = uParams.uHiZLevelCount;
    let minLevel = uParams.uHiZMinLevel;
    if (lod >= levels || lod < minLevel) {
        return false;
    }
    // the levels from \`minLevel\` side by side: level i of size 2^(levels - 1 - i)
    var offset = 0;
    for (var i = minLevel; i < lod; i++) {
        offset += 1 << u32(levels - 1 - i);
    }
    let ts = 1 << u32(levels - 1 - lod);
    // the texels the rectangle covers: at most two per axis at this level (its size at most a texel)
    let x0 = i32(floor(aabb.x * f32(ts)));
    let x1 = i32(floor(aabb.z * f32(ts)));
    let y0 = i32(floor(aabb.y * f32(ts)));
    let y1 = i32(floor(aabb.w * f32(ts)));
    if (x0 < 0 || y0 < 0 || x1 >= ts || y1 >= ts || x1 - x0 > 1 || y1 - y0 > 1) {
        return false;
    }
    return z <= hiZDepth(x0 + offset, y0) && z <= hiZDepth(x1 + offset, y0) && z <= hiZDepth(x0 + offset, y1) &&
        z <= hiZDepth(x1 + offset, y1);
}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let cell = id.x;
    if (cell >= u32(uParams.uCellCount)) {
        return;
    }
    let sphere = cellSpheres[cell];
    if (outside(uParams.uPlane0, sphere) || outside(uParams.uPlane1, sphere) || outside(uParams.uPlane2, sphere) ||
        outside(uParams.uPlane3, sphere) || outside(uParams.uPlane4, sphere) || outside(uParams.uPlane5, sphere)) {
        cellLevels[cell] = 0u;
        return;
    }
    let d = dot(uParams.uCameraPlane.xyz, sphere.xyz) + uParams.uCameraPlane.w;
    let range = uParams.uLodRange;
    if ((range.z > 0.0 && (d + sphere.w < range.x || d - sphere.w > range.y)) || occluded(sphere)) {
        cellLevels[cell] = 0u;
        return;
    }
    let levels = u32(uParams.uLevelCount);
    if (levels == 0u) {
        cellLevels[cell] = 1u;
        return;
    }
    var mask = 0u;
    for (var j = 0u; j < levels; j++) {
        let level = uParams.uLevels[j];
        if (!(d + sphere.w < level.x || d - sphere.w > level.y)) {
            mask |= 1u << j;
        }
    }
    cellLevels[cell] = mask;
}
`;

/**
 * Per instance of a drawn cell: its data copied into a slot of its segment's region at each level of the cell. The
 * instances of a workgroup are consecutive, mostly of one segment: those of the workgroup's first segment take their
 * slots from a workgroup counter, reserved by one global atomic per level (100k instances on a counter per segment
 * contend otherwise); the others take theirs from the global counter
 */
export const instance_copy_wgsl = /* wgsl */ `
var<workgroup> wgCount: atomic<u32>;
var<workgroup> wgBase: u32;

fn copyInstance(o: u32, t: u32) {
    transforms[t * 4u] = instanceData[o];
    transforms[t * 4u + 1u] = instanceData[o + 1u];
    transforms[t * 4u + 2u] = instanceData[o + 2u];
    transforms[t * 4u + 3u] = instanceData[o + 3u];
    let is = instanceData[o + 4u];
    instances[t] = is.x;
    segments[t * 3u] = is.y;
    segments[t * 3u + 1u] = is.z;
    segments[t * 3u + 2u] = is.w;
    spheres[t] = instanceData[o + 5u];
    lods[t] = instanceData[o + 6u];
}

@compute @workgroup_size(WorkgroupSizeX)
fn main(
    @builtin(global_invocation_id) id: vec3<u32>,
    @builtin(local_invocation_index) local: u32,
    @builtin(workgroup_id) wid: vec3<u32>,
) {
    let instance = id.x;
    let instanceCount = u32(uParams.uInstanceCount);
    let segmentCount = u32(uParams.uSegmentCount);
    let valid = instance < instanceCount;
    var mask = 0u;
    var segment = 0u;
    if (valid) {
        let cs = instanceCells[instance];
        mask = cellLevels[cs.x];
        segment = cs.y;
    }
    // the segment of the workgroup's first instance (always valid)
    let wgSegment = instanceCells[wid.x * WorkgroupSizeX].y;
    let firstSegment = segment == wgSegment;
    let o = instance * ${InstanceDataVec4s}u;
    for (var j = 0u; j < 4u; j++) {
        let drawn = (mask & (1u << j)) != 0u;
        if (local == 0u) {
            atomicStore(&wgCount, 0u);
        }
        workgroupBarrier();
        var slot = 0u;
        if (drawn && firstSegment) {
            slot = atomicAdd(&wgCount, 1u);
        }
        workgroupBarrier();
        if (local == 0u) {
            let total = atomicLoad(&wgCount);
            if (total > 0u) {
                wgBase = atomicAdd(&counters[j * segmentCount + wgSegment], total);
            }
        }
        workgroupBarrier();
        if (drawn) {
            if (firstSegment) {
                slot += wgBase;
            } else {
                slot = atomicAdd(&counters[j * segmentCount + segment], 1u);
            }
            let s = select(segment, wgSegment, firstSegment);
            copyInstance(o, j * instanceCount + segmentBases[s] + slot);
        }
        workgroupBarrier();
    }
}
`;

/** Per (level, segment) record: the instances drawn as its instance count */
export const cull_args_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let k = id.x;
    if (k >= u32(uParams.uRecordCount)) {
        return;
    }
    args[k * u32(uParams.uRecordWords) + 1u] = counters[k];
}
`;
