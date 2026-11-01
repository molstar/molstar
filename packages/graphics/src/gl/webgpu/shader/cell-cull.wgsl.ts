/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { InstanceDataVec4s } from '../../renderable/merged.js';
import { clip_storage_wgsl } from './chunks/clip.wgsl.js';

/**
 * The occlusion test of a sphere (world space, `occluded`) against a pyramid (`HiZSchema` in `cell-cull.ts`), as
 * `HiZPass.isOccluded` tests on the CPU
 */
export const hiz_occlusion_wgsl = /* wgsl */ `
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
    let levels = uParams.uHiZLevelCount;
    let minLevel = uParams.uHiZMinLevel;
    var lod = max(i32(ceil(log2(pr / 2.0))), minLevel);
    var ts = 0;
    var x0 = 0;
    var x1 = 0;
    var y0 = 0;
    var y1 = 0;
    // the texels the rectangle covers, at most two per axis: a level up where it straddles three
    loop {
        if (lod >= levels) {
            return false;
        }
        ts = 1 << u32(levels - 1 - lod);
        x0 = i32(floor(aabb.x * f32(ts)));
        x1 = i32(floor(aabb.z * f32(ts)));
        y0 = i32(floor(aabb.y * f32(ts)));
        y1 = i32(floor(aabb.w * f32(ts)));
        if (x1 - x0 <= 1 && y1 - y0 <= 1) {
            break;
        }
        lod += 1;
    }
    if (x0 < 0 || y0 < 0 || x1 >= ts || y1 >= ts) {
        return false;
    }
    // the levels from \`minLevel\` side by side: level i of size 2^(levels - 1 - i)
    var offset = 0;
    for (var i = minLevel; i < lod; i++) {
        offset += 1 << u32(levels - 1 - i);
    }
    return z <= hiZDepth(x0 + offset, y0) && z <= hiZDepth(x1 + offset, y0) && z <= hiZDepth(x0 + offset, y1) &&
        z <= hiZDepth(x1 + offset, y1);
}
`;

/** The tests of a sphere (a cell's, an instance's) both kernels make (`SphereTestSchema` in `cell-cull.ts`) */
const sphere_test_wgsl = /* wgsl */ `
${hiz_occlusion_wgsl}
${clip_storage_wgsl('segmentData', 'u32(uParams.uClipBase)', 'uParams.uClipObjectCount')}
fn outside(plane: vec4<f32>, sphere: vec4<f32>) -> bool {
    return dot(plane.xyz, sphere.xyz) + plane.w < -sphere.w;
}

/**
 * The levels (bits) a sphere is drawn at: those whose distance range it overlaps, none outside the frustum or the
 * object's distance range, behind the pyramid, or clipped whole by the clip objects (one level without levels)
 */
fn sphereLevels(sphere: vec4<f32>) -> u32 {
    if (outside(uParams.uPlane0, sphere) || outside(uParams.uPlane1, sphere) || outside(uParams.uPlane2, sphere) ||
        outside(uParams.uPlane3, sphere) || outside(uParams.uPlane4, sphere) || outside(uParams.uPlane5, sphere)) {
        return 0u;
    }
    let d = dot(uParams.uCameraPlane.xyz, sphere.xyz) + uParams.uCameraPlane.w;
    let range = uParams.uLodRange;
    if ((range.z > 0.0 && (d + sphere.w < range.x || d - sphere.w > range.y)) || occluded(sphere)) {
        return 0u;
    }
    // clipped whole by the clip objects (without a clipping layer: none exempt)
    if (uParams.uClipObjectCount > 0 && clipSphereClass(sphere.xyz, sphere.w, 0) == 2u) {
        return 0u;
    }
    let levels = u32(uParams.uLevelCount);
    if (levels == 0u) {
        return 1u;
    }
    var mask = 0u;
    for (var j = 0u; j < levels; j++) {
        let level = uParams.uLevels[j];
        if (!(d + sphere.w < level.x || d - sphere.w > level.y)) {
            mask |= 1u << j;
        }
    }
    return mask;
}
`;

/**
 * Per cell: its sphere tested against the frustum and the object's distance range as `cullSegment` tests cells; the
 * levels whose distance range it overlaps as bits (see `cell-cull.ts`)
 */
export const cell_cull_wgsl = /* wgsl */ `
${sphere_test_wgsl}
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let cell = id.x;
    if (cell >= u32(uParams.uCellCount)) {
        return;
    }
    cellLevels[cell] = sphereLevels(cellSpheres[cell]);
}
`;

/**
 * Per instance of a drawn cell: its own sphere tested as the cell's (finer: the instances of a cell are drawn at fewer
 * levels, or culled), its data copied into a slot of its segment's region at each level of both. The
 * instances of a workgroup are consecutive, mostly of one segment: those of the workgroup's first segment take their
 * slots from a workgroup counter, reserved by one global atomic per level (100k instances on a counter per segment
 * contend otherwise); the others take theirs from the global counter
 */
export const instance_copy_wgsl = /* wgsl */ `
${sphere_test_wgsl}
/** in an instance's cell: the first phase drew it */
const DrawnBit = 0x80000000u;

var<workgroup> wgCount: atomic<u32>;
var<workgroup> wgBase: u32;

/** the instance's sphere at \`o\` (\`instanceData\`) of \`segment\`: its invariant sphere transformed, by the largest scale */
fn instanceSphere(o: u32, segment: u32) -> vec4<f32> {
    let s = segmentData[segment * 3u];
    // the transform's rows (\`affine-rows.ts\`)
    let r0 = instanceData[o];
    let r1 = instanceData[o + 1u];
    let r2 = instanceData[o + 2u];
    let center = vec3<f32>(dot(r0.xyz, s.xyz) + r0.w, dot(r1.xyz, s.xyz) + r1.w, dot(r2.xyz, s.xyz) + r2.w);
    let c0 = vec3<f32>(r0.x, r1.x, r2.x);
    let c1 = vec3<f32>(r0.y, r1.y, r2.y);
    let c2 = vec3<f32>(r0.z, r1.z, r2.z);
    return vec4<f32>(center, s.w * sqrt(max(dot(c0, c0), max(dot(c1, c1), dot(c2, c2)))));
}

/**
 * Whether the sphere raster draws the instance with sphere \`sphere\` of \`segment\` at level \`j\` (\`uRasterRadius\` above
 * 0): in front of the near plane, its segment's largest sphere at the level's scale at most \`uRasterRadius\` pixels
 * across at its nearest point; it then goes into the level's raster region, not the impostors'
 */
fn rasterized(sphere: vec4<f32>, segment: u32, j: u32) -> bool {
    if (uParams.uRasterRadius <= 0.0) {
        return false;
    }
    let s = segmentData[segment * 3u];
    if (s.w <= 0.0) {
        return false;
    }
    let r = segmentData[segment * 3u + 1u].y * (sphere.w / s.w) * segmentData[segment * 3u + 2u][min(j, 3u)];
    // its spheres' centers at least at the level's minimum distance (closer ones the level culls) and in the instance
    let d = max(dot(uParams.uCameraPlane.xyz, sphere.xyz) + uParams.uCameraPlane.w - sphere.w, uParams.uLevels[j].x);
    if (d - r <= uParams.uRasterNear) {
        return false;
    }
    // the tangents from the eye, as the raster's bounds
    return r * uParams.uRasterScale <= uParams.uRasterRadius * sqrt(max(d * d - r * r, 0.0));
}

fn copyInstance(o: u32, t: u32) {
    // the transform's rows (\`affine-rows.ts\`)
    transforms[t * 3u] = instanceData[o];
    transforms[t * 3u + 1u] = instanceData[o + 1u];
    transforms[t * 3u + 2u] = instanceData[o + 2u];
    let is = instanceData[o + 3u];
    instances[t] = is.x;
    segments[t] = vec4<f32>(is.yzw, instanceData[o + 4u].x);
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
    var sphere = vec4<f32>(0.0);
    if (valid) {
        let cs = instanceCells[instance];
        let drawnFirst = (cs.x & DrawnBit) != 0u;
        mask = cellLevels[cs.x & ~DrawnBit];
        segment = cs.y;
        // the second phase: those the first left out
        if (uParams.uPhase != 0 && drawnFirst) {
            mask = 0u;
        }
        if (mask != 0u) {
            sphere = instanceSphere(instance * ${InstanceDataVec4s}u, segment);
            mask &= sphereLevels(sphere);
            // clipped by its center, as the vertex stage clips instances (\`clip_instance_wgsl\`)
            if (uParams.uClipInstance != 0 && clipTest(sphere.xyz, 0)) {
                mask = 0u;
            }
        }
        if (uParams.uPhase == 0) {
            instanceCells[instance].x = (cs.x & ~DrawnBit) | select(0u, DrawnBit, mask != 0u);
        }
    }
    let rasterBase = u32(uParams.uRasterLevels);
    // the segment of the workgroup's first instance (always valid)
    let wgSegment = instanceCells[wid.x * WorkgroupSizeX].y;
    let firstSegment = segment == wgSegment;
    let o = instance * ${InstanceDataVec4s}u;
    for (var j = 0u; j < 4u; j++) {
        var drawn = (mask & (1u << j)) != 0u;
        // the raster's instances into its regions after the levels' (by the global counters)
        if (drawn && rasterized(sphere, segment, j)) {
            let k = rasterBase + j;
            let slot = atomicAdd(&counters[k * segmentCount + segment], 1u);
            copyInstance(o, k * instanceCount + bitcast<u32>(segmentData[segment * 3u + 1u].x) + slot);
            drawn = false;
        }
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
            copyInstance(o, j * instanceCount + bitcast<u32>(segmentData[s * 3u + 1u].x) + slot);
        }
        workgroupBarrier();
    }
}
`;

/**
 * Per (level, segment) record: the instances drawn as its instance count. The second phase (\`uPhase\`): those of the
 * first phase's on in a record of \`second\`, the record the sum for the frame's later passes
 */
export const cull_args_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let k = id.x;
    if (k >= u32(uParams.uRecordCount)) {
        return;
    }
    let words = u32(uParams.uRecordWords);
    let o = k * words;
    if (uParams.uPhase != 0) {
        let first = args[o + 1u];
        second[o + 1u] = counters[k] - first;
        second[o + words - 1u] = args[o + words - 1u] + first;
    }
    args[o + 1u] = counters[k];
}
`;
