/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { pack_int_wgsl } from './chunks/common.wgsl.js';
import { hiz_occlusion_wgsl } from './cell-cull.wgsl.js';
import { clip_storage_wgsl } from './chunks/clip.wgsl.js';
import { lod_at_wgsl } from './chunks/lod.wgsl.js';
import { size_wgsl } from './chunks/size.wgsl.js';
import { sphere_shared_wgsl } from './chunks/sphere.wgsl.js';
import { SphereRasterLargePixels } from '../../compute/sphere-raster.js';

/** What the kernels share (`SphereRasterParams` in `sphere-raster.ts`) */
const sphere_raster_common_wgsl = /* wgsl */ `
${lod_at_wgsl}
${sphere_shared_wgsl}

/** whether level \`j\` fades and culls by distance (\`assign_lod_factor\`: a scale and a range) */
fn levelActive(j: u32) -> bool {
    let lod = uParams.uLevels[j];
    return lod.w != 0.0 && (lod.x != 0.0 || lod.y != 0.0);
}

/** the distance from the camera plane of model-space \`p\`, as \`assign_lod_factor\` */
fn cameraDistance(p: vec3<f32>) -> f32 {
    return (dot(uParams.uCameraPlane.xyz, p) + uParams.uCameraPlane.w) / uParams.uModelScale;
}
`;

/** The transforms of the instances (`transforms`) */
const sphere_raster_transform_wgsl = /* wgsl */ `
/** the transform of instance \`slot\`: its affine rows, three vec4s (\`affine-rows.ts\`) */
fn slotTransform(slot: u32) -> mat4x4<f32> {
    let o = slot * 3u;
    let r0 = transforms[o];
    let r1 = transforms[o + 1u];
    let r2 = transforms[o + 2u];
    return mat4x4<f32>(
        vec4<f32>(r0.x, r1.x, r2.x, 0.0),
        vec4<f32>(r0.y, r1.y, r2.y, 0.0),
        vec4<f32>(r0.z, r1.z, r2.z, 0.0),
        vec4<f32>(r0.w, r1.w, r2.w, 1.0)
    );
}
`;

/**
 * Appending an instance's chunks of 64 spheres: two vec4s each, (slot, first, end, level) and (\`aInstance\`, its
 * group-instance offset, its group offset, the level's scale), in the object's region of the frame's chunks. An instance
 * of one sphere takes a row from the region's end instead, read by a thread of its own (\`uSingles\`): a workgroup per
 * such sphere is mostly idle threads; a row for a chunk, so the capacity holds
 */
const sphere_raster_append_wgsl = /* wgsl */ `
fn appendChunks(slot: u32, first: u32, count: u32, level: u32, instance: vec4<u32>) {
    // the rows asked for reserve them from one count, for both ends of the region: what is written never overlaps, and
    // each end's own count stays contiguous for the raster
    let object = u32(uParams.uObjectIndex) * 5u;
    let capacity = u32(uParams.uChunkCapacity);
    if (count == 1u) {
        if (atomicAdd(&counter[object + 4u], 1u) >= capacity) { return; }
        let single = atomicAdd(&counter[object + 1u], 1u);
        let c = (u32(uParams.uChunkBase) + u32(uParams.uChunkCapacity) - 1u - single) * 2u;
        chunks[c] = vec4<u32>(slot, first, first + 1u, level);
        chunks[c + 1u] = instance;
        return;
    }
    let n = (count + 63u) / 64u;
    let reserved = atomicAdd(&counter[object + 4u], n);
    if (reserved >= capacity) { return; }
    let written = min(n, capacity - reserved);
    let base = atomicAdd(&counter[object], written);
    for (var k = 0u; k < written; k++) {
        let c = (u32(uParams.uChunkBase) + base + k) * 2u;
        chunks[c] = vec4<u32>(slot, first + k * 64u, first + min(k * 64u + 64u, count), level);
        chunks[c + 1u] = instance;
    }
}
`;

/**
 * Unmerged objects: a thread per instance (slot), the frustum and the pyramid (\`occluded\`, as merged objects' cells and
 * instances are tested) by its bounding sphere and, per level, the level's range against the sphere's distances; its
 * chunks appended per drawn level. Instances clipped by their center as the vertex stage clips them
 * (\`clip_instance_wgsl\`)
 */
export const sphere_raster_cull_wgsl = /* wgsl */ `
${sphere_raster_common_wgsl}
${hiz_occlusion_wgsl}
${sphere_raster_transform_wgsl}
${sphere_raster_append_wgsl}
${clip_storage_wgsl('clipObjects', '0u', 'uParams.uClipObjectCount')}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3u, @builtin(num_workgroups) groups: vec3u) {
    let slot = gid.x + gid.y * groups.x * WorkgroupSizeX;
    if (slot >= u32(uParams.uInstanceCount)) { return; }
    let m = slotTransform(slot);
    let bs = uParams.uBoundingSphere;
    let model = (m * vec4<f32>(bs.xyz, 1.0)).xyz;
    if (clipTest(model, 0)) { return; }
    // the instance's scale (its largest axis)
    let s = max(length(m[0].xyz), max(length(m[1].xyz), length(m[2].xyz)));
    let radius = bs.w * s * uParams.uModelScale;
    let view = (uParams.uModelView * vec4<f32>(model, 1.0)).xyz;
    if (-view.z + radius < uParams.uNear) { return; }
    if (-view.z - radius > uParams.uNear) {
        let b = spherePixelBounds(view, radius, uParams.uProjection, uParams.uViewport);
        if (b.z < 0.0 || b.w < 0.0 || b.x > uParams.uViewport.x || b.y > uParams.uViewport.y) { return; }
    }
    if (occluded(vec4<f32>(model, radius / uParams.uModelScale))) { return; }
    let d = cameraDistance(model);
    let dr = radius / uParams.uModelScale;
    let aInstance = u32(instances[slot]);
    let groupInstance = aInstance * u32(uParams.uGroupCount);
    for (var j = 0u; j < u32(uParams.uLevelCount); j++) {
        let lod = uParams.uLevels[j];
        if (levelActive(j) && uParams.uModelScale == 1.0 && (d + dr < lod.x || d - dr > lod.y)) { continue; }
        let instance = vec4<u32>(aInstance, groupInstance, 0u, bitcast<u32>(lod.w));
        appendChunks(slot, 0u, u32(uParams.uLevelCounts[j]), j, instance);
    }
}
`;

/**
 * Merged objects culled on the GPU (\`merged-culler.ts\`): a thread per slot of the culler's compacted instances, a
 * region per level and segment; the slots within a region's count (its indirect record) append the chunks of the
 * record's spheres, at the segment's scale of the level (\`tSegment\`, as \`segmentLod\`)
 */
export const sphere_raster_regions_wgsl = /* wgsl */ `
${sphere_raster_common_wgsl}
${sphere_raster_append_wgsl}

fn record(r: u32, word: u32) -> u32 {
    return records[r * u32(uParams.uRecordWords) + word];
}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3u, @builtin(num_workgroups) groups: vec3u) {
    // the slots of the raster regions, after the levels' (\`uSlotBase\`, their records after \`uRecordBase\`)
    let local = gid.x + gid.y * groups.x * WorkgroupSizeX;
    let instanceCount = u32(uParams.uInstanceCount);
    let j = local / instanceCount;
    if (j >= u32(uParams.uLevelCount)) { return; }
    let slot = u32(uParams.uSlotBase) + local;
    let recordBase = u32(uParams.uRecordBase);
    // the segment whose region holds the slot: the last whose first instance is at or below (the record's last word)
    let segmentCount = u32(uParams.uSegmentCount);
    let firstWord = u32(uParams.uRecordWords) - 1u;
    var lo = 0u;
    var hi = segmentCount - 1u;
    while (lo < hi) {
        let mid = (lo + hi + 1u) / 2u;
        if (record(recordBase + j * segmentCount + mid, firstWord) <= slot) { lo = mid; } else { hi = mid - 1u; }
    }
    let r = recordBase + j * segmentCount + lo;
    if (slot - record(r, firstWord) >= record(r, 1u)) { return; }
    // vertices and the first vertex (or index) of the quads, 6 a sphere
    let count = record(r, 0u) / 6u;
    let first = record(r, 2u) / 6u;
    let segment = segments[slot];
    let aInstance = instances[slot];
    var scale = uParams.uLevels[0].w;
    if (u32(uParams.uSegmentLevels) > 0u) {
        let o = u32(uParams.uSegmentOffset) + (u32(segment.w) * 2u + 1u) * 4u + min(j, 3u);
        scale = bitcast<f32>(segmentPage[o]);
    }
    let groupInstance = u32(segment.y + aInstance * segment.x);
    appendChunks(slot, first, count, j, vec4<u32>(u32(aInstance), groupInstance, u32(segment.z), bitcast<u32>(scale)));
}
`;

/**
 * The indirect dispatch of an object's chunks not drawn yet (65535 workgroups a row): those after the counters' words 2
 * and 3. \`uMark\`: before a phase appending more, those so far marked as drawn instead
 */
export const sphere_raster_args_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main() {
    let o = u32(uParams.uObjectIndex) * 5u;
    if (uParams.uMark != 0) {
        counter[o + 2u] = min(counter[o], u32(uParams.uChunkCapacity));
        counter[o + 3u] = counter[o + 1u];
        return;
    }
    let n = min(counter[o], u32(uParams.uChunkCapacity)) - counter[o + 2u];
    let d = u32(uParams.uObjectIndex) * 8u;
    args[d] = min(n, 65535u);
    args[d + 1u] = (n + 65534u) / 65535u;
    args[d + 2u] = 1u;
    // the instances of one sphere, 64 a workgroup
    let singles = (counter[o + 1u] - counter[o + 3u] + 63u) / 64u;
    args[d + 4u] = min(singles, 65535u);
    args[d + 5u] = (singles + 65534u) / 65535u;
    args[d + 6u] = 1u;
}
`;

/**
 * A workgroup per chunk, a thread per sphere: as the impostors size and level it (\`assign_size\`,
 * \`assign_lod_factor\`), where \`rasterizesSphere\` the pixels whose centers its front hit covers, each written with
 * \`atomicStoreMax\` of its id (\`SphereRasterTable\`: its chunk and lane) and its depth (reversed, as \`frag_depth\`).
 * Pixel clipping as the impostors' (\`sphereImpostor\`, \`dClipSphere\`): a sphere whole by its center
 * (\`dClipPrimitive\`), else per pixel the point it shows, its front, the cap of its solid interior where the front is
 * clipped (\`clipExit\`), its back where double sided, or none, the depth of that point written
 */
export const sphere_raster_wgsl = /* wgsl */ `
enable atomic_vec2u_min_max;

${pack_int_wgsl}
${size_wgsl}
${sphere_raster_common_wgsl}
${sphere_raster_transform_wgsl}
${clip_storage_wgsl('clipObjects', '0u', 'uParams.uClipObjectCount')}

/** model space (\`clip_pixel_wgsl\`'s, the model scale 1 in the raster's views) of view-space \`p\` */
fn clipModel(p: vec3<f32>) -> vec3<f32> {
    return (uParams.uInvModelView * vec4<f32>(p, 1.0)).xyz;
}

/**
 * The clip objects sphere \`i\` of \`group\` of the chunk's instance \`e\` is exempt from (\`clip_flag_wgsl\`'s flag): from
 * the clipping layer at the impostors' index (\`assign_data_index_wgsl\`), its values a byte each; 0 without one
 */
fn clipFlag(group: u32, e: vec4<u32>) -> i32 {
    if (uParams.uClippingOffset < 0) { return 0; }
    let index = select(e.y + group, e.x, uParams.uClippingInstance != 0);
    let word = clipping[u32(uParams.uClippingOffset) + (index >> 2u)];
    return i32((word >> ((index & 3u) * 8u)) & 0xffu);
}

fn readPositionGroup(i: u32) -> vec4<f32> {
    let o = u32(uParams.uPositionOffset) + i * 4u;
    return bitcast<vec4<f32>>(vec4<u32>(positions[o], positions[o + 1u], positions[o + 2u], positions[o + 3u]));
}

/** the size of sphere \`i\` of group \`group\` of the chunk's instance \`e\`, as \`assign_size\` */
fn sphereSize(i: u32, group: u32, e: vec4<u32>) -> f32 {
    var size = uParams.uSize;
    if (uParams.uSizeIndexType >= 0) {
        let instance = i32(e.x);
        let index = sizeIndex(
            uParams.uSizeIndexType,
            instance,
            i32(e.z + group),
            i32(e.y + group),
            i32(i),
            instance * uParams.uVertexCount + i32(i)
        );
        size = unpackSize(unpack4x8unorm(sizes[u32(uParams.uSizeOffset) + u32(index)]));
    }
    return size * uParams.uSizeFactor;
}

/**
 * Writes the point sphere \`id\` (at view-space \`c\` of \`radius\`) shows at pixel \`px\`, if any: its front, or with
 * \`clipPixels\` its cap, its back or none (\`clipExit\`)
 */
fn rasterPixel(px: vec2<i32>, c: vec3<f32>, radius: f32, id: u32, flag: i32, clipPixels: bool) {
    let projection = uParams.uProjection;
    let viewport = uParams.uViewport;
    // the view ray through the pixel's center
    let ndc = vec2<f32>((f32(px.x) + 0.5) / viewport.x * 2.0 - 1.0, 1.0 - (f32(px.y) + 0.5) / viewport.y * 2.0);
    let ray = normalize(vec3<f32>(
        (ndc.x + projection[2].x) / projection[0][0],
        (ndc.y + projection[2].y) / projection[1][1],
        -1.0
    ));
    let hits = raySphere(ray, c, radius);
    if (hits.x < 0.0) { return; }
    var t = hits.y;
    if (clipPixels) {
        let front = clipModel(ray * hits.y);
        if (clipTest(front, flag)) {
            let back = clipModel(ray * hits.z);
            if (uParams.uSolidInterior != 0) {
                // the first point inside the sphere no clip object removes (the front is before the near plane)
                let s = clipExit(front, back, 0.0, flag);
                if (s < 0.0) { return; }
                t = mix(hits.y, hits.z, s);
            } else if (uParams.uDoubleSided != 0 && !clipTest(back, flag)) {
                t = hits.z;
            } else {
                return;
            }
        }
    }
    let hit = ray * t;
    let clip = projection * vec4<f32>(hit, 1.0);
    // reversed depth, as \`frag_depth\`: nearer is larger, and positive floats order as their bits; beyond the far plane
    // clipped (a negative depth's bits would win)
    let depth = clip.z / clip.w;
    if (depth <= 0.0 || depth > 1.0) { return; }
    atomicStoreMax(&vis[u32(px.y) * u32(viewport.x) + u32(px.x)], vec2<u32>(id, bitcast<u32>(depth)));
}

/** Spheres covering more pixels than this are drawn by the whole workgroup, their pixels strided over its threads */
const LargeSpherePixels = ${SphereRasterLargePixels}u;

/** the workgroup's large spheres: their view-space center and radius, pixel bounds, id, clip flag and mode */
var<workgroup> wgLargeCount: atomic<u32>;
var<workgroup> wgLargeSphere: array<vec4<f32>, WorkgroupSizeX>;
var<workgroup> wgLargeBounds: array<vec4<i32>, WorkgroupSizeX>;
var<workgroup> wgLargeInfo: array<vec4<u32>, WorkgroupSizeX>;

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(workgroup_id) wid: vec3u, @builtin(local_invocation_index) lid: u32) {
    let local = wid.x + wid.y * 65535u;
    let singles = uParams.uSingles != 0;
    let objectWord = u32(uParams.uObjectIndex) * 5u;
    // after those of an earlier phase (\`sphere_raster_args_wgsl\`)
    let chunkIndex = local + counter[objectWord + 2u];
    if (!singles && chunkIndex >= min(counter[objectWord], u32(uParams.uChunkCapacity))) { return; }
    if (lid == 0u) {
        atomicStore(&wgLargeCount, 0u);
    }
    workgroupBarrier();

    // the rows of single spheres from the region's end, a row per thread (lane 0 of its id, as a chunk's first sphere)
    let single = counter[objectWord + 3u] + local * WorkgroupSizeX + lid;
    let singleValid = singles && single < min(counter[objectWord + 1u], u32(uParams.uChunkCapacity));
    let chunk = select(u32(uParams.uChunkBase) + chunkIndex, u32(uParams.uChunkBase) + u32(uParams.uChunkCapacity) - 1u - single, singles);
    let e = select(vec4<u32>(0u), chunks[chunk * 2u], !singles || singleValid);
    let i = select(e.y + lid, e.y, singles);
    // the thread's sphere, drawn here where small, queued for the workgroup where large (no return before the barrier)
    var valid = i < e.z && (!singles || singleValid);
    let instance = chunks[chunk * 2u + 1u];
    var pg = vec4<f32>(0.0);
    if (valid) { pg = readPositionGroup(i); }
    let model = (slotTransform(e.x) * vec4<f32>(pg.xyz, 1.0)).xyz;

    let clipped = uParams.uClipObjectCount > 0;
    var flag = 0;
    if (valid && clipped) { flag = clipFlag(u32(pg.w), instance); }
    if (valid && clipped && uParams.uClipPrimitive != 0 && clipTest(model, flag)) { valid = false; }
    var radius = 0.0;
    if (valid) { radius = sphereSize(i, u32(pg.w), instance) * uParams.uModelScale; }
    let scale = bitcast<f32>(instance.w);
    if (valid && levelActive(e.w)) {
        if (uParams.uModelScale != 1.0) {
            radius *= scale;
        } else {
            let level = lodAt(cameraDistance(model), uParams.uLevels[e.w].xyz, scale);
            if (level.culled) { valid = false; }
            radius *= level.factor;
        }
    }

    // clipped per pixel only where the sphere reaches a clip object's surface
    var clipPixels = false;
    if (valid && clipped && uParams.uClipPrimitive == 0) {
        let k = clipSphereClass(model, radius / uParams.uModelScale, flag);
        if (k == 2u) { valid = false; }
        clipPixels = k == 1u;
    }
    let c = (uParams.uModelView * vec4<f32>(model, 1.0)).xyz;
    let viewport = uParams.uViewport;
    var lo = vec2<i32>(0);
    var hi = vec2<i32>(-1);
    if (valid) {
        let b = spherePixelBounds(c, radius, uParams.uProjection, viewport);
        if (rasterizesSphere(c, radius, uParams.uNear, b, uParams.uMaxRadius)) {
            // the pixels whose centers lie within the bounds
            lo = max(vec2<i32>(ceil(b.xy - 0.5)), vec2<i32>(0));
            hi = min(vec2<i32>(floor(b.zw - 0.5)), vec2<i32>(viewport) - 1);
        }
    }
    let id = chunk * 64u + select(lid, 0u, singles) + 1u;
    if (all(lo <= hi)) {
        let size = hi - lo + 1;
        if (u32(size.x * size.y) > LargeSpherePixels) {
            let q = atomicAdd(&wgLargeCount, 1u);
            wgLargeSphere[q] = vec4<f32>(c, radius);
            wgLargeBounds[q] = vec4<i32>(lo, hi);
            wgLargeInfo[q] = vec4<u32>(id, bitcast<u32>(flag), select(0u, 1u, clipPixels), 0u);
        } else {
            for (var y = lo.y; y <= hi.y; y++) {
                for (var x = lo.x; x <= hi.x; x++) {
                    rasterPixel(vec2<i32>(x, y), c, radius, id, flag, clipPixels);
                }
            }
        }
    }

    // the large spheres by all threads, each a pixel at a time
    workgroupBarrier();
    let n = atomicLoad(&wgLargeCount);
    for (var q = 0u; q < n; q++) {
        let sphere = wgLargeSphere[q];
        let bounds = wgLargeBounds[q];
        let info = wgLargeInfo[q];
        let width = u32(bounds.z - bounds.x + 1);
        let total = width * u32(bounds.w - bounds.y + 1);
        for (var p = lid; p < total; p += WorkgroupSizeX) {
            let px = bounds.xy + vec2<i32>(i32(p % width), i32(p / width));
            rasterPixel(px, sphere.xyz, sphere.w, info.x, bitcast<i32>(info.y), info.z != 0u);
        }
    }
}
`;
