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
import {
  SphereRasterArgsWords,
  SphereRasterChunksPerWorkgroup,
  SphereRasterCounterWords,
  SphereRasterLargePixels,
  SphereRasterWideChunksPerWorkgroup,
  SphereRasterWideSpheres,
} from '../../compute/sphere-raster.js';
import { CulledData } from '../../culled-data.js';
import { SegmentTexels } from '../../renderable/schema.js';

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

/** The transforms of the instances (`transforms`, `CulledData`: a culler's, or an unmerged object's own) */
const sphere_raster_transform_wgsl = /* wgsl */ `
${CulledData.wgsl('transforms', 'slot')}

/**
 * the transform of instance \`slot\`: its affine rows (\`affine-rows.ts\`, a pose's built), the object's transient
 * transform after it
 */
fn slotTransform(slot: u32) -> mat4x4<f32> {
    let rows = slotRows(slot);
    let r0 = rows[0];
    let r1 = rows[1];
    let r2 = rows[2];
    return uParams.uObjectTransform * mat4x4<f32>(
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
/** the chunks of \`count\` items from \`first\`, every \`stride\`-th (a level's), a chunk 64 of them */
fn appendChunks(slot: u32, first: u32, count: u32, level: u32, stride: u32, instance: vec4<u32>) {
    let levelStride = level | (stride << 2u);
    // the rows asked for reserve them from one count, for both ends of the region: what is written never overlaps, and
    // each end's own count stays contiguous for the raster
    let object = u32(uParams.uObjectIndex) * ${SphereRasterCounterWords}u;
    let capacity = u32(uParams.uChunkCapacity);
    if (count == 1u) {
        if (atomicAdd(&counter[object + 4u], 1u) >= capacity) { return; }
        let single = atomicAdd(&counter[object + 1u], 1u);
        let c = (u32(uParams.uChunkBase) + u32(uParams.uChunkCapacity) - 1u - single) * 2u;
        chunks[c] = vec4<u32>(slot, first, 1u, levelStride);
        chunks[c + 1u] = instance;
        return;
    }
    let n = (count + 63u) / 64u;
    let reserved = atomicAdd(&counter[object + 4u], n);
    if (reserved >= capacity) { return; }
    let written = min(n, capacity - reserved);
    let base = atomicAdd(&counter[object], written);
    atomicAdd(&counter[object + 5u], count);
    for (var k = 0u; k < written; k++) {
        let c = (u32(uParams.uChunkBase) + base + k) * 2u;
        chunks[c] = vec4<u32>(slot, first + k * 64u * stride, min(64u, count - k * 64u), levelStride);
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
    let aInstance = instances[slot];
    let groupInstance = aInstance * u32(uParams.uGroupCount);
    for (var j = 0u; j < u32(uParams.uLevelCount); j++) {
        let lod = uParams.uLevels[j];
        if (levelActive(j) && (d + dr < lod.x || d - dr > lod.y)) { continue; }
        let instance = vec4<u32>(aInstance, groupInstance, 0u, bitcast<u32>(lod.w));
        appendChunks(slot, 0u, u32(uParams.uLevelCounts[j]), j, u32(uParams.uLevelStrides[j]), instance);
    }
}
`;

/**
 * Merged objects culled on the GPU (\`merged-culler.ts\`): a workgroup a block of a raster record's instances (its
 * slots of the culler's regions of instance indices; \`blockFirsts\` a record's first block, the blocks' total after
 * them, \`sphere-raster-blocks\`), not a thread a slot of every level's regions (most in no raster region: the
 * yeast's 126 M threads 2.7 ms a phase for nothing); its slots append the chunks of the record's spheres, at the
 * segment's scale of the level (\`tSegment\`, as \`segmentLod\`), by the instance the slot holds (the chunk's slot: an
 * index of \`CulledData\`, whose transforms come first, as the raster reads them)
 */
export const sphere_raster_regions_wgsl = /* wgsl */ `
${sphere_raster_common_wgsl}
${sphere_raster_append_wgsl}
${CulledData.wgsl('data')}

fn record(r: u32, word: u32) -> u32 {
    return records[r * u32(uParams.uRecordWords) + word];
}

/** the workgroup's raster record (from the first, \`uRecordBase\`) */
var<workgroup> wgRecord: u32;

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(workgroup_id) wid: vec3u, @builtin(num_workgroups) groups: vec3u,
        @builtin(local_invocation_index) lid: u32) {
    let w = wid.x + wid.y * groups.x;
    let n = u32(uParams.uRecordCount);
    // its record: the last whose first block is at or below (one of no blocks the next's first: past it)
    if (lid == 0u) {
        var lo = 0u;
        var hi = n - 1u;
        while (lo < hi) {
            let mid = (lo + hi + 1u) / 2u;
            if (blockFirsts[mid] <= w) { lo = mid; } else { hi = mid - 1u; }
        }
        wgRecord = lo;
    }
    let rr = workgroupUniformLoad(&wgRecord);
    if (w >= blockFirsts[n]) { return; }
    let r = u32(uParams.uRecordBase) + rr;
    let k = (w - blockFirsts[rr]) * WorkgroupSizeX + lid;
    if (k >= record(r, 1u)) { return; }
    let slot = record(r, u32(uParams.uRecordWords) - 1u) + k;
    let j = rr / u32(uParams.uSegmentCount);
    // vertices and the first vertex (or index) of the quads, 6 a sphere
    let count = record(r, 0u) / 6u;
    // the first vertex id: an indexed record's base vertex (its first index 0, merged-culler.ts), else its first vertex
    let first = select(record(r, 2u), record(r, 3u), u32(uParams.uRecordWords) == 5u) / 6u;
    let i = indices[slot];
    let segment = culledSegment(i);
    let aInstance = culledInstance(i);
    var scale = uParams.uLevels[0].w;
    // the segment's every stride-th item at the level (1 without levels, or drawn whole), from its first
    var stride = 1u;
    if (u32(uParams.uSegmentLevels) > 0u) {
        let o = u32(uParams.uSegmentOffset) + (u32(segment.w) * ${SegmentTexels}u + 1u) * 4u + min(j, 3u);
        scale = bitcast<f32>(segmentPage[o]);
        stride = u32(bitcast<f32>(segmentPage[o + 4u]));
    }
    let groupInstance = u32(segment.y + i32(aInstance) * segment.x);
    appendChunks(i, first, count, j, stride, vec4<u32>(aInstance, groupInstance, u32(segment.z), bitcast<u32>(scale)));
}
`;

/**
 * The indirect dispatches of an object's chunks not drawn yet (65535 workgroups a row; \`SphereRasterArgsWords\` words
 * each): those after the counters' words 2 and 3, by one of the two raster kernels (\`SphereRasterChunksPerWorkgroup\`),
 * the other none; the singles. \`uMark\`: before a phase appending more, those so far marked as drawn instead
 */
export const sphere_raster_args_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main() {
    let o = u32(uParams.uObjectIndex) * ${SphereRasterCounterWords}u;
    if (uParams.uMark != 0) {
        counter[o + 2u] = min(counter[o], u32(uParams.uChunkCapacity));
        counter[o + 3u] = counter[o + 1u];
        return;
    }
    // the wide kernel where the chunks hold few spheres on average (the spheres of all appended so far)
    let chunkCount = min(counter[o], u32(uParams.uChunkCapacity));
    let wide = counter[o + 5u] <= ${SphereRasterWideSpheres}u * chunkCount;
    let g = select(${SphereRasterChunksPerWorkgroup}u, ${SphereRasterWideChunksPerWorkgroup}u, wide);
    let n = (chunkCount - counter[o + 2u] + g - 1u) / g;
    let d = u32(uParams.uObjectIndex) * ${SphereRasterArgsWords}u;
    let narrow = select(n, 0u, wide);
    args[d] = min(narrow, 65535u);
    args[d + 1u] = (narrow + 65534u) / 65535u;
    args[d + 2u] = 1u;
    let many = select(0u, n, wide);
    args[d + 8u] = min(many, 65535u);
    args[d + 9u] = (many + 65534u) / 65535u;
    args[d + 10u] = 1u;
    // the instances of one sphere, 64 a workgroup
    let singles = (counter[o + 1u] - counter[o + 3u] + 63u) / 64u;
    args[d + 4u] = min(singles, 65535u);
    args[d + 5u] = (singles + 65534u) / 65535u;
    args[d + 6u] = 1u;
}
`;

/**
 * A workgroup per \`chunksPerWorkgroup\` chunks, their spheres dealt over its threads (a round per 64; the singles' a
 * row a thread), a thread per sphere: as the impostors size and level it (\`assign_size\`,
 * \`assign_lod_factor\`), where \`rasterizesSphere\` the pixels whose centers its front hit covers, each written with
 * \`atomicStoreMax\` of its id (\`SphereRasterTable\`: its chunk and lane) and its depth (reversed, as \`frag_depth\`).
 * Pixel clipping as the impostors' (\`sphereImpostor\`, \`dClipSphere\`): a sphere whole by its center
 * (\`dClipPrimitive\`), else per pixel the point it shows, its front, the cap of its solid interior where the front is
 * clipped (\`clipExit\`), its back where double sided, or none, the depth of that point written
 */
export const sphere_raster_wgsl = (chunksPerWorkgroup: number) => /* wgsl */ `
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

/**
 * the center and group of sphere \`i\` of instance \`slot\`: from \`tPositionGroup\`, or quantized (\`tSite\`, int16 x, y, z
 * in two words) by the object's origin and scale or (merged) the slot's segment's, the group the item (the segment's)
 */
fn readPositionGroup(i: u32, slot: u32) -> vec4<f32> {
    if (uParams.uSites != 0) {
        let w = u32(uParams.uPositionOffset) + i * 2u;
        let a = positions[w];
        let b = positions[w + 1u];
        let q = vec3<f32>(f32(bitcast<i32>(a << 16u) >> 16u), f32(bitcast<i32>(a) >> 16u), f32(bitcast<i32>(b << 16u) >> 16u));
        var originScale = vec4<f32>(uParams.uSiteOrigin, uParams.uSiteScale);
        var base = 0u;
        if (uParams.uSites == 2) {
            let segment = u32(slotSegment(slot).w);
            originScale = slotSegmentSites(segment);
            base = slotSegmentItem(segment);
        }
        return vec4<f32>(originScale.xyz + q * originScale.w, f32(i - base));
    }
    // the center's words (\`rgb32f-i32\`), the group's not read as a float (\`sphereGroup\`)
    let o = u32(uParams.uPositionOffset) + i * 4u;
    return vec4<f32>(bitcast<vec3<f32>>(vec3<u32>(positions[o], positions[o + 1u], positions[o + 2u])), 0.0);
}

/**
 * the group of sphere \`i\` of instance \`slot\` as an integer (exact past 2^24): quantized, the item (the
 * segment's); else \`tPositionGroup\`'s w
 */
fn sphereGroup(i: u32, slot: u32) -> u32 {
    if (uParams.uSites != 0) {
        var base = 0u;
        if (uParams.uSites == 2) { base = slotSegmentItem(u32(slotSegment(slot).w)); }
        return i - base;
    }
    // the fourth word an integer (\`rgb32f-i32\`)
    return positions[u32(uParams.uPositionOffset) + i * 4u + 3u];
}

/**
 * the index of quantized center \`i\`'s type in the size table (\`siteType\` sizes: its \`tSite\` element's uint16; merged,
 * from the slot's segment's table)
 */
fn siteSizeIndex(i: u32, slot: u32) -> i32 {
    let siteType = positions[u32(uParams.uPositionOffset) + i * 2u + 1u] >> 16u;
    var base = 0u;
    if (uParams.uSites == 2) { base = slotSegmentSizeBase(u32(slotSegment(slot).w)); }
    return i32(base + siteType);
}

/** the size of sphere \`i\` of group \`group\` of the chunk's instance \`e\` (at \`slot\`), as \`assign_size\` */
fn sphereSize(i: u32, group: u32, e: vec4<u32>, slot: u32) -> f32 {
    var size = uParams.uSize;
    if (uParams.uSizeIndexType == 5) {
        size = unpackSize(unpack4x8unorm(sizes[u32(uParams.uSizeOffset) + u32(siteSizeIndex(i, slot))]));
    } else if (uParams.uSizeIndexType >= 0) {
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

/** Spheres covering more pixels than this are drawn one after another, each strided over the threads */
const LargeSpherePixels = ${SphereRasterLargePixels}u;

/**
 * the workgroup's spheres of more than a pixel: their view-space center and radius, pixel bounds, id, clip flag, mode
 * and first pixel in the row of their pixels; the large ones from the end (\`wgLargeCount\`)
 */
var<workgroup> wgLargeCount: atomic<u32>;
/** the others': their count (the high 7 bits) and pixels (the low 25), one add each, so their first pixels rise with their slot */
var<workgroup> wgSharedCount: atomic<u32>;
var<workgroup> wgLargeSphere: array<vec4<f32>, WorkgroupSizeX>;
var<workgroup> wgLargeBounds: array<vec4<i32>, WorkgroupSizeX>;
var<workgroup> wgLargeInfo: array<vec4<u32>, WorkgroupSizeX>;

/** Chunks a workgroup draws (\`SphereRasterChunksPerWorkgroup\`), their spheres dealt over its threads */
const ChunksPerWorkgroup = ${chunksPerWorkgroup}u;

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(workgroup_id) wid: vec3u, @builtin(local_invocation_index) lid: u32) {
    let local = wid.x + wid.y * 65535u;
    let singles = uParams.uSingles != 0;
    let objectWord = u32(uParams.uObjectIndex) * ${SphereRasterCounterWords}u;
    let chunkEnd = min(counter[objectWord], u32(uParams.uChunkCapacity));
    // the workgroup's chunks, after those of an earlier phase (\`sphere_raster_args_wgsl\`)
    let chunkFirst = local * ChunksPerWorkgroup + counter[objectWord + 2u];
    if (!singles && chunkFirst >= chunkEnd) { return; }
    // their spheres (the same reads in every thread: uniform), a round of the workgroup for each 64
    var counts: array<u32, ChunksPerWorkgroup>;
    var total = WorkgroupSizeX;
    if (!singles) {
        total = 0u;
        for (var k = 0u; k < ChunksPerWorkgroup; k++) {
            let c = chunkFirst + k;
            counts[k] = select(0u, chunks[(u32(uParams.uChunkBase) + c) * 2u].z, c < chunkEnd);
            total += counts[k];
        }
    }
    let rounds = (total + WorkgroupSizeX - 1u) / WorkgroupSizeX;
    for (var r = 0u; r < rounds; r++) {
        // the queue read by all before the next round clears it
        if (r > 0u) { workgroupBarrier(); }
        rasterRound(r, lid, local, singles, objectWord, chunkFirst, counts);
    }
}

/**
 * The workgroup's spheres \`r * 64\` to \`r * 64 + 63\`: thread \`lid\`'s is the \`j\`-th of its chunk (singles: the row
 * \`local * 64 + lid\`)
 */
fn rasterRound(r: u32, lid: u32, local: u32, singles: bool, objectWord: u32, chunkFirst: u32, counts: array<u32, ChunksPerWorkgroup>) {
    if (lid == 0u) {
        atomicStore(&wgLargeCount, 0u);
        atomicStore(&wgSharedCount, 0u);
    }
    workgroupBarrier();

    // the thread's chunk and its place in it
    var chunkIndex = chunkFirst;
    var j = 0u;
    var inChunk = false;
    if (!singles) {
        let s = r * WorkgroupSizeX + lid;
        var before = 0u;
        for (var k = 0u; k < ChunksPerWorkgroup; k++) {
            if (!inChunk && s < before + counts[k]) {
                chunkIndex = chunkFirst + k;
                j = s - before;
                inChunk = true;
            }
            before += counts[k];
        }
    }
    // the rows of single spheres from the region's end, a row per thread (lane 0 of its id, as a chunk's first sphere)
    let single = counter[objectWord + 3u] + local * WorkgroupSizeX + lid;
    let singleValid = singles && single < min(counter[objectWord + 1u], u32(uParams.uChunkCapacity));
    let chunk = select(u32(uParams.uChunkBase) + chunkIndex, u32(uParams.uChunkBase) + u32(uParams.uChunkCapacity) - 1u - single, singles);
    let e = select(vec4<u32>(0u), chunks[chunk * 2u], inChunk || singleValid);
    // the chunk's level and stride: its \`j\`-th sphere is item \`e.y + j * stride\`
    let level = e.w & 3u;
    let i = select(e.y + j * (e.w >> 2u), e.y, singles);
    // the thread's sphere, drawn here where small, queued for the workgroup where large (no return before the barrier)
    var valid = inChunk || singleValid;
    let instance = chunks[chunk * 2u + 1u];
    var pg = vec4<f32>(0.0);
    var group = 0u;
    if (valid) {
        pg = readPositionGroup(i, e.x);
        group = sphereGroup(i, e.x);
    }
    let model = (slotTransform(e.x) * vec4<f32>(pg.xyz, 1.0)).xyz;
    // the instance's own scale (a pose's: a coarse node's particles drawn larger), its spheres' radii with it
    let rows = slotRows(e.x);
    let poseScale = length(vec3<f32>(rows[0].x, rows[1].x, rows[2].x));

    let clipped = uParams.uClipObjectCount > 0;
    var flag = 0;
    if (valid && clipped) { flag = clipFlag(group, instance); }
    if (valid && clipped && uParams.uClipPrimitive != 0 && clipTest(model, flag)) { valid = false; }
    var radius = 0.0;
    if (valid) { radius = sphereSize(i, group, instance, e.x) * poseScale * uParams.uModelScale; }
    let scale = bitcast<f32>(instance.w);
    if (valid && levelActive(level)) {
        let l = lodAt(cameraDistance(model), uParams.uLevels[level].xyz, scale);
        if (l.culled) { valid = false; }
        radius *= l.factor;
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
    let id = chunk * 64u + j + 1u;
    if (all(lo <= hi)) {
        let size = hi - lo + 1;
        let pixels = u32(size.x * size.y);
        if (pixels > LargeSpherePixels) {
            // a slot from the queue's end
            let q = WorkgroupSizeX - 1u - atomicAdd(&wgLargeCount, 1u);
            wgLargeSphere[q] = vec4<f32>(c, radius);
            wgLargeBounds[q] = vec4<i32>(lo, hi);
            wgLargeInfo[q] = vec4<u32>(id, bitcast<u32>(flag), select(0u, 1u, clipPixels), 0u);
        } else if (pixels > 1u) {
            // their pixels in a row (at most 64 * 4096: 25 bits)
            let packed = atomicAdd(&wgSharedCount, (1u << 25u) | pixels);
            let q = packed >> 25u;
            wgLargeSphere[q] = vec4<f32>(c, radius);
            wgLargeBounds[q] = vec4<i32>(lo, hi);
            wgLargeInfo[q] = vec4<u32>(id, bitcast<u32>(flag), select(0u, 1u, clipPixels), packed & 0x1ffffffu);
        } else {
            rasterPixel(lo, c, radius, id, flag, clipPixels);
        }
    }

    // the queued spheres by all threads, a pixel at a time: those in the row, strided over the threads (each pixel's
    // sphere the last whose first is at or before it), then the large
    workgroupBarrier();
    let sharedWord = atomicLoad(&wgSharedCount);
    let sharedSpheres = sharedWord >> 25u;
    let sharedPixels = sharedWord & 0x1ffffffu;
    var a = 0u;
    for (var p = lid; p < sharedPixels; p += WorkgroupSizeX) {
        // the pixels rise, so the sphere does: from the last one on
        while (a + 1u < sharedSpheres && wgLargeInfo[a + 1u].w <= p) { a++; }
        let sphere = wgLargeSphere[a];
        let bounds = wgLargeBounds[a];
        let info = wgLargeInfo[a];
        let width = u32(bounds.z - bounds.x + 1);
        let k = p - info.w;
        let px = bounds.xy + vec2<i32>(i32(k % width), i32(k / width));
        rasterPixel(px, sphere.xyz, sphere.w, info.x, bitcast<i32>(info.y), info.z != 0u);
    }
    let n = atomicLoad(&wgLargeCount);
    for (var l = 0u; l < n; l++) {
        let q = WorkgroupSizeX - 1u - l;
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

/**
 * A raster record's blocks of \`WorkgroupSizeX\` instances (\`uRecordCount\` records from \`uRecordBase\`; a record of no
 * region, its first none, none): what \`sphere-raster-regions\` dispatches a workgroup each
 */
export const sphere_raster_blocks_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3u) {
    let r = gid.x;
    if (r >= u32(uParams.uRecordCount)) { return; }
    let words = u32(uParams.uRecordWords);
    let o = (u32(uParams.uRecordBase) + r) * words;
    let count = select(records[o + 1u], 0u, records[o + words - 1u] == 0xffffffffu);
    blocks[r] = (count + 63u) / 64u;
}
`;

/** The blocks' total after their firsts (`blockFirsts`, the prefix sum's), and the regions' dispatch of them */
export const sphere_raster_blocks_args_wgsl = /* wgsl */ `
@compute @workgroup_size(1)
fn main() {
    let total = totals[0];
    blockFirsts[u32(uParams.uRecordCount)] = total;
    dispatchArgs[0] = min(total, 65535u);
    dispatchArgs[1] = (total + 65534u) / 65535u;
    dispatchArgs[2] = 1u;
}
`;
