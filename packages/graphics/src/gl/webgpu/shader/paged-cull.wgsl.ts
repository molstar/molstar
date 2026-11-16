/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PagedCullPageSize, PagedCullWorkgroupKeys } from '../../compute/paged-cull.js';
import { hiz_occlusion_wgsl } from './cell-cull.wgsl.js';
import { clip_storage_wgsl } from './chunks/clip.wgsl.js';

/** What the passes share: the store's constants */
const common_wgsl = /* wgsl */ `
const PageSize = ${PagedCullPageSize}u;
const WgKeys = ${PagedCullWorkgroupKeys}u;
const NoIndex = 0xffffffffu;
const PerParticle = 0x80000000u;

/** a particle's member, its scale, its run (NoIndex where its members are a particle); member NoIndex between runs */
struct Membership {
    member: u32,
    scale: f32,
    run: u32,
};

/** particle j of a page (its record's info, its extras' first): its run by a binary search of the page's, or its extras' */
fn membershipOf(info: vec4<u32>, extras: u32, j: u32) -> Membership {
    var m: Membership;
    m.member = NoIndex;
    m.scale = 1.0;
    m.run = NoIndex;
    if (info.y == PerParticle) {
        m.member = extrasWords[(extras + j) * 2u];
        m.scale = bitcast<f32>(extrasWords[(extras + j) * 2u + 1u]);
    } else if (info.y > 0u) {
        // the run holding j: the last whose first is at most j
        let base = u32(uParams.uRunsBase);
        var lo = info.x;
        var hi = info.x + info.y - 1u;
        while (lo < hi) {
            let mid = (lo + hi + 1u) / 2u;
            if (tables[base + mid].x <= j) { lo = mid; } else { hi = mid - 1u; }
        }
        let r = tables[base + lo];
        m.run = lo;
        // before the first run or between runs: a target the store's placement left out
        if (j >= r.x && j < r.x + r.y) {
            m.member = r.z;
            m.scale = bitcast<f32>(r.w);
            if (m.scale < 0.0) {
                m.scale = bitcast<f32>(extrasWords[(extras + j) * 2u + 1u]);
            }
        }
    }
    if (!(m.scale > 0.0)) {
        m.scale = 1.0;
    }
    return m;
}

/** the position of store index \`index\` (its page's grid \`grid\`), as stored */
fn positionOf(index: u32, grid: vec4<f32>) -> vec3<f32> {
    if (uParams.uPositions == 0) {
        // on a grid (\`GridPositions\`): 21 bits an axis in two words
        let w0 = positions[index * 2u];
        let w1 = positions[index * 2u + 1u];
        let steps = vec3<u32>(w0 & 0x1fffffu, (w0 >> 21u) | ((w1 & 0x3ffu) << 11u), (w1 >> 10u) & 0x1fffffu);
        return grid.xyz + grid.w * vec3<f32>(steps);
    }
    return bitcast<vec3<f32>>(vec3<u32>(positions[index * 3u], positions[index * 3u + 1u], positions[index * 3u + 2u]));
}
`;

/**
 * A page store's particle composed (see `paged-cull.ts`), when its page is written: its run by a binary search of the
 * page's (or its member a particle, its extras'), its position and rotation as stored (`composePages`: a page and its
 * first pose slot each), its pose (`t = p - s R c`, the unit quaternion)
 * at its slot: `t` and the scale, the quaternion snorm16, its store index, its member (`NoIndex` none: a slot past its
 * page's count, a target the placement left out); a workgroup a group of a page's particles (`WorkgroupSizeX`, the
 * cull's), its instances' centers' box into `groupBoxes` (two vec4s a group; none: min above max)
 */
export const paged_cull_compose_wgsl = /* wgsl */ `
${common_wgsl}
/** the particle of thread g composed (its pose written), its center and 1 where it has an instance, else 0 */
fn compose(g: u32) -> vec4<f32> {
    let b = g / PageSize;
    if (b >= u32(uParams.uComposeCount)) {
        return vec4<f32>(0.0);
    }
    let page = composePages[b * 2u];
    let poseFirst = composePages[b * 2u + 1u];
    let j = g % PageSize;
    let info = tables[page * 4u + 1u];
    if (poseFirst == NoIndex || j >= info.z || poseFirst + j >= u32(uParams.uPoseRoom)) {
        return vec4<f32>(0.0);
    }
    let o = u32(uParams.uPoseBase) + (poseFirst + j) * 2u;
    let index = page * PageSize + j;
    let grid = bitcast<vec4<f32>>(tables[page * 4u]);
    let m = membershipOf(info, tables[page * 4u + 2u].w, j);
    let member = m.member;
    let scale = m.scale;
    if (member == NoIndex) {
        poses[o + 1u] = vec4<u32>(0u, 0u, index, NoIndex);
        return vec4<f32>(0.0);
    }
    let p = positionOf(index, grid);
    let v = index;
    var q = vec4<f32>(0.0, 0.0, 0.0, 1.0);
    if (uParams.uRotations == 1) {
        q = bitcast<vec4<f32>>(vec4<u32>(rotations[v * 4u], rotations[v * 4u + 1u], rotations[v * 4u + 2u],
            rotations[v * 4u + 3u]));
    } else if (uParams.uRotations == 2) {
        q = vec4<f32>(unpack2x16snorm(rotations[v * 2u]), unpack2x16snorm(rotations[v * 2u + 1u]));
    }
    let n = length(q);
    q = select(vec4<f32>(0.0, 0.0, 0.0, 1.0), q / n, n > 0.0);
    // R c: c + 2 w (q × c) + 2 q × (q × c)
    let c = bitcast<vec4<f32>>(tables[u32(uParams.uMembersBase) + member]).xyz;
    let tc = 2.0 * cross(q.xyz, c);
    let rc = c + q.w * tc + cross(q.xyz, tc);
    let t = p - scale * rc;
    poses[o] = vec4<u32>(bitcast<vec3<u32>>(t), bitcast<u32>(scale));
    poses[o + 1u] = vec4<u32>(pack2x16snorm(q.xy), pack2x16snorm(q.zw), index, member);
    return vec4<f32>(p, 1.0);
}

var<workgroup> wgLo: array<vec3<f32>, WorkgroupSizeX>;
var<workgroup> wgHi: array<vec3<f32>, WorkgroupSizeX>;

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_index) lid: u32) {
    let g = gid.x + gid.y * u32(uParams.uRowThreads);
    let c = compose(g);
    // the group's box (a workgroup's particles, those of an instance): min and max by halves
    wgLo[lid] = select(vec3<f32>(3.4e38), c.xyz, c.w > 0.0);
    wgHi[lid] = select(vec3<f32>(-3.4e38), c.xyz, c.w > 0.0);
    for (var k = WorkgroupSizeX / 2u; k > 0u; k = k / 2u) {
        workgroupBarrier();
        if (lid < k) {
            wgLo[lid] = min(wgLo[lid], wgLo[lid + k]);
            wgHi[lid] = max(wgHi[lid], wgHi[lid + k]);
        }
    }
    let b = g / PageSize;
    if (lid == 0u && b < u32(uParams.uComposeCount)) {
        let group = composePages[b * 2u] * (PageSize / WorkgroupSizeX) + (g % PageSize) / WorkgroupSizeX;
        groupBoxes[group * 2u] = vec4<f32>(wgLo[0], 0.0);
        groupBoxes[group * 2u + 1u] = vec4<f32>(wgHi[0], 0.0);
    }
}
`;

/**
 * A page store's particle culled (see `paged-cull.ts`): its page from the active pages, its run (its member and
 * scale) and its position as stored (8 B a particle on a grid, not its pose's 32: a pass of the store's particles a
 * little faster), culled by the frustum at its member's reach
 * times its scale, its level by its distance from the camera's plane
 */
const instance_wgsl = /* wgsl */ `
${hiz_occlusion_wgsl}
${clip_storage_wgsl((i) => `bitcast<vec4<f32>>(tables[${i}])`, 'u32(uParams.uClipBase)', 'uParams.uClipObjectCount')}
${common_wgsl}

fn activePage(a: u32) -> u32 {
    return tables[u32(uParams.uActiveBase) + a / 4u][a % 4u];
}

fn poseFirst(page: u32) -> u32 {
    return tables[u32(uParams.uPoseFirstsBase) + page / 4u][page % 4u];
}

struct Instance {
    /** its run's index in the runs' section (the workgroup's keys), NoIndex where its members are a particle */
    run: u32,
    /** at its level (\`level\`) */
    visible: bool,
    /** the ghost's too (its member's region after the levels', \`uGhostBase\`) */
    ghost: bool,
    /** drawn by the sphere raster: its level's raster region (\`uRasterBase\`), not the level's */
    raster: bool,
    member: u32,
    level: u32,
    /** its store index (its page's first and its index): its bit in \`occludedBits\` */
    index: u32,
    /** its pose's slot: what its regions hold */
    slot: u32,
};

fn levelCount() -> u32 {
    return max(1u, u32(uParams.uLevelCount));
}

fn outside(plane: vec4<f32>, p: vec3<f32>, r: f32) -> bool {
    return dot(plane.xyz, p) + plane.w < -r;
}

/** a sphere of the object as drawn: moved by its transient transform (\`uObjectTransform\`, the identity at rest) */
fn moved(sphere: vec4<f32>) -> vec4<f32> {
    let center = uParams.uObjectTransform * vec4<f32>(sphere.xyz, 1.0);
    return vec4<f32>(center.xyz, sphere.w * uParams.uObjectScale);
}

/**
 * Whether the sphere raster draws an instance at \`p\` of reach \`radius\` (its scale \`scale\`) at \`level\` (uRasterRadius
 * above 0, as the merged culler's \`rasterized\`): in front of the near plane, its member's largest sphere at the level's
 * scale (the raster rows: a vec4 a member) at most uRasterRadius pixels across at its nearest point
 */
fn rasterized(p: vec3<f32>, radius: f32, scale: f32, member: u32, level: u32) -> bool {
    if (uParams.uRasterRadius <= 0.0) {
        return false;
    }
    let r = scale * bitcast<vec4<f32>>(tables[u32(uParams.uRasterRowsBase) + member])[min(level, 3u)];
    if (!(r > 0.0)) {
        return false;
    }
    // its spheres' centers at least at the level's minimum distance (closer ones the level culls) and in the instance
    let d = max(dot(uParams.uCameraPlane.xyz, p) + uParams.uCameraPlane.w - radius, uParams.uLevels[level].x);
    if (d - r <= uParams.uRasterNear) {
        return false;
    }
    // the tangents from the eye, as the raster's bounds; orthographic its radius (\`isRasterSize\`)
    let extent = select(sqrt(max(d * d - r * r, 0.0)), 1.0, uParams.uRasterOrtho != 0.0);
    return r * uParams.uRasterScale <= uParams.uRasterRadius * extent;
}

/** the region of a visible instance: its member's at its level, or its level's raster region's */
fn regionOf(raster: bool, member: u32, level: u32) -> u32 {
    if (raster) {
        return u32(uParams.uRasterBase) + level * u32(uParams.uMemberCount) + member;
    }
    return member * levelCount() + level;
}

/** a workgroup's page: its particles tested as they are (0), none (1, occluded whole), all in the second phase (2) */
var<workgroup> wgPageState: u32;

/**
 * The group of the workgroup's particles (all its threads, before \`instance\`; a page's group of \`WgKeys\`, its
 * centers' box \`groupBoxes\`, the compose's) tested whole, its sphere padded by the largest particle sphere
 * (\`uMaxReach\`): outside the frustum none of its particles tested (1); against the pyramid in the first phase
 * occluded none (1), the group marked (a bit a group after the particles'); in the second a marked group occluded
 * again none (1), else all its particles tested (2: none of them marked in the first). Not with the ghost (what the
 * clip objects remove drawn whatever the frustum). A group's box, not its page's: a page a tile's particles by
 * target, its box the tile's (the yeast drawn whole a frame's GPU 8.0 -> 7.4 ms)
 */
fn pageState(g: u32, lid: u32) -> u32 {
    if (lid == 0u) {
        var state = 0u;
        let a = g / PageSize;
        if (uParams.uGhost == 0 && a < u32(uParams.uActiveCount)) {
            let page = activePage(a);
            // the workgroup's group of the page (its particles' centers' box in groupBoxes)
            let group = page * (PageSize / WgKeys) + (g % PageSize) / WgKeys;
            let lo = groupBoxes[group * 2u].xyz;
            let hi = groupBoxes[group * 2u + 1u].xyz;
            let sphere = moved(vec4<f32>((lo + hi) * 0.5, length(hi - lo) * 0.5 + uParams.uMaxReach));
            let word = u32(uParams.uPageBitsBase) + group / 32u;
            let bit = 1u << (group % 32u);
            let marked = (atomicLoad(&occludedBits[word]) & bit) != 0u;
            if (all(lo <= hi) && (outside(uParams.uPlane0, sphere.xyz, sphere.w) ||
                outside(uParams.uPlane1, sphere.xyz, sphere.w) || outside(uParams.uPlane2, sphere.xyz, sphere.w) ||
                outside(uParams.uPlane3, sphere.xyz, sphere.w) || outside(uParams.uPlane4, sphere.xyz, sphere.w) ||
                outside(uParams.uPlane5, sphere.xyz, sphere.w))) {
                // outside the frustum whole
                state = 1u;
            } else if (uParams.uHiZ != 0u && (uParams.uPhase == 0 || marked)) {
                if (all(lo <= hi) && occluded(sphere)) {
                    state = 1u;
                    if (uParams.uPhase == 0) {
                        atomicOr(&occludedBits[word], bit);
                    }
                } else if (uParams.uPhase == 1) {
                    state = 2u;
                }
            }
        }
        wgPageState = state;
    }
    return workgroupUniformLoad(&wgPageState);
}

// in the first phase an instance the pyramid occludes marked (occludedBits, its store index's bit) for the second,
// which tests the marked alone (uPhase 1, the pyramid of the frame's depth so far); \`pageOf\` its page's state
// (\`pageState\`)
fn instance(g: u32, pageOf: u32) -> Instance {
    var r: Instance;
    r.run = NoIndex;
    r.visible = false;
    r.ghost = false;
    r.raster = false;
    let a = g / PageSize;
    if (a >= u32(uParams.uActiveCount) || pageOf == 1u) {
        return r;
    }
    let page = activePage(a);
    let j = g % PageSize;
    let info = tables[page * 4u + 1u];
    let first = poseFirst(page);
    if (j >= info.z || first == NoIndex || first + j >= u32(uParams.uPoseRoom)) {
        return r;
    }
    let m = membershipOf(info, tables[page * 4u + 2u].w, j);
    r.run = m.run;
    let member = m.member;
    if (member == NoIndex) {
        return r;
    }
    let scale = m.scale;
    let index = page * PageSize + j;
    let slot = first + j;
    let md = bitcast<vec4<f32>>(tables[u32(uParams.uMembersBase) + member]);
    // as drawn (the vertex stage moves the poses by the transient transform): the tests of the moved sphere
    let s = moved(vec4<f32>(positionOf(index, bitcast<vec4<f32>>(tables[page * 4u])), md.w * scale));
    let p = s.xyz;
    let radius = s.w;

    if (outside(uParams.uPlane0, p, radius) || outside(uParams.uPlane1, p, radius) ||
        outside(uParams.uPlane2, p, radius) || outside(uParams.uPlane3, p, radius) ||
        outside(uParams.uPlane4, p, radius) || outside(uParams.uPlane5, p, radius)) {
        return r;
    }
    // the clip objects: a sphere they clip whole not drawn, with instance clipping one clipped by its center (as the
    // vertex stage clips instances)
    var clipClass = 0u;
    if (uParams.uClipObjectCount > 0) {
        clipClass = clipSphereClass(p, radius, 0);
    }
    let centerClipped = uParams.uClipInstance != 0 && clipTest(p, 0);
    // the ghost's (uGhost): reached by the clip objects (with none classified all), at the coarsest level, not by the
    // pyramid; in the first phase alone
    r.ghost = uParams.uGhost != 0 && uParams.uPhase == 0 && (uParams.uClipObjectCount == 0 || clipClass != 0u) &&
        (uParams.uClipInstance == 0 || centerClipped);
    var visible = clipClass != 2u && !centerClipped;
    // its level: the first whose distance range holds its center's (none: drawn at none)
    var level = 0u;
    if (visible && uParams.uLevelCount > 0) {
        let d = dot(uParams.uCameraPlane.xyz, p) + uParams.uCameraPlane.w;
        level = NoIndex;
        for (var k = 0u; k < u32(uParams.uLevelCount); k++) {
            let range = uParams.uLevels[k];
            if (d >= range.x && d < range.y) {
                level = k;
                break;
            }
        }
        visible = level != NoIndex;
    }
    // occluded: in the first phase by the last frame's pyramid (marked for the second), in the second the marked alone
    // by the pyramid of the frame's depth so far
    let bit = 1u << (index % 32u);
    if (visible && uParams.uPhase == 1 && pageOf != 2u) {
        visible = (atomicLoad(&occludedBits[index / 32u]) & bit) != 0u;
    }
    if (visible && uParams.uHiZ != 0u && occluded(vec4<f32>(p, radius))) {
        if (uParams.uPhase == 0) {
            atomicOr(&occludedBits[index / 32u], bit);
        }
        visible = false;
    }
    if (!visible && !r.ghost) {
        return r;
    }
    r.visible = visible;
    r.raster = visible && rasterized(p, radius, scale, member, level);
    r.member = member;
    r.level = level;
    r.index = index;
    r.slot = slot;
    return r;
}

/**
 * a workgroup's key of an instance: its run's index from the first thread's times the levels and its level, twice
 * (the raster's odd), its particles consecutive in one page (a few runs); none (WgKeys) where its members are a
 * particle or past the keys
 */
fn workgroupKey(r: Instance, first: u32) -> u32 {
    if (!r.visible || r.run == NoIndex || first == NoIndex || r.run < first) {
        return WgKeys;
    }
    let key = ((r.run - first) * levelCount() + r.level) * 2u + select(0u, 1u, r.raster);
    return select(WgKeys, key, key < WgKeys);
}

/** the region of a workgroup's key: its run's member and its level, the raster's or not */
fn keyRegion(first: u32, key: u32) -> u32 {
    let k = key / 2u;
    let run = tables[u32(uParams.uRunsBase) + first + k / levelCount()];
    return regionOf(key % 2u == 1u, run.z, k % levelCount());
}
`;

/**
 * Per particle, the visible ones' pose slots written into their regions (a member and level,
 * \`member * levels + level\`; the raster's, the ghost's), a workgroup's by key first (its runs' and levels': one global
 * atomic a key). A level's and its raster's regions share a
 * member's slots of the level (\`regions\`, laid out on the CPU: a member's room its instances; a particle at one
 * level, drawn by one), the level's from its first up, the raster's from its end (\`regions\` of a raster region) down;
 * in the second phase after the first's (\`regions\` from \`uRegionCount\`: the first phase's counts, copied there)
 */
export const paged_cull_wgsl = /* wgsl */ `
${instance_wgsl}
var<workgroup> wgFirst: u32;
var<workgroup> wgCounts: array<atomic<u32>, WgKeys>;
var<workgroup> wgSlots: array<u32, WgKeys>;

fn isRaster(region: u32) -> bool {
    return region >= u32(uParams.uRasterBase) && region < u32(uParams.uGhostBase);
}

/**
 * the first of \`c\` slots of \`region\` taken at its cursor \`at\`: a raster region's down from its end, the others' up;
 * in the second phase past the first's (\`regions\` from \`uRegionCount\`)
 */
fn regionSlot(region: u32, at: u32, c: u32) -> u32 {
    let past = select(0u, regions[u32(uParams.uRegionCount) + region], uParams.uPhase == 1) + at;
    if (isRaster(region)) {
        return regions[region] - past - c;
    }
    return regions[region] + past;
}

/** pose slot v at index slot k (none past the indices' room, the device's buffer limit) */
fn writeIndex(k: u32, v: u32) {
    if (k < u32(uParams.uIndexRoom)) {
        indices[k] = v;
    }
}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_index) lid: u32) {
    let g = gid.x + gid.y * u32(uParams.uRowThreads);
    let r = instance(g, pageState(g, lid));
    atomicStore(&wgCounts[lid], 0u);
    if (lid == 0u) {
        wgFirst = r.run;
    }
    workgroupBarrier();
    let first = workgroupUniformLoad(&wgFirst);
    let key = workgroupKey(r, first);
    var local = 0u;
    if (r.visible && key < WgKeys) {
        local = atomicAdd(&wgCounts[key], 1u);
    }
    workgroupBarrier();
    let c = atomicLoad(&wgCounts[lid]);
    if (c > 0u) {
        let region = keyRegion(first, lid);
        wgSlots[lid] = regionSlot(region, atomicAdd(&counts[region], c), c);
    }
    workgroupBarrier();
    if (r.ghost) {
        let region = u32(uParams.uGhostBase) + r.member;
        writeIndex(regionSlot(region, atomicAdd(&counts[region], 1u), 1u), r.slot);
    }
    if (!r.visible) {
        return;
    }
    if (key < WgKeys) {
        writeIndex(wgSlots[key] + local, r.slot);
    } else {
        let region = regionOf(r.raster, r.member, r.level);
        writeIndex(regionSlot(region, atomicAdd(&counts[region], 1u), 1u), r.slot);
    }
}
`;

/**
 * Per draw record (a level's, a member's: level-major, `j * segmentCount + s`), its instance count and first instance
 * from its region's count and first (`counts[s * levels + j]`, `bases`: a raster region's its end, its slots down from
 * it; the second phase's past the first's, `prior`), the first its last word; a member past the regions none. The
 * raster's records after the levels' (`uRasterLevels` blocks, their regions level-major after the levels'; a member
 * past the regions' first past all, as the raster finds a slot's record by the first), the ghost's after them
 */
export const paged_cull_args_wgsl = /* wgsl */ `
const NoIndex = 0xffffffffu;
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let r = gid.x;
    if (r >= u32(uParams.uRecordCount)) {
        return;
    }
    let segments = u32(uParams.uSegmentCount);
    let levels = u32(uParams.uLevels);
    let s = r % segments;
    let j = r / segments;
    let words = u32(uParams.uRecordWords);
    let members = u32(uParams.uMemberCount);
    let rasterLevels = u32(uParams.uRasterLevels);
    let raster = j >= levels && j < levels + rasterLevels;
    var count = 0u;
    var first = select(0u, NoIndex, raster);
    if (s < members) {
        var region = s * levels + j;
        if (raster) {
            region = members * levels + (j - levels) * members + s;
        } else if (j >= levels) {
            region = members * (levels + rasterLevels) + s;
        }
        let past = select(0u, prior[region], uParams.uSecond != 0);
        count = counts[region];
        if (raster) {
            first = bases[region] - past - count;
        } else {
            first = bases[region] + past;
        }
        // cut at the indices' room: none written past it
        let room = u32(uParams.uIndexRoom);
        count = min(count, room - min(room, first));
    }
    args[r * words + 1u] = count;
    args[r * words + words - 1u] = first;
}
`;
