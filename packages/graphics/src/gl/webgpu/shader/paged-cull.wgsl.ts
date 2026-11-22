/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import {
  PagedCullGroupArgs,
  PagedCullGroupsWorkgroup,
  PagedCullMaxRow,
  PagedCullOrderBuckets,
  PagedCullPageSize,
  PagedCullWorkgroupKeys,
  PagedEntryLevelShift,
  PagedEntryLevelsShift,
  PagedEntryPageMask,
} from '../../compute/paged-cull.js';
import { PagedCopyShift, PagedSegmentMask } from '../../culled-data.js';
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

/** word k of the extras (two a particle: its member, its scale), in the tables' section from \`uExtrasBase\` */
fn extraWord(k: u32) -> u32 {
    return tables[u32(uParams.uExtrasBase) + k / 4u][k % 4u];
}

/** particle j of a page (its record's info, its extras' first): its run by a binary search of the page's, or its extras' */
fn membershipOf(info: vec4<u32>, extras: u32, j: u32) -> Membership {
    var m: Membership;
    m.member = NoIndex;
    m.scale = 1.0;
    m.run = NoIndex;
    if (info.y == PerParticle) {
        m.member = extraWord((extras + j) * 2u);
        m.scale = bitcast<f32>(extraWord((extras + j) * 2u + 1u));
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
                m.scale = bitcast<f32>(extraWord((extras + j) * 2u + 1u));
            }
        }
    }
    if (!(m.scale > 0.0)) {
        m.scale = 1.0;
    }
    return m;
}

/**
 * a particle's words of position (a grid's two, float32's three) and of rotation (none, float32's four, snorm16's
 * two)
 */
fn positionWords() -> u32 {
    return select(3u, 2u, uParams.uPositions == 0);
}
fn rotationWords() -> u32 {
    return select(select(0u, 4u, uParams.uRotations == 1), 2u, uParams.uRotations == 2);
}

/**
 * store index \`index\`'s first word of its page's block in \`positions\` (a page its positions, then its rotations):
 * the block a written page has (the tables' section from \`uValueBlocksBase\`, a word a page), not the page's own place
 */
fn pageBase(index: u32) -> u32 {
    let page = index / PageSize;
    let block = tables[u32(uParams.uValueBlocksBase) + page / 4u][page % 4u];
    return block * PageSize * (positionWords() + rotationWords());
}

/** the position of store index \`index\` (its page's grid \`grid\`), as stored */
fn positionOf(index: u32, grid: vec4<f32>) -> vec3<f32> {
    let j = index % PageSize;
    if (uParams.uPositions == 0) {
        // on a grid (\`GridPositions\`): 21 bits an axis in two words
        let o = pageBase(index) + j * 2u;
        let w0 = positions[o];
        let w1 = positions[o + 1u];
        let steps = vec3<u32>(w0 & 0x1fffffu, (w0 >> 21u) | ((w1 & 0x3ffu) << 11u), (w1 >> 10u) & 0x1fffffu);
        return grid.xyz + grid.w * vec3<f32>(steps);
    }
    let o = pageBase(index) + j * 3u;
    return bitcast<vec3<f32>>(vec3<u32>(positions[o], positions[o + 1u], positions[o + 2u]));
}

/** the unit rotation of store index \`index\`, as stored (the identity where none) */
fn rotationOf(index: u32) -> vec4<f32> {
    let o = pageBase(index) + PageSize * positionWords() + (index % PageSize) * rotationWords();
    var q = vec4<f32>(0.0, 0.0, 0.0, 1.0);
    if (uParams.uRotations == 1) {
        q = bitcast<vec4<f32>>(vec4<u32>(positions[o], positions[o + 1u], positions[o + 2u], positions[o + 3u]));
    } else if (uParams.uRotations == 2) {
        q = vec4<f32>(unpack2x16snorm(positions[o]), unpack2x16snorm(positions[o + 1u]));
    }
    let n = length(q);
    return select(vec4<f32>(0.0, 0.0, 0.0, 1.0), q / n, n > 0.0);
}
`;

/**
 * A page store's particle composed (see `paged-cull.ts`), when its page is written: its run by a binary search of the
 * page's (or its member a particle, its extras'), its position and rotation as stored (`composePages`: a page and its
 * first pose slot each, \`NoIndex\` none: a template's page, its copies' poses composed by the cull from its values),
 * its pose (`t = p - s R c`, the unit quaternion) at its slot: `t` and the scale, the quaternion snorm16, its store
 * index, its member (`NoIndex` none: a slot past its page's count, a target the placement left out); a workgroup a
 * group of a page's particles (`WorkgroupSizeX`, the cull's), its instances' centers' box into `groupBoxes` (two vec4s
 * a group; none: min above max), a page's with pose slots or not
 */
export const paged_cull_compose_wgsl = /* wgsl */ `
${common_wgsl}
/**
 * the particle of thread g composed (its pose written where it has a slot), its center and its scale where it has an
 * instance, else 0
 */
fn compose(g: u32) -> vec4<f32> {
    let b = g / PageSize;
    if (b >= u32(uParams.uComposeCount)) {
        return vec4<f32>(0.0);
    }
    let page = composePages[b * 2u];
    let poseFirst = composePages[b * 2u + 1u];
    let j = g % PageSize;
    let info = tables[page * 4u + 1u];
    if (j >= info.z) {
        return vec4<f32>(0.0);
    }
    let posed = poseFirst != NoIndex && poseFirst + j < u32(uParams.uPoseRoom);
    let o = u32(uParams.uPoseBase) + (poseFirst + j) * 2u;
    let index = page * PageSize + j;
    let grid = bitcast<vec4<f32>>(tables[page * 4u]);
    let m = membershipOf(info, tables[page * 4u + 2u].w, j);
    let member = m.member;
    let scale = m.scale;
    if (member == NoIndex) {
        if (posed) {
            poses[o + 1u] = vec4<u32>(0u, 0u, index, NoIndex);
        }
        return vec4<f32>(0.0);
    }
    let p = positionOf(index, grid);
    if (!posed) {
        return vec4<f32>(p, scale);
    }
    let q = rotationOf(index);
    // R c: c + 2 w (q × c) + 2 q × (q × c)
    let c = bitcast<vec4<f32>>(tables[u32(uParams.uMembersBase) + member]).xyz;
    let tc = 2.0 * cross(q.xyz, c);
    let rc = c + q.w * tc + cross(q.xyz, tc);
    let t = p - scale * rc;
    poses[o] = vec4<u32>(bitcast<vec3<u32>>(t), bitcast<u32>(scale));
    poses[o + 1u] = vec4<u32>(pack2x16snorm(q.xy), pack2x16snorm(q.zw), index, member);
    return vec4<f32>(p, scale);
}

var<workgroup> wgLo: array<vec3<f32>, WorkgroupSizeX>;
var<workgroup> wgHi: array<vec3<f32>, WorkgroupSizeX>;
var<workgroup> wgScale: array<f32, WorkgroupSizeX>;

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_index) lid: u32) {
    let g = gid.x + gid.y * u32(uParams.uRowThreads);
    let c = compose(g);
    // the group's box (a workgroup's particles, those of an instance): min and max by halves
    wgLo[lid] = select(vec3<f32>(3.4e38), c.xyz, c.w > 0.0);
    wgHi[lid] = select(vec3<f32>(-3.4e38), c.xyz, c.w > 0.0);
    wgScale[lid] = max(c.w, 0.0);
    for (var k = WorkgroupSizeX / 2u; k > 0u; k = k / 2u) {
        workgroupBarrier();
        if (lid < k) {
            wgLo[lid] = min(wgLo[lid], wgLo[lid + k]);
            wgHi[lid] = max(wgHi[lid], wgHi[lid + k]);
            wgScale[lid] = max(wgScale[lid], wgScale[lid + k]);
        }
    }
    let b = g / PageSize;
    if (lid == 0u && b < u32(uParams.uComposeCount)) {
        let group = composePages[b * 2u] * (PageSize / WorkgroupSizeX) + (g % PageSize) / WorkgroupSizeX;
        // its largest scale in w: its sphere padded by the largest reach times it (\`groupState\`)
        groupBoxes[group * 2u] = vec4<f32>(wgLo[0], wgScale[0]);
        groupBoxes[group * 2u + 1u] = vec4<f32>(wgHi[0], 0.0);
    }
}
`;

/**
 * The active entries and the copies (the select's and the cull's): an entry's page word and copy, a copy's pose,
 * sphere and level, the selected entries (\`uSelectedBase\`: a header vec4, their count in x, then an entry a word)
 */
const entries_wgsl = /* wgsl */ `
/** active entry a's page word (two words an entry: its page and its levels, its copy) */
fn activeWord(a: u32) -> u32 {
    return tables[u32(uParams.uActiveBase) + a / 2u][(a % 2u) * 2u];
}

/** active entry a's page */
fn activePage(a: u32) -> u32 {
    return activeWord(a) & ${PagedEntryPageMask}u;
}

/** active entry a's copy: NoIndex a page drawn as it is, else a template's page placed by copy \`activeCopy(a)\` */
fn activeCopy(a: u32) -> u32 {
    return tables[u32(uParams.uActiveBase) + a / 2u][(a % 2u) * 2u + 1u];
}

/** copy c's pose: \`t\` (its scale 1 in w), its unit quaternion (f32 bits); its center and its template's spacing */
fn copyT(c: u32) -> vec3<f32> {
    return bitcast<vec4<f32>>(tables[u32(uParams.uCopiesBase) + 4u * c]).xyz;
}
fn copyQ(c: u32) -> vec4<f32> {
    return bitcast<vec4<f32>>(tables[u32(uParams.uCopiesBase) + 4u * c + 1u]);
}
fn copySphere(c: u32) -> vec4<f32> {
    return bitcast<vec4<f32>>(tables[u32(uParams.uCopiesBase) + 4u * c + 2u]);
}
/** copy c's template's radius about its center and its block's largest scale (below 0: none known) */
fn copyBounds(c: u32) -> vec2<f32> {
    return bitcast<vec4<f32>>(tables[u32(uParams.uCopiesBase) + 4u * c + 3u]).xy;
}

/** copy c's finest level drawn (its template's leaf pages not resident: 1, else 0) */
fn copyLeastLevel(c: u32) -> u32 {
    return u32(bitcast<vec4<f32>>(tables[u32(uParams.uCopiesBase) + 4u * c + 3u]).z);
}

/**
 * copy c's group: its template's copy of the same top particle in every block (its leaf block's and its levels'),
 * what a copy lacking a leaf page holds coarse in all of them (\`lacking\`, \`paged-cull-wants\`)
 */
fn copyGroup(c: u32) -> u32 {
    return u32(bitcast<vec4<f32>>(tables[u32(uParams.uCopiesBase) + 4u * c + 3u]).w);
}

/** page p's value block (\`NoIndex\` none: not written, a leaf page not resident) */
fn valueBlockOf(p: u32) -> u32 {
    return tables[u32(uParams.uValueBlocksBase) + p / 4u][p % 4u];
}

/** v turned by the unit quaternion q: v + 2 w (q × v) + 2 q × (q × v) */
fn turnByQuat(q: vec4<f32>, v: vec3<f32>) -> vec3<f32> {
    let t = 2.0 * cross(q.xyz, v);
    return v + q.w * t + cross(q.xyz, t);
}

/** the point v of a template's frame where copy c places it */
fn copyPoint(c: u32, v: vec3<f32>) -> vec3<f32> {
    return copyT(c) + turnByQuat(copyQ(c), v);
}

fn inFrustum(s: vec4<f32>) -> bool {
    return !(outside(uParams.uPlane0, s.xyz, s.w) || outside(uParams.uPlane1, s.xyz, s.w) ||
        outside(uParams.uPlane2, s.xyz, s.w) || outside(uParams.uPlane3, s.xyz, s.w) ||
        outside(uParams.uPlane4, s.xyz, s.w) || outside(uParams.uPlane5, s.xyz, s.w));
}

/**
 * the coarse level copy c draws its template at (of \`levels\`, 0 the template whole): its template's spacing at its
 * center's distance in pixels (\`uCoarseScale\`), a level each halving below \`uCoarsePixels\`; not finer than its
 * least (\`copyLeastLevel\`: a copy whose leaf pages are not resident drawn at its first coarse level meanwhile)
 */
fn copyLevel(c: u32, levels: u32) -> u32 {
    if (levels == 0u) {
        return 0u;
    }
    let least = min(levels, copyLeastLevel(c));
    if (uParams.uCoarsePixels <= 0.0) {
        return least;
    }
    let s = copySphere(c);
    let center = moved(vec4<f32>(s.xyz, s.w));
    let d = select(max(dot(uParams.uCameraPlane.xyz, center.xyz) + uParams.uCameraPlane.w, 1e-3), 1.0,
        uParams.uCoarseOrtho != 0);
    let px = center.w * uParams.uCoarseScale / d;
    if (!(px < uParams.uCoarsePixels)) {
        return least;
    }
    return max(least, min(levels, u32(floor(log2(uParams.uCoarsePixels / max(px, 1e-9))))));
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
 * the selected entries the first phase's select chose (the header's y): those after it the second phase's select
 * appended (copies the first left out, occluded whole by the last frame's pyramid, shown by the frame's)
 */
fn firstSelected() -> u32 {
    return tables[u32(uParams.uSelectedBase)].y;
}

/** the k-th selected entry (an active entry's index), 0xffffffff past their count */
fn selectedEntry(k: u32) -> u32 {
    let base = u32(uParams.uSelectedBase);
    if (k >= tables[base].x) {
        return 0xffffffffu;
    }
    let i = 4u + k;
    return tables[base + i / 4u][i % 4u];
}
`;

/**
 * the order bucket of an entry of copy c (\`uOrderNear\`, \`uOrderScale\`: the copies' depth range in the view over the
 * buckets, near first); a page drawn as it is the first's
 */
const order_wgsl = /* wgsl */ `
const OrderBuckets = ${PagedCullOrderBuckets}u;
fn bucketOf(c: u32) -> u32 {
    if (c == NoIndex) {
        return 0u;
    }
    let s = moved(vec4<f32>(copySphere(c).xyz, 0.0));
    let d = dot(uParams.uCameraPlane.xyz, s.xyz) + uParams.uCameraPlane.w;
    return u32(clamp((d - uParams.uOrderNear) * uParams.uOrderScale, 0.0, f32(OrderBuckets - 1u)));
}
`;

/** the phases proper (not the fallback's pass, \`uFallback\`: it marks nothing, has no ghost, tests every particle) */
const phases_wgsl = /* wgsl */ `
fn firstPhase() -> bool {
    return uParams.uPhase == 0 && uParams.uFallback == 0;
}
fn secondPhase() -> bool {
    return uParams.uPhase == 1 && uParams.uFallback == 0;
}
`;

/**
 * Group \`local\` of selected entry \`k\` (a page's group of \`WgKeys\` particles, its centers' box \`groupBoxes\`, the
 * compose's) tested whole, its sphere padded by the largest particle sphere (\`uMaxReach\`): outside the frustum none
 * of its particles tested (1); against the pyramid in the first phase occluded none (1), the group marked (a bit a
 * group after the particles'); in the second a marked group occluded again none (1), else all its particles tested (2:
 * none of them marked in the first), one not occluded whole none where none of its particles is marked (1). Not with
 * the ghost (what the clip objects remove drawn whatever the frustum). A group's box, not its page's: a page a tile's
 * particles by target, its box the tile's (the yeast drawn whole a frame's GPU 8.0 -> 7.4 ms)
 */
const group_state_wgsl = /* wgsl */ `
fn groupState(k: u32, local: u32) -> u32 {
    var state = 0u;
    // the group's selected entry (its copy at its page's level: \`paged-cull-select\`)
    let a = selectedEntry(k);
    if (uParams.uGhost == 0 && a != NoIndex) {
        let page = activePage(a);
        let copy = activeCopy(a);
        // the group of the page (its particles' centers' box in groupBoxes: a template's in its frame, its sphere
        // placed by the copy), its bit by its entry's (a template's page drawn once a copy)
        let group = page * (PageSize / WgKeys) + local;
        let lo4 = groupBoxes[group * 2u];
        let lo = lo4.xyz;
        let hi = groupBoxes[group * 2u + 1u].xyz;
        var center = (lo + hi) * 0.5;
        if (copy != NoIndex) {
            center = copyPoint(copy, center);
        }
        // padded by its particles' largest sphere: the largest reach times the group's largest scale
        let sphere = moved(vec4<f32>(center, length(hi - lo) * 0.5 + uParams.uMaxReach * lo4.w));
        // an entry past the bits' room (\`uBitsRoom\`) not tested against the pyramid: drawn in the first phase
        let tested = k < u32(uParams.uBitsRoom);
        let entryGroup = k * (PageSize / WgKeys) + local;
        let word = u32(uParams.uPageBitsBase) + entryGroup / 32u;
        let bit = 1u << (entryGroup % 32u);
        var marked = false;
        if (tested) {
            marked = (atomicLoad(&occludedBits[word]) & bit) != 0u;
        }
        if (all(lo <= hi) && (outside(uParams.uPlane0, sphere.xyz, sphere.w) ||
            outside(uParams.uPlane1, sphere.xyz, sphere.w) || outside(uParams.uPlane2, sphere.xyz, sphere.w) ||
            outside(uParams.uPlane3, sphere.xyz, sphere.w) || outside(uParams.uPlane4, sphere.xyz, sphere.w) ||
            outside(uParams.uPlane5, sphere.xyz, sphere.w))) {
            // outside the frustum whole
            state = 1u;
        } else if (secondPhase() && k >= firstSelected()) {
            // a copy the first phase left out (occluded whole), shown by the frame's pyramid: its groups tested
            // whole against it, the rest all their particles (none marked: the first culled none of them)
            if (uParams.uHiZ != 0u && all(lo <= hi) && occluded(sphere)) {
                state = 1u;
            } else {
                state = 2u;
            }
        } else if (!tested) {
            if (secondPhase()) {
                state = 1u;
            }
        } else if (uParams.uHiZ != 0u && (!secondPhase() || marked)) {
            if (all(lo <= hi) && occluded(sphere)) {
                state = 1u;
                if (firstPhase()) {
                    atomicOr(&occludedBits[word], bit);
                }
            } else if (secondPhase()) {
                state = 2u;
            }
        } else if (secondPhase()) {
            // the group not occluded whole: none of its particles marked none tested again (their words, before a
            // particle's values are read: the second phase a pass of the first's cost else)
            state = 1u;
            let first = (k * PageSize + local * WgKeys) / 32u;
            for (var w = 0u; w < WgKeys / 32u; w++) {
                if (atomicLoad(&occludedBits[first + w]) != 0u) {
                    state = 0u;
                    break;
                }
            }
        }
    }
    return state;
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
${entries_wgsl}
${phases_wgsl}

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
    /** its pose's slot: what its regions hold (a copy's: its frame slot, its pose composed there by \`composeCopy\`) */
    slot: u32,
    /** its entry's copy (\`activeCopy\`), NoIndex a page's drawn as it is */
    copy: u32,
    /** a copy's: its template particle's position (in its template's frame) and its scale */
    local: vec3<f32>,
    scale: f32,
};

fn levelCount() -> u32 {
    return max(1u, u32(uParams.uLevelCount));
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

// in the first phase an instance the pyramid occludes marked (occludedBits, its store index's bit) for the second,
// which tests the marked alone (uPhase 1, the pyramid of the frame's depth so far); \`pageOf\` its group's state
// (\`groupState\`, 0 or 2: a group culled whole none of its workgroup)
fn instance(g: u32, pageOf: u32) -> Instance {
    var r: Instance;
    r.run = NoIndex;
    r.visible = false;
    r.ghost = false;
    r.raster = false;
    r.copy = NoIndex;
    let k = g / PageSize;
    let a = selectedEntry(k);
    if (a == NoIndex || pageOf == 1u) {
        return r;
    }
    let page = activePage(a);
    let copy = activeCopy(a);
    let j = g % PageSize;
    let info = tables[page * 4u + 1u];
    let first = poseFirst(page);
    // a page drawn as it is by its pose slots; a template's (a copy's) none: composed from its values (\`composeCopy\`)
    if (j >= info.z || (copy == NoIndex && (first == NoIndex || first + j >= u32(uParams.uPoseRoom)))) {
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
    // a template's particle where its copy places it
    var stored = positionOf(index, bitcast<vec4<f32>>(tables[page * 4u]));
    r.local = stored;
    r.scale = scale;
    if (copy != NoIndex) {
        stored = copyPoint(copy, stored);
    }
    // as drawn (the vertex stage moves the poses by the transient transform): the tests of the moved sphere
    let s = moved(vec4<f32>(stored, md.w * scale));
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
    r.ghost = uParams.uGhost != 0 && firstPhase() && (uParams.uClipObjectCount == 0 || clipClass != 0u) &&
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
    // by the pyramid of the frame's depth so far; its bit its entry's particle's (a template's drawn once a copy)
    let b = k * PageSize + j;
    let bit = 1u << (b % 32u);
    if (visible && secondPhase() && pageOf != 2u) {
        visible = k < u32(uParams.uBitsRoom) && (atomicLoad(&occludedBits[b / 32u]) & bit) != 0u;
    }
    if (visible && uParams.uHiZ != 0u && k < u32(uParams.uBitsRoom) && occluded(vec4<f32>(p, radius))) {
        if (firstPhase()) {
            atomicOr(&occludedBits[b / 32u], bit);
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
    r.copy = copy;
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
 * atomic a key). A region's slots its room (\`regions\`, laid out on the CPU by what culls needed, their rooms from
 * twice \`uRegionCount\`; a slot past it not written, its count still taken: read back, the room grown), a level's from
 * its first up, a raster region's from its end (\`regions\` of a raster region) down; in the second phase after the
 * first's (\`regions\` from \`uRegionCount\`: the first phase's counts, copied there)
 */
export const paged_cull_wgsl = /* wgsl */ `
${instance_wgsl}
const GroupsPerPage = PageSize / WgKeys;
var<workgroup> wgFirst: u32;
var<workgroup> wgCounts: array<atomic<u32>, WgKeys>;
var<workgroup> wgSlots: array<u32, WgKeys>;
var<workgroup> wgRegions: array<u32, WgKeys>;
var<workgroup> wgFrameCount: atomic<u32>;
var<workgroup> wgFrameBase: u32;

/**
 * The pose of a copy's instance \`r\` at frame slot \`frame\`: its template particle's composed from its values as
 * stored (\`t_x = p - s_x R_x c\`, its member's center \`c\`, as the compose kernel's), placed by its copy (\`t_c + R_c
 * t_x\`, the scale \`s_x\`, \`q_c q_x\`; its store index and member the template particle's). Why from its values: a
 * template's particles' poses kept a slot each grew the poses' buffer by its size as it was made resident (the yeast
 * cell's 31.5 M: 1 GB, a second or more a frame of a zoom)
 */
fn composeCopy(r: Instance, frame: u32) {
    let qx = rotationOf(r.index);
    let c = bitcast<vec4<f32>>(tables[u32(uParams.uMembersBase) + r.member]).xyz;
    let tx = r.local - r.scale * turnByQuat(qx, c);
    let copy = r.copy;
    let qc = copyQ(copy);
    let t = copyT(copy) + turnByQuat(qc, tx);
    let q = vec4<f32>(qc.w * qx.xyz + qx.w * qc.xyz + cross(qc.xyz, qx.xyz), qc.w * qx.w - dot(qc.xyz, qx.xyz));
    let f = u32(uParams.uPoseBase) + (u32(uParams.uFrameBase) + frame) * 2u;
    poses[f] = vec4<u32>(bitcast<vec3<u32>>(t), bitcast<u32>(r.scale));
    // the template particle's store index and member, the copy plus one above the member (\`PagedCopyShift\`: a pick's
    // copy; none past 2^16 - 2 copies)
    let tag = select(0u, (copy + 1u) << ${PagedCopyShift}u, copy + 1u <= ${PagedSegmentMask}u);
    let word = (r.member & ${PagedSegmentMask}u) | tag;
    poses[f + 1u] = vec4<u32>(pack2x16snorm(q.xy), pack2x16snorm(q.zw), r.index, word);
}

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

/**
 * pose slot v at index slot k of \`region\`: none past its room (\`regions\` from twice \`uRegionCount\`: a raster
 * region's below its end) nor past the indices' (the device's buffer limit)
 */
fn writeIndex(region: u32, k: u32, v: u32) {
    let room = regions[2u * u32(uParams.uRegionCount) + region];
    let lo = select(regions[region], regions[region] - room, isRaster(region));
    if (k >= lo && k - lo < room && k < u32(uParams.uIndexRoom)) {
        indices[k] = v;
    }
}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_index) lid: u32) {
    // a workgroup a group the group pass kept (\`paged-cull-groups\`: its entry and group, all its particles tested in
    // the bit above; their count the selected header's z), not a workgroup a group of every selected page. Why: a
    // dense view's groups mostly culled whole, their workgroups' warps stalled at the barriers below (the isg 16x16x16
    // near: a third of the issue stalls, 87 % of warp slots compute)
    let w = wid.x + wid.y * ${PagedCullMaxRow}u;
    if (w >= tables[u32(uParams.uSelectedBase)].z) {
        return;
    }
    let item = regions[u32(uParams.uSurvivorsBase) + 1u + ${PagedCullOrderBuckets}u + w];
    let group = item & 0x7fffffffu;
    let pageOf = select(0u, 2u, (item >> 31u) != 0u);
    let g = (group / GroupsPerPage) * PageSize + (group % GroupsPerPage) * WgKeys + lid;
    var r = instance(g, pageOf);
    // a copy's visible instance (and its ghost) a frame slot, its pose composed there: the workgroup's taken at once
    // (one global atomic a workgroup: its particles one entry's), those past the frame's room not drawn, their copy's
    // group dropped (\`uDroppedBase\`: the fallback's, drawn a level coarser into its slots after the phases')
    let fallback = uParams.uFallback != 0;
    let frameRoom = u32(select(uParams.uFrameRoom, uParams.uFallbackRoom, fallback));
    if (lid == 0u) {
        atomicStore(&wgFrameCount, 0u);
    }
    workgroupBarrier();
    let framed = r.copy != NoIndex && (r.visible || r.ghost);
    var frameLocal = 0u;
    if (framed) {
        frameLocal = atomicAdd(&wgFrameCount, 1u);
    }
    workgroupBarrier();
    if (lid == 0u) {
        let n = atomicLoad(&wgFrameCount);
        let cursor = u32(uParams.uFrameCursor) + select(0u, 1u, fallback);
        wgFrameBase = select(0u, atomicAdd(&counts[cursor], n), n > 0u);
        if (!fallback && n > 0u && wgFrameBase + n > frameRoom) {
            atomicStore(&occludedBits[u32(uParams.uDroppedBase) + copyGroup(activeCopy(selectedEntry(g / PageSize)))], 1u);
        }
    }
    let frameBase = workgroupUniformLoad(&wgFrameBase);
    if (framed) {
        let frame = frameBase + frameLocal;
        if (frame < frameRoom) {
            let at = frame + select(0u, u32(uParams.uFrameRoom), fallback);
            composeCopy(r, at);
            r.slot = u32(uParams.uFrameBase) + at;
        } else {
            r.visible = false;
            r.ghost = false;
            r.raster = false;
        }
    }
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
        wgRegions[lid] = region;
    }
    workgroupBarrier();
    if (r.ghost) {
        let region = u32(uParams.uGhostBase) + r.member;
        writeIndex(region, regionSlot(region, atomicAdd(&counts[region], 1u), 1u), r.slot);
    }
    if (!r.visible) {
        return;
    }
    if (key < WgKeys) {
        writeIndex(wgRegions[key], wgSlots[key] + local, r.slot);
    } else {
        let region = regionOf(r.raster, r.member, r.level);
        writeIndex(region, regionSlot(region, atomicAdd(&counts[region], 1u), 1u), r.slot);
    }
}
`;

/**
 * The entries the cull runs over (\`paged-cull-select\`, a thread an active entry, see \`paged-cull.ts\`): a page drawn
 * as it is; a template's page at its copy's level, the copy's sphere (its template's about its center, padded by the
 * largest reach) in the frustum (not with the ghost: what clip objects remove drawn whatever the frustum)
 */
export const paged_cull_select_wgsl = /* wgsl */ `
const NoIndex = 0xffffffffu;
${hiz_occlusion_wgsl}
${entries_wgsl}
${order_wgsl}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let a = gid.x + gid.y * u32(uParams.uRowThreads);
    if (a >= u32(uParams.uActiveCount)) {
        return;
    }
    let c = activeCopy(a);
    // the fallback's and the second phase's: copies alone
    if (c == NoIndex && (uParams.uFallback != 0 || uParams.uSelectAppend != 0)) {
        return;
    }
    if (c != NoIndex) {
        let w = activeWord(a);
        // a leaf page not written: not drawn (its copy held coarse where it lacks it in view, \`paged-cull-wants\`)
        if (valueBlockOf(w & ${PagedEntryPageMask}u) == NoIndex) {
            return;
        }
        let levels = (w >> ${PagedEntryLevelsShift}u) & 7u;
        var level = copyLevel(c, levels);
        if (levels > 0u && lacking[copyGroup(c)] != 0u) {
            level = max(level, 1u);
        }
        // the fallback's: a dropped copy group's entries a level coarser, its coarsest at that again (what of it the
        // phases drew drawn twice); none without coarse levels
        if (uParams.uFallback != 0) {
            if (levels == 0u || atomicLoad(&marks[u32(uParams.uDroppedBase) + copyGroup(c)]) == 0u) {
                return;
            }
            level = min(levels, level + 1u);
        }
        if ((w >> ${PagedEntryLevelShift}u) != level) {
            return;
        }
        let b = copyBounds(c);
        let known = uParams.uGhost == 0 && b.x >= 0.0 && b.y >= 0.0;
        if (uParams.uSelectAppend != 0 && !known) {
            return;
        }
        if (known) {
            let s = moved(vec4<f32>(copySphere(c).xyz, b.x + uParams.uMaxReach * b.y));
            if (outside(uParams.uPlane0, s.xyz, s.w) || outside(uParams.uPlane1, s.xyz, s.w) ||
                outside(uParams.uPlane2, s.xyz, s.w) || outside(uParams.uPlane3, s.xyz, s.w) ||
                outside(uParams.uPlane4, s.xyz, s.w) || outside(uParams.uPlane5, s.xyz, s.w)) {
                return;
            }
            // a copy the pyramid hides whole: the first phase's left out, marked (\`uCopyMarksBase\`); the second
            // phase's those marked the frame's pyramid shows, appended to the first's
            if (uParams.uFallback == 0 && uParams.uHiZ != 0u) {
                let mark = u32(uParams.uCopyMarksBase) + c;
                if (uParams.uSelectAppend != 0) {
                    if (atomicLoad(&marks[mark]) == 0u || occluded(s)) {
                        return;
                    }
                } else if (occluded(s)) {
                    atomicStore(&marks[mark], 1u);
                    return;
                }
            }
        }
    }
    // its bucket's count (the first pass), else its place from its bucket's cursor (the second: the counts made the
    // buckets' firsts by \`paged-cull-select-args\`)
    let bucket = 4u + bucketOf(c);
    if (uParams.uSelectPass == 0) {
        atomicAdd(&selectArgs[bucket], 1u);
        return;
    }
    let k = atomicAdd(&selectArgs[bucket], 1u);
    let i = 4u + k;
    tables[u32(uParams.uSelectedBase) + i / 4u][i % 4u] = a;
}
`;

/**
 * The leaf pages the view wants (\`paged-cull-wants\`, a thread an active entry, before the select): a leaf entry
 * (its template's with coarse levels, level 0) of a copy at level 0 whose page is not written, its record's box (in
 * the template's frame) placed by the copy, padded by the largest reach times the block's largest scale: in the frustum
 * its page wanted (\`wanted\`, a bit a page) and its copy's group lacking it (\`lacking\`: held at its first coarse
 * level, no holes); grown by \`uWantPadding\` of its radius its page wanted alone (written before it comes into
 * view). Read back, the CPU writes the wanted pages (see \`paged-cull.ts\`). Not by the pyramid: a page's sphere hid
 * little (the yeast 2x2's a quarter 539 MB either way), and a copy held coarse drew a pyramid of a sample's gaps,
 * wanting more pages (a feedback)
 */
export const paged_cull_wants_wgsl = /* wgsl */ `
const NoIndex = 0xffffffffu;
${entries_wgsl}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let a = gid.x + gid.y * u32(uParams.uRowThreads);
    if (a >= u32(uParams.uActiveCount)) {
        return;
    }
    let c = activeCopy(a);
    if (c == NoIndex) {
        return;
    }
    let w = activeWord(a);
    let levels = (w >> ${PagedEntryLevelsShift}u) & 7u;
    if (levels == 0u || (w >> ${PagedEntryLevelShift}u) != 0u || copyLevel(c, levels) != 0u) {
        return;
    }
    let page = w & ${PagedEntryPageMask}u;
    // a page written wanted too (kept: \`PagedLeafReleaseMs\`), its copy lacking none
    let written = valueBlockOf(page) != NoIndex;
    let lo = bitcast<vec4<f32>>(tables[page * 4u + 2u]).xyz;
    let hi = bitcast<vec4<f32>>(tables[page * 4u + 3u]).xyz;
    if (!all(lo <= hi)) {
        return;
    }
    let scale = copyBounds(c).y;
    let pad = uParams.uMaxReach * select(1.0, max(scale, 1.0), scale >= 0.0);
    let s = moved(vec4<f32>(copyPoint(c, (lo + hi) * 0.5), length(hi - lo) * 0.5 + pad));
    let bit = 1u << (page % 32u);
    if (inFrustum(s)) {
        atomicOr(&wanted[page / 32u], bit);
        if (!written) {
            atomicStore(&lacking[copyGroup(c)], 1u);
        }
    } else if (inFrustum(vec4<f32>(s.xyz, s.w * (1.0 + uParams.uWantPadding)))) {
        atomicOr(&wanted[page / 32u], bit);
    }
}
`;

/**
 * The groups the cull runs over (\`paged-cull-groups\`, a thread a selected entry's group, see \`paged-cull.ts\`):
 * tested whole (\`groupState\`), those not culled whole appended to the survivors (\`uSurvivorsBase\` in \`bases\`: their
 * count, then an entry's group \`k * GroupsPerPage + local\` each, the top bit where all its particles are tested)
 */
export const paged_cull_groups_wgsl = /* wgsl */ `
const NoIndex = 0xffffffffu;
const PageSize = ${PagedCullPageSize}u;
const WgKeys = ${PagedCullWorkgroupKeys}u;
const GroupsPerPage = PageSize / WgKeys;
${hiz_occlusion_wgsl}
${entries_wgsl}
${order_wgsl}
${phases_wgsl}
${group_state_wgsl}

@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let t = gid.x + gid.y * ${PagedCullMaxRow * PagedCullGroupsWorkgroup}u;
    let k = t / GroupsPerPage;
    if (k >= tables[u32(uParams.uSelectedBase)].x) {
        return;
    }
    let local = t % GroupsPerPage;
    let state = groupState(k, local);
    if (state == 1u) {
        return;
    }
    // by its entry's bucket (the select's order: a pass counting, the args their firsts, a pass placing), the counts
    // and cursors the region's first words (not \`selectArgs\`: this pass's indirect dispatch). Why: appended at one
    // cursor the survivors came in the order the workgroups ran, the near copies first no more (the isg 16x16x16's
    // first frames of a pyramid 17 M frame slots, not 9 M: the rooms grown by them)
    let base = u32(uParams.uSurvivorsBase);
    let bucket = base + 1u + bucketOf(activeCopy(selectedEntry(k)));
    if (uParams.uGroupsPass == 0) {
        atomicAdd(&survivors[bucket], 1u);
        return;
    }
    let i = atomicAdd(&survivors[bucket], 1u);
    // past the room (more entries selected than it holds, read back and grown: \`readSelected\`) none
    if (i < u32(uParams.uSurvivorsRoom)) {
        let item = (k * GroupsPerPage + local) | select(0u, 0x80000000u, state == 2u);
        atomicStore(&survivors[base + 1u + OrderBuckets + i], item);
    }
}
`;

/**
 * The selected entries' count into the tables' header and the cull's indirect dispatch (\`selectArgs\` from its second
 * word): a workgroup a group of a selected entry's page, rows of \`PagedCullMaxRow\`
 */
export const paged_cull_select_args_wgsl = /* wgsl */ `
const row = ${PagedCullMaxRow}u;
@compute @workgroup_size(WorkgroupSizeX)
fn main() {
    if (uParams.uArgsMode != 0) {
        cullArgs();
        return;
    }
    // the buckets' counts their firsts (the select's second pass places an entry at its bucket's), their sum the count;
    // a second phase's after the first's (\`uSelectAppend\`: the header's y the first's count)
    let base = u32(uParams.uSelectedBase);
    let first = select(0u, tables[base].x, uParams.uSelectAppend != 0);
    var count = first;
    for (var b = 0u; b < ${PagedCullOrderBuckets}u; b++) {
        let n = atomicLoad(&selectArgs[4u + b]);
        atomicStore(&selectArgs[4u + b], count);
        count += n;
    }
    atomicStore(&selectArgs[0], count);
    tables[base] = vec4<u32>(count, select(count, first, uParams.uSelectAppend != 0), 0u, 0u);
    // the group pass's: a thread a selected entry's group
    let threads = count * ${PagedCullPageSize / PagedCullWorkgroupKeys}u;
    let workgroups = (threads + ${PagedCullGroupsWorkgroup - 1}u) / ${PagedCullGroupsWorkgroup}u;
    atomicStore(&selectArgs[1], min(workgroups, row));
    atomicStore(&selectArgs[2], max(1u, (workgroups + row - 1u) / row));
    atomicStore(&selectArgs[3], 1u);
}

/**
 * the cull's dispatch (\`uArgsMode\` 1, between the group pass's counting and placing): the buckets' counts their firsts,
 * their sum the survivors' count (the header's z), a workgroup a survivor
 */
fn cullArgs() {
    let sb = u32(uParams.uSurvivorsBase);
    var total = 0u;
    for (var b = 0u; b < ${PagedCullOrderBuckets}u; b++) {
        let c = survivors[sb + 1u + b];
        survivors[sb + 1u + b] = total;
        total += c;
    }
    survivors[sb] = total;
    let n = min(total, u32(uParams.uSurvivorsRoom));
    let base = u32(uParams.uSelectedBase);
    let h = tables[base];
    tables[base] = vec4<u32>(h.x, h.y, n, 0u);
    atomicStore(&selectArgs[${PagedCullGroupArgs + 1}u], min(n, row));
    atomicStore(&selectArgs[${PagedCullGroupArgs + 2}u], max(1u, (n + row - 1u) / row));
    atomicStore(&selectArgs[${PagedCullGroupArgs + 3}u], 1u);
}
`;

/**
 * Per draw record (a level's, a member's: level-major, `j * segmentCount + s`), its instance count and first instance
 * from its region's count and first (`counts[s * levels + j]`, `bases`: a raster region's its end, its slots down from
 * it; the second phase's past the first's, `prior`; cut at its room, `bases` from `uRoomsBase`), the first its last
 * word; a member past the regions none. The
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
        // cut at the region's room (\`bases\` from \`uRoomsBase\`): none written past it
        let roomOf = bases[u32(uParams.uRoomsBase) + region];
        count = min(counts[region], roomOf - min(roomOf, past));
        if (raster) {
            first = bases[region] - min(past, roomOf) - count;
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
