/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import { Quat } from '@molstar/core/math/linear-algebra/3d/quat';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { ValueCell } from '@molstar/core/util/value-cell';
import {
  type ComputeContext,
  type ComputeKernel,
  KernelCode,
  type KernelRequest,
  type StorageBuffer,
  StorageBufferSpec,
} from '../compute.js';
import { UniformSpec } from '../renderable/schema.js';
import { getDispatchSize } from './prefix-sum.js';
import { createEmptyHiZ, createHiZValues, HiZSchema, type OcclusionPyramid, setHiZValues } from './cell-cull.js';
import type { GraphicsResources } from '../context.js';
import { ClipStorageVec4s } from '../clip-storage.js';
import { Clip } from '@molstar/core/util/clip';
import { createRangeAllocator } from '../range-allocator.js';
import { PagedSegmentMask } from '../culled-data.js';

// The culling of a page store's particles (`ParticlePageStore`): a thread a particle, nothing of
// the instances laid out on the CPU. A particle's pose composed on the GPU when its page is written (`compose`, at
// the next cull), at a slot of its page's range (a range a page, its count's: a pose a particle),
// not each frame; each frame a pass a phase: a particle's run and position as stored tested (8 B a particle, not its
// pose's 32: the yeast's pass 2.6-3.0 ms against 3.0-3.3), a
// visible instance's pose slot written into its region (a member and level, the raster's, the ghost's) at its cursor;
// the regions laid out on the CPU (a region's room by what culls needed, read back: `PagedRegionRoomShare`; a raster
// region's from its end), so no count and no prefix sum before. The draws read a pose by the slot (`CulledData`, `CulledPagedFlag`). The store's tables in one
// buffer (the records, the runs, the members, the active pages, the pages' pose slots: the storage buffers a stage has
// are few: the extras there too), its values in one (a block a page: its positions, then its rotations). A template's
// pages have no pose slots: its copies' visible instances' poses composed by the cull from the values.

/** a store page's particles as bits (`StorePageBits`): 8k, a tile its own pages (a 16³ grid's yeast tile 11-20k) */
export const PagedCullPageBits = 13;
/** particles a store page */
export const PagedCullPageSize = 1 << PagedCullPageBits;
/** a workgroup's invocations and its table of regions (the instances of a region a workgroup counted at once) */
export const PagedCullWorkgroupKeys = 256;
/** words of a page's record (`PageRecordWords`), and its runs' count's flag of its members a particle */
export const PagedRecordWords = 16;
export const PagedPerParticleMembers = 0x80000000;
/** a record's flags word, and its flag of a template's page (drawn by its copies, not an active page of its own) */
export const PagedFlagsWord = 15;
export const PagedTemplatePage = 1;
/** a page's groups (a cull workgroup's particles each), their boxes the compose's (`groupBoxes`) */
export const PagedCullGroupsPerPage = PagedCullPageSize / PagedCullWorkgroupKeys;
/** vec4s of a page's record in the tables */
const RecordVec4s = PagedRecordWords / 4;
/**
 * u32s of a composed pose (a slot's): `t` and the scale (f32), the quaternion (two snorm16x2), its store index (its
 * page's first and its index), its member (none: no instance)
 */
export const PagedPoseWords = 8;
/** a page without pose slots (the device's buffer limit reached) */
const NoSlots = 0xffffffff;
/** an active entry's copy word of a page drawn as it is (not a template's) */
export const NoCopy = 0xffffffff;
/**
 * an active entry's page word: its page in the low bits, a template's page's coarse level (0 the template, from 1 its
 * levels' samples: `ParticleNesting.getTemplateLevels`) above `PagedEntryLevelShift`, its template's coarse levels'
 * count above `PagedEntryLevelsShift` (a copy drawn at one of them by its distance, `PagedCullView.coarse`)
 */
export const PagedEntryPageMask = 0x03ffffff;
export const PagedEntryLevelsShift = 26;
export const PagedEntryLevelShift = 29;
/** a template's coarse levels at most (three bits) */
export const PagedMaxCoarseLevels = 7;
/**
 * vec4s of a copy in the tables: `t` and its scale, its quaternion, its center (where it places its template's) and
 * its template's spacing, its template's radius about that center
 */
const CopyVec4s = 4;
/**
 * the frame's slots of the copies' visible instances at most (32 B each: 1 GB): those past it not drawn (far copies
 * by a template's coarse level not yet)
 */
export const PagedFrameMaxSlots = 1 << 25;
/**
 * the frame's slots at first (64 MB), grown to what a frame's cull wrote past them (read back): the visible copies'
 * instances, not all of them (the isg cell on a 4x4x4 grid whole, 29.4 M instances: 0.56 M a still frame's, at most
 * 1.2 M rotating, the pyramid's occlusion hiding the rest)
 */
export const PagedFrameInitialSlots = 1 << 21;
/**
 * the share of the frame's slots kept for the fallback (a room below its copies' instances, a template with coarse
 * levels): a copy group whose instances found no slot in the cull's phases (`dropped`) drawn once more at a level
 * coarser into them, after the last. Why: a zoom step's first cull, of the last view's pyramid, wanted more slots than
 * the room's cap (the isg cell's 16x16x16 grid 87.6 M of 33.5 M); the near copies first (`PagedCullOrderBuckets`) its
 * first frame lacked about 220 px along the cells' seams, none with this share. Without that order the dropped copies
 * were any (a share of a third their need, at its edge)
 */
export const PagedFallbackShare = 1 / 8;
/**
 * a region's room at first (a member's at a level, its raster's): its member's instances times this share, at least
 * `PagedRegionMinRoom` of them (all where fewer), grown to what a cull wrote past it (read back). Why: a region a member's
 * every instance at every level held the instances a level a view never draws (the isg cell on a 4x4x4 grid: 4 B times
 * 29.4 M times four levels, 449 MB)
 */
export const PagedRegionRoomShare = 1 / 16;
/**
 * a region's room at first at least this (all where fewer): a floor by the regions, not the view. Why 512, not 4,096:
 * a list of thousands of members (the yeast's 5,639, four levels, the raster's regions too) held 92 M slots whatever a
 * view drew (its indices 609 MiB; the sphere raster's chunks sized by its rooms, 1.6 GB), their making a second of a
 * first near frame; its frames the same with 512 (the grown rooms what views need)
 */
export const PagedRegionMinRoom = 512;
/**
 * the regions' rooms at first at most this many slots in all (64 MB): the shares scaled down where more (each at least
 * `PagedRegionMinRoom`), grown by what culls need. Why: a share of every instance outgrew the device's buffer limit
 * with many copies (the isg cell on a 16x16x16 grid, 1.9 G instances: a sixteenth at nine regions a member past 2^29
 * slots, the regions past the limit not written, nothing drawn)
 */
export const PagedRegionInitialSlots = 1 << 24;
/**
 * the selected entries the occlusion bits hold at first (1 KB an entry: 4 MB), grown to half again a cull's count past
 * them (read back); an entry past them not tested against the pyramid (drawn in the first phase: culled less, none
 * hidden). Why: bits for every active entry (every copy's pages at every level) were most of the cull's memory with
 * many copies, while the cull runs over the selected alone (the isg cell on a 16x16x16 grid: 285 MB)
 */
export const PagedOcclusionInitialEntries = 4096;
/**
 * a template's leaf pages (a template with coarse levels, `level` 0) written where the cull wants them
 * (`paged-cull-wants`: a page of a copy at level 0 in the frustum, read back); kept `PagedLeafReleaseMs` after the last
 * cull that wanted them. Why: written whole at the first sync, the yeast cell's 31.5 M particles' positions and
 * rotations (576 MB, 0.8 s of its load) while its whole view draws its coarse levels alone
 */
/** a leaf page let go this long after the last cull that wanted it (its value block) */
export const PagedLeafReleaseMs = 10_000;
/** the cull's wants read back at most this often where none are pending (`readWants`) */
export const PagedWantsReadMs = 250;
/**
 * a leaf page wanted where its sphere, grown by this share of its radius, is in view (written before it is), its copy
 * held coarse where it is in view itself
 */
export const PagedLeafPagePadding = 0.5;
/**
 * a wanted leaf block's pages written at most this many bytes a cull (their positions and rotations), the view culled
 * again till all are, its copies drawn at their first coarse level meanwhile (`copyLeastLevel`); resident once all are.
 * Why: written in one cull, the yeast cell's leaves (576 MB) a frame of about a second as a zoom came near
 */
export const PagedLeafUploadBytes = 32 * 1024 * 1024;
/**
 * the value blocks' bytes at most (every page's written: leaf and not), within the device's storage buffer limit: past
 * it the leaf pages no read of the wants wanted last let go for room (the least recently wanted first), a wanted page
 * without room left unwritten (its copies drawn coarse). Why: the values' buffer one binding, the device's limit its
 * size (grown past it, the writes out of bounds)
 */
export const PagedValueBytes = Number.POSITIVE_INFINITY;
/** the transient transform of an object at rest */
const IdentityTransform = Mat4.identity();

/** a workgroup a group of a page's particles (the cull's): its centers' box reduced (`groupBoxes`) */
export const PagedCullComposeKernel = KernelCode('paged-cull-compose', [PagedCullWorkgroupKeys]);
export const PagedCullKernel = KernelCode('paged-cull', [PagedCullWorkgroupKeys]);
export const PagedCullArgsKernel = KernelCode('paged-cull-args', [64]);
/**
 * the entries the cull runs over (`paged-cull-select`, a thread an active entry): a page drawn as it is, a template's
 * page at its copy's level (`PagedCullView.coarse`) where the copy's sphere is in the frustum, compacted into the
 * tables' selected section (`uSelectedBase`), counted in `selectArgs`; then (`paged-cull-select-args`, one thread) the
 * count into the section's header and the cull's indirect dispatch (`selectArgs` from its second word). Why: the cull
 * a workgroup per group of every copy's every level's page, those of another level or out of view leaving at once,
 * was most of a frame with many copies (the isg cell on a 16x16x16 grid: 66 of 69 ms rotating)
 */
export const PagedCullSelectKernel = KernelCode('paged-cull-select', [64]);
/**
 * the selected entries' order (`paged-cull-select`, a counting sort: a pass counting per bucket, one thread their
 * firsts, a pass placing): a copy's by its center's depth over the copies' range in the view (`copyDepthSphere`),
 * this many buckets, near first. Why: a cull's workgroups take the frame's slots about in their order, a room
 * overflowed then drops the far copies (most hidden by the near), not whichever came late (a zoom step's first cull of
 * the isg cell's 16x16x16 grid: 87.6 M of 33.5 M, a third of the copies in view drawn, holes)
 */
export const PagedCullOrderBuckets = 256;
/**
 * the groups the cull runs over (`paged-cull-groups`, a thread a selected entry's group: tested whole, those not culled
 * whole appended to the survivors, the cull a workgroup a survivor), this many threads a workgroup. Why: the cull a
 * workgroup a group of every selected page, those culled whole leaving at once, still launched them and their
 * barriers (the isg 16x16x16 near: 1.05 M a phase, 630 k of them culled whole)
 */
export const PagedCullGroupsKernel = KernelCode('paged-cull-groups', [64]);
export const PagedCullGroupsWorkgroup = 64;
/** `selectArgs`' words of the cull's indirect dispatch less one (after the order buckets') */
export const PagedCullGroupArgs = 4 + PagedCullOrderBuckets;
/** the leaf pages the view wants (`paged-cull-wants`, a thread an active entry, before the select) */
export const PagedCullWantsKernel = KernelCode('paged-cull-wants', [64]);
export const PagedCullSelectArgsKernel = KernelCode('paged-cull-select-args', [1]);
/** the cull's rows of workgroups (its dispatch's x where in 2d, `paged-cull-select-args`) */
export const PagedCullMaxRow = 65535;

const PagedCullUniforms = {
  /**
   * the object's transient transform (`uObjectTransform`, a drag's preview) and its largest scale: the spheres tested
   * moved by it, as the draws move the object
   */
  uObjectTransform: UniformSpec('m4'),
  uObjectScale: UniformSpec('f'),
  uPlane0: UniformSpec('v4'),
  uPlane1: UniformSpec('v4'),
  uPlane2: UniformSpec('v4'),
  uPlane3: UniformSpec('v4'),
  uPlane4: UniformSpec('v4'),
  uPlane5: UniformSpec('v4'),
  uCameraPlane: UniformSpec('v4'),
  /** per level (column) its distance range in x and y */
  uLevels: UniformSpec('m4'),
  /** 0 without levels (all at level 0) */
  uLevelCount: UniformSpec('i'),
  /**
   * the active entries: the dispatch's pages, a page's particles a thread each; an entry its page and its copy (two
   * words: `NoCopy` a page drawn as it is, else a template's page placed by its copy, `uCopiesBase`)
   */
  uActiveCount: UniformSpec('i'),
  /** threads a row of workgroups (the dispatch in 2d past the per-dimension limit) */
  uRowThreads: UniformSpec('i'),
  /** the tables' sections, in vec4s: the runs, the members (center, reach), the active pages and the pages' first pose slots (four a vec4) */
  uRunsBase: UniformSpec('i'),
  uMembersBase: UniformSpec('i'),
  uActiveBase: UniformSpec('i'),
  /** a page's value block in `positions` (a word a page, `NoSlots` none: not written), the tables' section */
  uValueBlocksBase: UniformSpec('i'),
  uPoseFirstsBase: UniformSpec('i'),
  /** the extras' section (two words a particle: member, scale), in vec4s */
  uExtrasBase: UniformSpec('i'),
  /**
   * the selected entries (`PagedCullSelectKernel`): a header vec4 (their count in x), then an active entry's index a
   * word, in vec4s; the cull's workgroups by them (a group of a selected entry's page each)
   */
  uSelectedBase: UniformSpec('i'),
  /**
   * the copies (four vec4s a copy: `t` and its scale, its quaternion, its center and its template's spacing, its
   * template's radius about its center and its block's largest scale (below 0: none known); f32 bits), in vec4s
   */
  uCopiesBase: UniformSpec('i'),
  /**
   * a copy's coarse level (`PagedCullView.coarse`): a template's spacing at its center's distance this many pixels a
   * unit (at a distance of 1; orthographic at any, `uCoarseOrtho`), a level coarser each halving below
   * `uCoarsePixels` (0: none)
   */
  uCoarseScale: UniformSpec('f'),
  uCoarseOrtho: UniformSpec('i'),
  uCoarsePixels: UniformSpec('f'),
  /**
   * the frame's slots of a copy's visible instances (their poses composed by the cull: the copy's and its template
   * particle's), from slot `uFrameBase` (the persistent slots' room), `uFrameRoom` of them (past it not drawn); their
   * cursor `counts[uFrameCursor]` (after the regions')
   */
  uFrameBase: UniformSpec('i'),
  uFrameRoom: UniformSpec('i'),
  uFrameCursor: UniformSpec('i'),
  /**
   * the fallback's pass (1, `PagedFallbackShare`): the dropped copy groups' entries (the select's), its frame slots
   * `uFallbackRoom` after `uFrameRoom`, their cursor `counts[uFrameCursor + 1]`; its regions' slots after the last
   * phase's (its counts not cleared), not by the ghost, nothing marked
   */
  uFallback: UniformSpec('i'),
  uFallbackRoom: UniformSpec('i'),
  /** a word a copy group in `occludedBits` from this word: its instances past the frame's room in a phase (`dropped`) */
  uDroppedBase: UniformSpec('i'),
  /** positions: 0 a grid's words (two a particle), 1 float32 (three) */
  uPositions: UniformSpec('i'),
  /** rotations: 0 none, 1 float32 (four words a particle), 2 snorm16 (two); after a page's positions (`positions`) */
  uRotations: UniformSpec('i'),
  /** 0 the first phase (the last frame's pyramid), 1 the second (the marked alone, the frame's pyramid so far) */
  uPhase: UniformSpec('i'),
  /** the clip objects (`clip_storage_wgsl`, in the tables from vec4 `uClipBase`): spheres they clip whole culled */
  uClipObjectCount: UniformSpec('i'),
  uClipBase: UniformSpec('i'),
  /** instance clipping: an instance's center tested too */
  uClipInstance: UniformSpec('i'),
  /** the ghost of what the clip objects remove (`Renderer.renderClipGhost`): a region a member from `uGhostBase` */
  uGhost: UniformSpec('i'),
  uGhostBase: UniformSpec('i'),
  /** the members of the regions */
  uMemberCount: UniformSpec('i'),
  /** the poses' first vec4 (after their `CulledData`'s header and segments, `setPosesFront`) */
  uPoseBase: UniformSpec('i'),
  /** the pose slots (the device's buffer limit): those past it none, not drawn */
  uPoseRoom: UniformSpec('i'),
  /** the indices' slots (the device's buffer limit): those past it not written, not drawn */
  uIndexRoom: UniformSpec('i'),
  /**
   * the sphere raster's split (`MergedRasterParams`, a radius of 0 none): an instance it draws into its level's raster
   * region (from `uRasterBase`, level-major; none where it is `uGhostBase`) by its member's raster row (the tables from
   * vec4 `uRasterRowsBase`)
   */
  uRasterRadius: UniformSpec('f'),
  uRasterScale: UniformSpec('f'),
  uRasterNear: UniformSpec('f'),
  /** 1 orthographic: a sphere's size on screen the same at every distance (`MergedRasterParams.ortho`) */
  uRasterOrtho: UniformSpec('f'),
  uRasterBase: UniformSpec('i'),
  uRasterRowsBase: UniformSpec('i'),
  /**
   * the largest member's reach: a group's sphere padded by it times the group's largest scale (`groupBoxes`' w, the
   * compose's), a copy's by it times its block's (the copies' fourth vec4's y), what a group or a copy tested whole
   * takes (none hidden whose particles show). Not times the largest scale of all: a template's coarse level's (its
   * particles drawn up to 63 times larger) padded every group and copy by it, none culled whole
   */
  uMaxReach: UniformSpec('f'),
  /** the groups' bits in `occludedBits` from this word: a group occluded whole in the first phase */
  uPageBitsBase: UniformSpec('i'),
  /** the selected entries with occlusion bits (`PagedOcclusionInitialEntries`): those past it not tested by the pyramid */
  uBitsRoom: UniformSpec('i'),
  /** the regions (`regions`: their firsts, then the first phase's counts) */
  uRegionCount: UniformSpec('i'),
  /**
   * the survivors in `regions` (`bases`) from this word (`paged-cull-groups`: their count, the order buckets' cursors,
   * then a word each), at most `uSurvivorsRoom` of them
   */
  uSurvivorsBase: UniformSpec('i'),
  uSurvivorsRoom: UniformSpec('i'),
  /** the group pass's: 0 counting by bucket, 1 placing */
  uGroupsPass: UniformSpec('i'),
  /** a copy's bucket: its center's depth less `uOrderNear`, times `uOrderScale` (the select's and the group pass's) */
  uOrderNear: UniformSpec('f'),
  uOrderScale: UniformSpec('f'),
  ...HiZSchema,
};

/** The poses of the pages in `composePages` (a page and its first pose slot each), from their values as stored */
export const PagedCullComposeSchema = {
  uRowThreads: UniformSpec('i'),
  uComposeCount: UniformSpec('i'),
  uRunsBase: UniformSpec('i'),
  uMembersBase: UniformSpec('i'),
  uExtrasBase: UniformSpec('i'),
  uPoseBase: UniformSpec('i'),
  uPoseRoom: UniformSpec('i'),
  uValueBlocksBase: UniformSpec('i'),
  /** positions: 0 a grid's words (two a particle), 1 float32 (three) */
  uPositions: UniformSpec('i'),
  /** rotations: 0 none, 1 float32 (four words a particle), 2 snorm16 (two) */
  uRotations: UniformSpec('i'),
  tables: StorageBufferSpec('vec4<u32>', 'read'),
  /** a page's values as stored: its positions, then its rotations (`uRotations`) */
  positions: StorageBufferSpec('u32', 'read'),
  composePages: StorageBufferSpec('u32', 'read'),
  poses: StorageBufferSpec('vec4<u32>', 'read-write'),
  /** per group (`PagedCullGroupsPerPage` a page) its instances' centers' box, two vec4s (min, max) */
  groupBoxes: StorageBufferSpec('vec4<f32>', 'read-write'),
};

export const PagedCullSchema = {
  ...PagedCullUniforms,
  tables: StorageBufferSpec('vec4<u32>', 'read'),
  positions: StorageBufferSpec('u32', 'read'),
  /** the regions' cursors (the phase's) */
  counts: StorageBufferSpec('atomic<u32>', 'read-write'),
  /**
   * the regions' firsts (a raster region's its end, laid out on the CPU), then the first phase's counts (copied there
   * for the second, its slots past them): one binding, the stage's eight
   */
  regions: StorageBufferSpec('u32', 'read'),
  /** the groups' boxes (the compose's) */
  groupBoxes: StorageBufferSpec('vec4<f32>', 'read'),
  /**
   * a bit an entry's particle: occluded in the first phase (`uPhase` 0), what the second tests; a bit its group's; a
   * word a copy group dropped (`uDroppedBase`)
   */
  occludedBits: StorageBufferSpec('atomic<u32>', 'read-write'),
  /** the visible instances' pose slots, by region */
  indices: StorageBufferSpec('u32', 'read-write'),
  /** the poses: a copy's visible instances' composed into the frame's slots (the persistent read) */
  poses: StorageBufferSpec('vec4<u32>', 'read-write'),
};

/**
 * The draw records of a page store's regions (`paged-cull-args`): a record a level and segment (level-major), its
 * instance count and first instance from its region's count and first
 */
export const PagedCullArgsSchema = {
  uRecordCount: UniformSpec('i'),
  uRecordWords: UniformSpec('i'),
  uSegmentCount: UniformSpec('i'),
  /** the regions a member (at least one) */
  uLevels: UniformSpec('i'),
  /** the members the regions are of (the segments past them none); the records past the levels' the ghost's */
  uMemberCount: UniformSpec('i'),
  /** the second phase's: its regions' firsts past the first's instances (`prior`) */
  uSecond: UniformSpec('i'),
  /** the indices' slots: a region's count cut where its slots would pass it (not written there) */
  uIndexRoom: UniformSpec('i'),
  /** the raster's blocks of records after the levels' (0 or the levels), the ghost's after them */
  uRasterLevels: UniformSpec('i'),
  /** the regions' rooms in `bases` (after their firsts and the first phase's counts: twice the regions) */
  uRoomsBase: UniformSpec('i'),
  counts: StorageBufferSpec('u32', 'read'),
  bases: StorageBufferSpec('u32', 'read'),
  prior: StorageBufferSpec('u32', 'read'),
  args: StorageBufferSpec('u32', 'read-write'),
};

/** The selected entries of the cull (`PagedCullSelectKernel`) */
export const PagedCullSelectSchema = {
  uRowThreads: UniformSpec('i'),
  uActiveCount: UniformSpec('i'),
  uActiveBase: UniformSpec('i'),
  uCopiesBase: UniformSpec('i'),
  uSelectedBase: UniformSpec('i'),
  uCoarseScale: UniformSpec('f'),
  uCoarseOrtho: UniformSpec('i'),
  uCoarsePixels: UniformSpec('f'),
  uCameraPlane: UniformSpec('v4'),
  uObjectTransform: UniformSpec('m4'),
  uObjectScale: UniformSpec('f'),
  uPlane0: UniformSpec('v4'),
  uPlane1: UniformSpec('v4'),
  uPlane2: UniformSpec('v4'),
  uPlane3: UniformSpec('v4'),
  uPlane4: UniformSpec('v4'),
  uPlane5: UniformSpec('v4'),
  uMaxReach: UniformSpec('f'),
  uGhost: UniformSpec('i'),
  uValueBlocksBase: UniformSpec('i'),
  tables: StorageBufferSpec('vec4<u32>', 'read-write'),
  selectArgs: StorageBufferSpec('atomic<u32>', 'read-write'),
  /** per copy group whether it lacks a leaf page in view (`paged-cull-wants`): held at its first coarse level */
  lacking: StorageBufferSpec('u32', 'read'),
  /** the fallback's (1): a dropped copy group's entries alone, a level coarser (`dropped`, from `uDroppedBase`) */
  uFallback: UniformSpec('i'),
  /**
   * the second phase's (1): the copies the first left out (occluded whole by the last frame's pyramid, marked from
   * `uCopyMarksBase`) that the frame's shows, appended to the first's entries (the header's y their first)
   */
  uSelectAppend: UniformSpec('i'),
  uCopyMarksBase: UniformSpec('i'),
  /** 0 the pass counting the buckets' entries, 1 the one placing them (`PagedCullOrderBuckets`) */
  uSelectPass: UniformSpec('i'),
  /** a copy's bucket: its center's depth less `uOrderNear`, times `uOrderScale` */
  uOrderNear: UniformSpec('f'),
  uOrderScale: UniformSpec('f'),
  uDroppedBase: UniformSpec('i'),
  /**
   * the cull's `occludedBits`: a word a copy group dropped from `uDroppedBase`, a word a copy occluded whole by the
   * last frame's pyramid from `uCopyMarksBase`
   */
  marks: StorageBufferSpec('atomic<u32>', 'read-write'),
  /** a copy tested whole against the pyramid (`uSelectAppend`) */
  ...HiZSchema,
};

/** The leaf pages a view wants (`PagedCullWantsKernel`) */
export const PagedCullWantsSchema = {
  uRowThreads: UniformSpec('i'),
  uActiveCount: UniformSpec('i'),
  uActiveBase: UniformSpec('i'),
  uCopiesBase: UniformSpec('i'),
  uSelectedBase: UniformSpec('i'),
  uValueBlocksBase: UniformSpec('i'),
  uCoarseScale: UniformSpec('f'),
  uCoarseOrtho: UniformSpec('i'),
  uCoarsePixels: UniformSpec('f'),
  uCameraPlane: UniformSpec('v4'),
  uObjectTransform: UniformSpec('m4'),
  uObjectScale: UniformSpec('f'),
  uPlane0: UniformSpec('v4'),
  uPlane1: UniformSpec('v4'),
  uPlane2: UniformSpec('v4'),
  uPlane3: UniformSpec('v4'),
  uPlane4: UniformSpec('v4'),
  uPlane5: UniformSpec('v4'),
  uMaxReach: UniformSpec('f'),
  /** a page wanted (not lacking) where grown by this share of its radius it is in view (`PagedLeafPagePadding`) */
  uWantPadding: UniformSpec('f'),
  tables: StorageBufferSpec('vec4<u32>', 'read'),
  /** a bit a store page: wanted */
  wanted: StorageBufferSpec('atomic<u32>', 'read-write'),
  /** per copy group: lacking a leaf page in view */
  lacking: StorageBufferSpec('atomic<u32>', 'read-write'),
};

/** The groups the cull runs over (`PagedCullGroupsKernel`) */
export const PagedCullGroupsSchema = {
  ...PagedCullUniforms,
  tables: StorageBufferSpec('vec4<u32>', 'read'),
  groupBoxes: StorageBufferSpec('vec4<f32>', 'read'),
  occludedBits: StorageBufferSpec('atomic<u32>', 'read-write'),
  /** `bases`: the survivors from `uSurvivorsBase` (their count, then a word each) */
  survivors: StorageBufferSpec('atomic<u32>', 'read-write'),
};

/** The selected entries' count and the group pass's indirect dispatch, or (`uArgsMode` 1) the cull's */
export const PagedCullSelectArgsSchema = {
  uSelectedBase: UniformSpec('i'),
  uSelectAppend: UniformSpec('i'),
  uArgsMode: UniformSpec('i'),
  /** the cull's (`uArgsMode` 1): the survivors' buckets' counts their firsts, their count, at most `uSurvivorsRoom` */
  uSurvivorsBase: UniformSpec('i'),
  uSurvivorsRoom: UniformSpec('i'),
  survivors: StorageBufferSpec('u32', 'read-write'),
  tables: StorageBufferSpec('vec4<u32>', 'read-write'),
  selectArgs: StorageBufferSpec('atomic<u32>', 'read-write'),
};

export const PagedCullKernels: readonly KernelRequest[] = [
  { code: PagedCullComposeKernel, schema: PagedCullComposeSchema },
  { code: PagedCullSelectKernel, schema: PagedCullSelectSchema },
  { code: PagedCullSelectArgsKernel, schema: PagedCullSelectArgsSchema },
  { code: PagedCullGroupsKernel, schema: PagedCullGroupsSchema },
  { code: PagedCullWantsKernel, schema: PagedCullWantsSchema },
  { code: PagedCullKernel, schema: PagedCullSchema },
  { code: PagedCullArgsKernel, schema: PagedCullArgsSchema },
];

/** What the culler reads of a page store (`ParticlePageStore` is one) */
export interface PagedCullSource {
  readonly format:
    | { readonly positions: 'grid' | 'float32'; readonly rotations: 'none' | 'float32' | 'snorm16' }
    | undefined;
  /** the page records, `PageRecordWords` words a page */
  readonly records: Uint32Array;
  readonly pageCapacity: number;
  /** the runs, four words a run */
  readonly runs: Uint32Array;
  readonly runCapacity: number;
  /** the extras, two words a particle */
  readonly extras: Uint32Array;
  readonly extraCapacity: number;
  /**
   * a page's values not read where it has none and is a deferred list's (a pack opened by range): its read begun
   * (once while in flight), undefined where it has its values or none to read
   */
  readonly requestValues?: (page: number) => Promise<void> | undefined;
  /** a page's values given up once written (a deferred list's: its copy on the GPU, read again when wanted) */
  readonly releaseValues?: (page: number) => void;
  /** a page's values as stored (none: free, or a deferred list's not read: `requestValues`) */
  readonly values: readonly (
    | { readonly positions: ArrayBufferView; readonly rotations?: ArrayBufferView }
    | undefined
  )[];
  readonly pageCount: number;
  /** the members: their centers (three a member) and their count */
  readonly pool: { readonly centers: Float64Array; readonly count: number; readonly counts: Int32Array };
  /**
   * the templates (a nested list's): their pages (`PagedTemplatePage`, not active of their own) and their copies'
   * poses (eight floats a copy: `t` and its scale, the quaternion), each page drawn once a copy; a template's coarse
   * level a template of its own (`level` from 1, the copies its template's), drawn by a copy at that level alone, of
   * `levels` levels (the template's and each of its levels'), the template's `center` (in its frame) and `spacing`
   * what picks a copy's level
   */
  readonly templates?: readonly {
    readonly pages: ArrayLike<number>;
    readonly copies: {
      readonly poses: Float32Array;
      readonly level?: number;
      readonly levels?: number;
      readonly center?: ArrayLike<number>;
      readonly spacing?: number;
      /** its template's particles' centers' distance from `center` at most (none: a copy not culled whole) */
      readonly radius?: number;
    };
  }[];
}

/** What changed in the source since the last sync (`PageStoreChanges`) */
export interface PagedCullChanges {
  readonly written: readonly number[];
  readonly released: readonly number[];
  readonly runs: readonly number[];
  readonly extras: readonly number[];
  /** the templates' copies set or let go */
  readonly copies?: readonly number[];
  /** the members whose instances or reach changed (none written, released nor grown with them) */
  readonly members?: ReadonlySet<number>;
  readonly grown: boolean;
}

/** The camera's view a cull is of, and the levels' distance ranges */
export interface PagedCullView {
  readonly frustum: Frustum3D;
  readonly cameraPlane: Plane3D;
  /** per level (column) its distance range in x and y */
  readonly levels: Mat4;
  /** 0 without levels */
  readonly levelCount: number;
  /** the ghost of what the clip objects remove too, into a region a member after the levels' and the raster's */
  readonly ghost?: boolean;
  /** the raster's regions (the levels' again, level-major after theirs): its split's instances into them (`raster`) */
  readonly rasterRegions?: boolean;
  /** the sphere raster's split (`MergedRasterParams`), with `rasterRegions` */
  readonly raster?: { readonly radius: number; readonly scale: number; readonly near: number; readonly ortho: boolean };
  /** the object's transient transform (`uObjectTransform`), the identity where none */
  readonly objectTransform?: Mat4;
  /**
   * a pyramid to come (Hi-Z on, none built yet: an app's first frames): the cull's frame slots overflowing not grown
   * for (nothing occluded: what it needs not what the view will), the view culled again once one is (`frameVersion`)
   */
  readonly pyramidToCome?: boolean;
  /**
   * a second phase to follow a cull with a pyramid (`cullSecond`, Hi-Z's two phases on): none, the first's slots the
   * frame's (its cursors read back, its fallback), its occluded left out
   */
  readonly secondPhase?: boolean;
  /**
   * a copy drawn at its template's coarse level (`PagedEntryLevelShift`) where its template's spacing at its center's
   * distance is below `pixels` (a level coarser each halving: its sample's spacing twice its level's before): `scale`
   * pixels a unit at a distance of 1 (orthographic at any); none: every copy the template whole
   */
  readonly coarse?: { readonly scale: number; readonly ortho: boolean; readonly pixels: number };
}

/** The clip objects a cull is of (`packStorageClipObjects`'s, `ClipStorageVec4s` a clip object) */
export interface PagedCullClip {
  readonly data: Float32Array;
  readonly count: number;
  /** instance clipping: an instance's center tested too */
  readonly instance: boolean;
}

/** What a paged cull's read-backs last found (`readCursors`, `readSelected`), counted since it was made */
export interface PagedCullReadBack {
  /** the cursors' reads, those that overflowed the frame's slots or a region's */
  reads: number;
  overflows: number;
  /** the frame's instances the last read needed, the phases' room then, the frame's room, whether at its cap */
  frameNeed: number;
  phasesRoom: number;
  frameRoom: number;
  frameCapped: boolean;
  /** the regions the last read found past their rooms, the most a region needed past its room */
  regionsOver: number;
  regionExcess: number;
  /** the selected entries last read, the occlusion bits' and the survivors' rooms then */
  selected: number;
  bitsRoom: number;
  survivorsEntries: number;
}

export interface PagedCuller {
  /**
   * The source's changes uploaded: its pages' records, runs and extras; its members' centers and `reaches` (a member's
   * geometry's reach from its center, a particle's sphere its scale times it); the pages written given pose slots and
   * composed at the next cull
   */
  sync: (source: PagedCullSource, changes: PagedCullChanges, reaches: ArrayLike<number>) => void;
  /**
   * The pages pending composed, then the visible instances of `view` written into their regions (`indices`, their
   * cursors `counts`): in compute passes of their own (ended)
   */
  cull: (view: PagedCullView, pyramid?: OcclusionPyramid) => void;
  /** the clip objects of the next culls (in the tables: kept through their layouts) */
  setClip: (clip: PagedCullClip) => void;
  /**
   * per member its raster row (a vec4 a member: its largest sphere at each level's scale, `rasterized`), kept through
   * the tables' layouts
   */
  setRasterRows: (rows: Float32Array) => void;
  /**
   * the vec4s before the poses in `poses` (their `CulledData`'s header and segments, its slots' first their count):
   * another count of them moves the poses (a copy)
   */
  setPosesFront: (front: Uint32Array) => void;
  /**
   * The second phase of the last `cull` (it had a pyramid): the instances the first occluded tested against `pyramid`
   * (of the frame's depth so far), the visible written into the same regions past the first's (`secondCounts` their
   * cursors from 0, `counts` the first's); false where nothing ran
   */
  cullSecond: (pyramid: OcclusionPyramid) => boolean;
  /**
   * The second phase of the last `cull` not to run (no pyramid of the frame, nothing drawn before it): the first's
   * cursors read back (its rooms grown, the view culled again where they overflowed); its fallback not (its records
   * dispatched before)
   */
  endFirst: () => void;
  /** the last select's entries in their order (`PagedCullOrderBuckets`): a page and its copy (`NoCopy` none) each */
  readSelectedEntries: () => Promise<{ page: number; copy: number }[]>;
  readonly counts: StorageBuffer;
  readonly secondCounts: StorageBuffer;
  /**
   * the regions' firsts in `indices` (a raster region's its end), the first phase's counts after them (copied there by
   * `cullSecond`), the regions' rooms after those (`regionRooms`)
   */
  readonly bases: StorageBuffer;
  /** the visible instances' pose slots by region (the regions' rooms one after another) */
  readonly indices: StorageBuffer;
  /** the slots of `indices` written (within the device's buffer limit), as of the last cull */
  readonly indexRoom: number;
  /** the composed poses, `PagedPoseWords` u32s a slot, after the front (`setPosesFront`) */
  readonly poses: StorageBuffer;
  /** the store's instances, as of the last sync */
  readonly instanceTotal: number;
  /**
   * per region its room (its slots in `indices`, `PagedRegionRoomShare`), as of the last cull; a region's instances
   * past it not written (its count still taken), the room grown and the view culled again (`frameVersion`)
   */
  readonly regionRooms: Uint32Array;
  /** the raster regions' rooms summed, as of the last cull */
  readonly rasterRoom: number;
  /** regions: members times levels (at least one), the raster's, and a member the ghost's after, as of the last cull */
  readonly regionCount: number;
  /** the members of the regions, as of the last cull */
  readonly memberCount: number;
  /** whether the last `cull` culled (a store with pages and members) */
  readonly culled: boolean;
  /**
   * advanced where a cull's instances overflowed the frame's slots or a region's (read back): its view to be culled
   * again (the room grown, or a pyramid of the frame drawn)
   */
  readonly frameVersion: number;
  /** what its read-backs last found (`PagedCullReadBack`): for diagnostics */
  readonly readBack: PagedCullReadBack;
  /** the pages of the templates' leaf blocks not written (`paged-cull-wants`), as of the last sync or cull */
  readonly deferredPages: number;
  /** the value blocks in use: the pages written (their positions and rotations in the values' buffer) */
  readonly valueBlocks: number;
  /**
   * the largest scale of the runs and extras uploaded so far (only grown): an instance's largest sphere at most its
   * member's raster row times it (a bound for a prediction on the CPU)
   */
  readonly maxScale: number;
  destroy: () => void;
}

export interface PagedCullerOptions {
  /** the frame's slots at first (`PagedFrameInitialSlots` by default) */
  readonly frameInitialSlots?: number;
  /** a region's room at first, its member's instances times this (`PagedRegionRoomShare` by default; 1: all of them) */
  readonly regionRoomShare?: number;
  /** the regions' rooms at first at most this many slots in all (`PagedRegionInitialSlots` by default) */
  readonly regionInitialSlots?: number;
  /** the selected entries the occlusion bits hold at first (`PagedOcclusionInitialEntries` by default) */
  readonly occlusionInitialEntries?: number;
  /** the selected entries whose groups the survivors hold at first (`PagedOcclusionInitialEntries` by default) */
  readonly survivorsInitialEntries?: number;
  /** a wanted leaf block's bytes written a cull at most (`PagedLeafUploadBytes` by default) */
  readonly leafUploadBytes?: number;
  /** the wants read back at most this often where none are pending (`PagedWantsReadMs` by default) */
  readonly wantsReadMs?: number;
  /** the value blocks' bytes at most (`PagedValueBytes` by default; the device's storage buffer limit at most) */
  readonly valueBytes?: number;
}

export function createPagedCuller(
  compute: ComputeContext,
  resources: GraphicsResources,
  options: PagedCullerOptions = {},
): PagedCuller {
  // a storage buffer's size a multiple of 4: of 16 (a vec4's) where that is within the device's binding limit, else of
  // 4 (the limit 2^31 - 4: rounded up to 16, 2^31, refused); a growth's 1.5 times capped at the limit by its caller
  const createBuffer = (bytes: number, label: string) => {
    const aligned = Math.ceil(bytes / 16) * 16;
    const size = aligned <= compute.maxStorageBufferSize ? aligned : Math.ceil(bytes / 4) * 4;
    return compute.createBuffer(Math.max(16, size), label);
  };
  const cell = (label: string) => ValueCell.create(createBuffer(16, label));
  const tables = cell('paged cull tables');
  /** a page's values as stored, one block a page (`valueWords`): its positions, then its rotations */
  const positions = cell('paged cull values');
  const composePages = cell('paged cull compose pages');
  const counts = cell('paged cull counts');
  const prior = cell('paged cull second counts');
  const bases = cell('paged cull bases');
  const poses = cell('paged cull poses');
  const groupBoxes = cell('paged cull group boxes');
  const indices = cell('paged cull indices');
  const occludedBits = cell('paged cull occluded');
  // the phases' cursors: the first's and the second's, swapped into `counts` and `prior` for the second
  let firstCounts = counts.ref.value;
  let secondCounts = prior.ref.value;
  let indexRoom = 0;
  const emptyHiZ = createEmptyHiZ(resources);

  /** the pyramid's values, the cull's and the select's (a copy tested whole) */
  const hiZValues = createHiZValues(emptyHiZ);
  const uniforms = {
    uObjectTransform: ValueCell.create(Mat4.identity()),
    uObjectScale: ValueCell.create(1),
    uPlane0: ValueCell.create(Vec4()),
    uPlane1: ValueCell.create(Vec4()),
    uPlane2: ValueCell.create(Vec4()),
    uPlane3: ValueCell.create(Vec4()),
    uPlane4: ValueCell.create(Vec4()),
    uPlane5: ValueCell.create(Vec4()),
    uCameraPlane: ValueCell.create(Vec4()),
    uLevels: ValueCell.create(Mat4()),
    uLevelCount: ValueCell.create(0),
    uActiveCount: ValueCell.create(0),
    uRowThreads: ValueCell.create(0),
    uRunsBase: ValueCell.create(0),
    uMembersBase: ValueCell.create(0),
    uActiveBase: ValueCell.create(0),
    uPoseFirstsBase: ValueCell.create(0),
    uValueBlocksBase: ValueCell.create(0),
    uExtrasBase: ValueCell.create(0),
    uSelectedBase: ValueCell.create(0),
    uCopiesBase: ValueCell.create(0),
    uCoarseScale: ValueCell.create(0),
    uCoarseOrtho: ValueCell.create(0),
    uCoarsePixels: ValueCell.create(0),
    uFrameBase: ValueCell.create(0),
    uFrameRoom: ValueCell.create(0),
    uFrameCursor: ValueCell.create(0),
    uFallback: ValueCell.create(0),
    uFallbackRoom: ValueCell.create(0),
    uDroppedBase: ValueCell.create(0),
    uPositions: ValueCell.create(0),
    uRotations: ValueCell.create(0),
    uPhase: ValueCell.create(0),
    uClipObjectCount: ValueCell.create(0),
    uClipBase: ValueCell.create(0),
    uClipInstance: ValueCell.create(0),
    uGhost: ValueCell.create(0),
    uGhostBase: ValueCell.create(0),
    uMemberCount: ValueCell.create(0),
    uPoseBase: ValueCell.create(0),
    uPoseRoom: ValueCell.create(0),
    uIndexRoom: ValueCell.create(0),
    uRasterRadius: ValueCell.create(0),
    uRasterScale: ValueCell.create(0),
    uRasterNear: ValueCell.create(0),
    uRasterOrtho: ValueCell.create(0),
    uRasterBase: ValueCell.create(0),
    uRasterRowsBase: ValueCell.create(0),
    uMaxReach: ValueCell.create(0),
    uPageBitsBase: ValueCell.create(0),
    uBitsRoom: ValueCell.create(0),
    uRegionCount: ValueCell.create(0),
    uSurvivorsBase: ValueCell.create(0),
    uSurvivorsRoom: ValueCell.create(0),
    uGroupsPass: ValueCell.create(0),
    uOrderNear: ValueCell.create(0),
    uOrderScale: ValueCell.create(0),
    ...hiZValues,
  };
  const composeValues = {
    uRowThreads: ValueCell.create(0),
    uComposeCount: ValueCell.create(0),
    uRunsBase: uniforms.uRunsBase,
    uMembersBase: uniforms.uMembersBase,
    uExtrasBase: uniforms.uExtrasBase,
    uPoseBase: uniforms.uPoseBase,
    uPoseRoom: uniforms.uPoseRoom,
    uValueBlocksBase: uniforms.uValueBlocksBase,
    uPositions: uniforms.uPositions,
    uRotations: uniforms.uRotations,
    tables,
    positions,
    composePages,
    poses,
    groupBoxes,
  };
  const cullValues = {
    ...uniforms,
    tables,
    positions,
    counts,
    regions: bases,
    groupBoxes,
    occludedBits,
    indices,
    poses,
  };
  /** the selected entries' count, then the cull's indirect dispatch (x, y, z) */
  const selectArgs = ValueCell.create(createBuffer((PagedCullGroupArgs + 4) * 4, 'paged cull selected'));
  /** a bit a store page the last cull wanted, a word a copy group lacking a page in view (`paged-cull-wants`) */
  const wanted = ValueCell.create(createBuffer(16, 'paged cull wanted'));
  const lacking = ValueCell.create(createBuffer(16, 'paged cull lacking'));
  const selectValues = {
    uRowThreads: ValueCell.create(0),
    uActiveCount: uniforms.uActiveCount,
    uActiveBase: uniforms.uActiveBase,
    uCopiesBase: uniforms.uCopiesBase,
    uSelectedBase: uniforms.uSelectedBase,
    uCoarseScale: uniforms.uCoarseScale,
    uCoarseOrtho: uniforms.uCoarseOrtho,
    uCoarsePixels: uniforms.uCoarsePixels,
    uCameraPlane: uniforms.uCameraPlane,
    uObjectTransform: uniforms.uObjectTransform,
    uObjectScale: uniforms.uObjectScale,
    uPlane0: uniforms.uPlane0,
    uPlane1: uniforms.uPlane1,
    uPlane2: uniforms.uPlane2,
    uPlane3: uniforms.uPlane3,
    uPlane4: uniforms.uPlane4,
    uPlane5: uniforms.uPlane5,
    uMaxReach: uniforms.uMaxReach,
    uGhost: uniforms.uGhost,
    tables,
    selectArgs,
    uValueBlocksBase: uniforms.uValueBlocksBase,
    lacking,
    uFallback: uniforms.uFallback,
    uDroppedBase: uniforms.uDroppedBase,
    marks: occludedBits,
    uSelectAppend: ValueCell.create(0),
    uCopyMarksBase: ValueCell.create(0),
    ...hiZValues,
    uSelectPass: ValueCell.create(0),
    uOrderNear: uniforms.uOrderNear,
    uOrderScale: uniforms.uOrderScale,
  };
  const selectArgsValues = {
    uSelectedBase: uniforms.uSelectedBase,
    uSelectAppend: selectValues.uSelectAppend,
    uArgsMode: ValueCell.create(0),
    uSurvivorsBase: uniforms.uSurvivorsBase,
    uSurvivorsRoom: uniforms.uSurvivorsRoom,
    survivors: bases,
    tables,
    selectArgs,
  };
  const groupsValues = { ...uniforms, tables, groupBoxes, occludedBits, survivors: bases };
  const wantsValues = {
    ...selectValues,
    uValueBlocksBase: uniforms.uValueBlocksBase,
    uWantPadding: ValueCell.create(PagedLeafPagePadding),
    wanted,
    lacking,
  };
  let wantsKernel: ComputeKernel | undefined;
  let composeKernel: ComputeKernel | undefined;
  let selectKernel: ComputeKernel | undefined;
  let selectArgsKernel: ComputeKernel | undefined;
  let groupsKernel: ComputeKernel | undefined;
  let cullKernel: ComputeKernel | undefined;

  /** the buffer of `c` at least `bytes`, its contents copied where kept (`keep` bytes of them) */
  const ensure = (c: ValueCell<StorageBuffer>, bytes: number, keep = 0) => {
    if (c.ref.value.byteLength >= bytes) return false;
    // grown by half again, within the device's limit
    const grown = Math.min(compute.maxStorageBufferSize, Math.ceil(c.ref.value.byteLength * 1.5));
    const next = createBuffer(Math.max(bytes, grown), 'paged cull');
    if (keep > 0) c.ref.value.copyTo(next, 0, 0, Math.min(keep, c.ref.value.byteLength));
    c.ref.value.destroy();
    ValueCell.update(c, next);
    return true;
  };

  // the tables' layout as written: its sections' capacities (a change of one the tables written whole)
  let tablesPages = -1;
  let tablesRuns = -1;
  let tablesMembers = -1;
  let tablesExtras = -1;
  let tablesActive = -1;
  let tablesCopies = -1;
  /** the active entries (page, copy: two words each) and the copies' poses, as of the last sync */
  let activeEntries = new Uint32Array(2);
  let copyPoses = new Float32Array(0);
  /**
   * the frame's slots for the copies' instances (`frameRoom`): their instances' at most (`frameInstances`), at first
   * `PagedFrameInitialSlots`, grown to what a cull needed (`frameWanted`, its cursor read back: `readCursors`)
   */
  let frameRoom = 0;
  let frameInstances = 0;
  let frameWanted = options.frameInitialSlots ?? PagedFrameInitialSlots;
  /** a read of the cursors pending (`readCursors`) */
  let cursorsRead = false;
  /** the last cull's view had a pyramid to come (`PagedCullView.pyramidToCome`): its overflow not grown for */
  let pyramidToCome = false;
  /** the frame's slots the last cull read back needed past a room at its cap (none: it fit, or the room not capped) */
  let cappedNeed = Number.POSITIVE_INFINITY;
  /** the last cull read back overflowed the frame's slots or a region: the next to fit culled once more (`readCursors`) */
  let lastOverflowed = false;
  const readBack: PagedCullReadBack = {
    reads: 0,
    overflows: 0,
    frameNeed: 0,
    phasesRoom: 0,
    frameRoom: 0,
    frameCapped: false,
    regionsOver: 0,
    regionExcess: 0,
    selected: 0,
    bitsRoom: 0,
    survivorsEntries: 0,
  };
  /**
   * per region the room a cull needed (`readCursors`: half again its count where it overflowed, 0 none asked yet), of
   * the layout's shape `regionShape` (members, levels, raster, ghost: another, none kept); the rooms laid out
   * (`layoutBases`: `regionRooms`, their sum `indexTotal`, the raster regions' `rasterRoom`), `budgeted` where one is
   * below its member's instances (the counts read back)
   */
  let regionWanted = new Float64Array(0);
  let regionShape = '';
  let regionRooms = new Uint32Array(0);
  let roomsVersion = 0;
  let indexTotal = 0;
  let rasterRoom = 0;
  let budgeted = false;
  const regionRoomShare = options.regionRoomShare ?? PagedRegionRoomShare;
  const regionInitialSlots = options.regionInitialSlots ?? PagedRegionInitialSlots;
  /**
   * the selected entries the occlusion bits hold (`uBitsRoom`), grown to half again a cull's count past them (its
   * `selectArgs` count read back: `readSelected`)
   */
  let bitsWanted = options.occlusionInitialEntries ?? PagedOcclusionInitialEntries;
  /**
   * the selected entries whose groups the survivors hold (`PagedCullGroupsKernel`: grown as the bits, by the selected
   * count read back), and those of the last cull (`survivorsEntries`)
   */
  let survivorsWanted = options.survivorsInitialEntries ?? PagedOcclusionInitialEntries;
  let survivorsEntries = 0;
  let bitsRoom = 0;
  let selectedRead = false;
  /** advanced where a cull overflowed the frame's slots or a region's: the view culled again (`PagedCuller.frameVersion`) */
  let frameVersion = 0;
  let positionWords = 2;
  let rotationWords = 0;
  let activeCount = 0;
  let memberCount = 0;
  let instanceTotal = 0;
  let regionCount = 1;
  /** the pages there were at the last sync: what a grown buffer of positions keeps */
  let keptPages = 0;
  let everSynced = false;
  /** the source and the members' reaches of the last sync: one changing nothing of them none (`sync`) */
  let syncedSource: PagedCullSource | undefined;
  let syncedReaches: number[] = [];
  const isUnchanged = (c: PagedCullChanges) =>
    c.written.length === 0 &&
    c.released.length === 0 &&
    (c.copies?.length ?? 0) === 0 &&
    c.runs.length === 0 &&
    c.extras.length === 0 &&
    (c.members?.size ?? 0) === 0 &&
    !c.grown;
  const sameReaches = (src: PagedCullSource, reaches: ArrayLike<number>) => {
    if (syncedReaches.length !== src.pool.count) return false;
    for (let m = 0; m < src.pool.count; ++m) if ((reaches[m] ?? 0) !== syncedReaches[m]) return false;
    return true;
  };
  let culled = false;
  let source: PagedCullSource | undefined;
  /** the pages written and not composed yet */
  const pendingPages = new Set<number>();
  /**
   * the leaf blocks' (a template with coarse levels, `level` 0) pages not written (no value block, no entries: written
   * where a near copy's view reaches them, `updateResidency`), and per leaf page when a cull last wanted it
   */
  const unwrittenPages = new Set<number>();
  const wantedAt = new Map<number, number>();
  /** per template's block (its pages) its largest scale as of the last sync (`collectActive`; -1 not known) */
  const blockScales = new WeakMap<ArrayLike<number>, number>();
  /** per page of a member a particle its extras' largest scale (by its first extra: written again, another) */
  const extraScaleOf = new Map<number, { first: number; scale: number }>();
  const getExtraScale = (src: PagedCullSource, p: number) => {
    const o = p * RecordVec4s * 4;
    const first = src.records[o + 11];
    const known = extraScaleOf.get(p);
    if (known && known.first === first) return known.scale;
    const n = src.records[o + 6];
    const f = new Float32Array(src.extras.buffer, src.extras.byteOffset, src.extras.length);
    let scale = 0;
    for (let j = 0; j < n; ++j) {
      const v = f[(first + j) * 2 + 1];
      if (v > scale) scale = v;
    }
    extraScaleOf.set(p, { first, scale });
    return scale;
  };
  /** the leaf pages the cull wanted and not written yet (read back, `readWants`), written by the next culls */
  const pendingWants = new Set<number>();
  /** the copy groups lacking a leaf page in view as last read back (the occlusion bits grown ahead by them) */
  let lackingCount = 0;
  /** a read of the wants pending (one at a time), and when the last began */
  let wantsRead = false;
  let wantsReadAt = Number.NEGATIVE_INFINITY;
  /**
   * the wants read back so far, and per leaf page the last read that wanted it: a page let go where reads since did
   * not want it, not by the time alone (a still view's culls none, its pages' wants not read: kept)
   */
  let wantsReads = 0;
  const wantedRead = new Map<number, number>();
  /** the copy groups (a template's copies of one top particle across its blocks), as of the last sync */
  let groupCount = 0;
  /** a template has coarse levels (the fallback's share of the frame's slots: none without, `PagedFallbackShare`) */
  let coarseCopies = false;
  /** a sync to run though nothing changed (leaf pages written or let go: the value blocks' table) */
  let forceSync = false;
  /** the last cull wrote leaf pages: its counts not growing the rooms (`cull`) */
  let residencyChanged = false;
  const isLeafBlock = (t: NonNullable<PagedCullSource['templates']>[number]) =>
    (t.copies.level ?? 0) === 0 && (t.copies.levels ?? 0) > 0 && t.pages.length > 0;
  /** the members' instances as of the last sync, and that sync's: the regions laid out again where another */
  let memberCounts = new Int32Array(0);
  let syncVersion = 0;
  let basesKey = '';

  // the pose slots: a range a page (its count's), its first in the tables' section (`NoSlots`: none)
  const poseBytes = PagedPoseWords * 4;
  /** the slots the device's limit holds after the front: the persistent ones' and the frame's */
  const limitSlots = () => Math.floor((compute.maxStorageBufferSize - posesFront.byteLength) / poseBytes);
  // the frame's room kept one at least: a frame's room grown later takes what the persistent leave (`sizeFrame`)
  const maxPoseSlots = () => Math.max(1, limitSlots() - Math.max(1, frameRoom));
  const slots = createRangeAllocator(0, 1);
  let pageFirsts = new Uint32Array(0);
  /**
   * per page its value block in the values' buffer (`NoSlots` none), the blocks let go (reused), the blocks made so
   * far: the buffer holds the written pages alone. Why: by the store's capacity it held a template's every page
   * (the yeast cell's leaves 576 MB) while a view drew its coarse levels alone
   */
  let pageBlocks = new Uint32Array(0);
  const freeBlocks: number[] = [];
  let blockCount = 0;
  let pageSlotCounts = new Int32Array(0);
  /** the poses' buffer for `capacity` slots after the front, what was there copied (the front's poses moved to `toFront`) */
  const sizePoses = (capacity: number, fromFront: number, toFront: number) => {
    // the frame's slots after the persistent ones (not kept: written each cull)
    const bytes = toFront + (capacity + frameRoom) * poseBytes;
    if (poses.ref.value.byteLength >= bytes && fromFront === toFront) return;
    const kept = slots.capacity * poseBytes;
    const next = createBuffer(Math.max(bytes, toFront + 16), 'paged cull poses');
    if (kept > 0)
      poses.ref.value.copyTo(next, fromFront, toFront, Math.min(kept, poses.ref.value.byteLength - fromFront));
    poses.ref.value.destroy();
    ValueCell.update(poses, next);
    // the copy before the front and the composes written after
    compute.submit();
    writeFront();
  };
  /** page `p`'s slots for `count` particles (its old ones let go): grown where full, within the device's limit */
  const allocateSlots = (p: number, count: number) => {
    if (pageFirsts.length <= p) {
      const firsts = new Uint32Array(Math.max(p + 1, pageFirsts.length * 2)).fill(NoSlots);
      firsts.set(pageFirsts);
      pageFirsts = firsts;
      const counts = new Int32Array(firsts.length);
      counts.set(pageSlotCounts);
      pageSlotCounts = counts;
    }
    freeSlots(p);
    if (count === 0) return;
    let first = slots.allocate(count);
    if (first < 0) {
      const max = maxPoseSlots();
      const capacity = Math.min(max, Math.max(slots.capacity + count, Math.ceil(slots.capacity * 1.5)));
      if (capacity > slots.capacity) {
        const front = posesFront.byteLength;
        sizePoses(capacity, front, front);
        slots.grow(capacity);
        ValueCell.updateIfChanged(uniforms.uPoseRoom, capacity);
        ValueCell.updateIfChanged(uniforms.uFrameBase, capacity);
        first = slots.allocate(count);
      }
    }
    if (first < 0) return;
    pageFirsts[p] = first;
    pageSlotCounts[p] = count;
  };
  function freeSlots(p: number) {
    if (p >= pageFirsts.length || pageFirsts[p] === NoSlots) return;
    slots.free(pageFirsts[p], pageSlotCounts[p]);
    pageFirsts[p] = NoSlots;
    pageSlotCounts[p] = 0;
  }
  const writePoseFirsts = () => {
    const n = Math.max(1, tablesPages);
    const words = new Uint32Array(Math.ceil(n / 4) * 4).fill(NoSlots);
    words.set(pageFirsts.subarray(0, Math.min(pageFirsts.length, n)));
    tables.ref.value.write(words, uniforms.uPoseFirstsBase.ref.value * 16);
    words.fill(NoSlots);
    words.set(pageBlocks.subarray(0, Math.min(pageBlocks.length, n)));
    tables.ref.value.write(words, uniforms.uValueBlocksBase.ref.value * 16);
  };

  /** the members' rows as written (center, reach), and whether the next are written whole (the tables made again) */
  let writtenMembers = new Float32Array(0);
  let membersWhole = true;
  /** the largest scale of the runs and extras synced (1: an instance's scale where none) */
  let maxScale = 1;
  let scannedScales = false;
  const writeMembers = (src: PagedCullSource, reaches: ArrayLike<number>, changed?: ReadonlySet<number>) => {
    const n = src.pool.count;
    const { centers } = src.pool;
    // the members that changed alone, their rows' span written (a step: a few of thousands), where the rows are as
    // many as last written
    if (changed && !membersWhole && writtenMembers.length === Math.max(1, n) * 4) {
      const rows = writtenMembers;
      let lo = n;
      let hi = -1;
      for (const m of changed) {
        if (m < 0 || m >= n) continue;
        const o = m * 4;
        const x = centers[m * 3];
        const y = centers[m * 3 + 1];
        const z = centers[m * 3 + 2];
        const w = reaches[m] ?? 0;
        if (rows[o] === Math.fround(x) && rows[o + 1] === Math.fround(y) && rows[o + 2] === Math.fround(z)) {
          if (rows[o + 3] === Math.fround(w)) continue;
        }
        rows[o] = x;
        rows[o + 1] = y;
        rows[o + 2] = z;
        rows[o + 3] = w;
        if (m < lo) lo = m;
        if (m > hi) hi = m;
      }
      if (hi >= lo)
        tables.ref.value.write(rows.subarray(lo * 4, (hi + 1) * 4), uniforms.uMembersBase.ref.value * 16 + lo * 16);
      return;
    }
    const data = new Float32Array(Math.max(1, n) * 4);
    for (let m = 0; m < n; ++m) {
      data[m * 4] = centers[m * 3];
      data[m * 4 + 1] = centers[m * 3 + 1];
      data[m * 4 + 2] = centers[m * 3 + 2];
      data[m * 4 + 3] = reaches[m] ?? 0;
    }
    // the rows that changed alone (a member's center and reach: set where it is given a target, its pages new: their
    // poses composed with it)
    const base = uniforms.uMembersBase.ref.value * 16;
    const before = writtenMembers;
    if (before.length !== data.length || membersWhole) {
      tables.ref.value.write(data, base);
    } else {
      for (let m = 0; m < n; ++m) {
        const o = m * 4;
        if (
          data[o] === before[o] &&
          data[o + 1] === before[o + 1] &&
          data[o + 2] === before[o + 2] &&
          data[o + 3] === before[o + 3]
        )
          continue;
        let e = m + 1;
        while (e < n) {
          const q = e * 4;
          if (
            data[q] === before[q] &&
            data[q + 1] === before[q + 1] &&
            data[q + 2] === before[q + 2] &&
            data[q + 3] === before[q + 3]
          )
            break;
          ++e;
        }
        tables.ref.value.write(data.subarray(o, e * 4), base + o * 4);
        m = e;
      }
    }
    writtenMembers = data;
    membersWhole = false;
  };

  let clip: PagedCullClip | undefined;
  const writeClip = () => {
    if (clip && clip.count > 0 && tablesPages >= 0) {
      const n = Math.min(clip.count, Clip.MaxObjectCount);
      tables.ref.value.write(clip.data.subarray(0, n * ClipStorageVec4s * 4), uniforms.uClipBase.ref.value * 16);
    }
  };

  let rasterRows: Float32Array | undefined;
  const writeRasterRows = () => {
    if (!rasterRows || tablesPages < 0) return;
    const n = Math.min(rasterRows.length / 4, tablesMembers);
    if (n > 0) tables.ref.value.write(rasterRows.subarray(0, n * 4), uniforms.uRasterRowsBase.ref.value * 16);
  };
  let posesFront: Uint32Array = new Uint32Array(0);
  const writeFront = () => {
    if (posesFront.length > 0) poses.ref.value.write(posesFront, 0);
  };

  /**
   * The active entries of `src` (`activeEntries`, two words each: a page, its copy or `NoCopy`) and its copies' poses
   * (`copyPoses`), the frame's room their instances: a page drawn as it is once, a template's (`PagedTemplatePage`,
   * in its frame) once a copy of it
   */
  const collectActive = (src: PagedCullSource) => {
    const entries: number[] = [];
    const recordWords = RecordVec4s * 4;
    const runScales = new Float32Array(src.runs.buffer, src.runs.byteOffset, src.runs.length);
    for (let p = 0; p < src.pageCount; ++p) {
      const o = p * recordWords;
      if (src.records[o + 6] > 0 && !(src.records[o + PagedFlagsWord] & PagedTemplatePage)) entries.push(p, NoCopy);
    }
    let copies = 0;
    let instances = 0;
    for (const t of src.templates ?? []) copies += t.copies.poses.length / 8;
    copyPoses = new Float32Array(copies * CopyVec4s * 4);
    let c = 0;
    const q = Quat();
    const center = Vec3();
    // the copy groups: a template's copies of one top particle across its blocks (their poses the same), what a copy
    // lacking a leaf page holds coarse in all of them (`copyGroup`)
    const groupBase = new Map<Float32Array, number>();
    groupCount = 0;
    coarseCopies = (src.templates ?? []).some((t) => (t.copies.levels ?? 0) > 0);
    for (const t of src.templates ?? []) {
      if (groupBase.has(t.copies.poses)) continue;
      groupBase.set(t.copies.poses, groupCount);
      groupCount += t.copies.poses.length / 8;
    }
    for (const t of src.templates ?? []) {
      const { poses, level = 0, levels = 0, spacing = 0, radius = -1 } = t.copies;
      const group = groupBase.get(poses)!;
      // its block's largest scale: its pages' runs', a page of a member a particle's its extras' (a tile-ordered
      // template's: thousands of targets a page), once a block (`extraScales`); -1 where not known
      let blockScale = 0;
      for (let i = 0; i < t.pages.length && blockScale >= 0; ++i) {
        const o = t.pages[i] * recordWords;
        const first = src.records[o + 4];
        const count = src.records[o + 5];
        if (count === PagedPerParticleMembers) {
          blockScale = Math.max(blockScale, getExtraScale(src, t.pages[i]));
          continue;
        }
        for (let k = first; k < first + count && blockScale >= 0; ++k) {
          const v = runScales[k * 4 + 3];
          blockScale = v > 0 ? Math.max(blockScale, v) : -1;
        }
      }
      const n = poses.length / 8;
      const tag =
        (Math.min(level, PagedMaxCoarseLevels) << PagedEntryLevelShift) |
        (Math.min(levels, PagedMaxCoarseLevels) << PagedEntryLevelsShift);
      for (let k = 0; k < n; ++k) {
        // its pose; its center the template's placed (`t + R c`), its template's spacing
        const o = (c + k) * CopyVec4s * 4;
        copyPoses.set(poses.subarray(k * 8, k * 8 + 8), o);
        Vec3.set(center, t.copies.center?.[0] ?? 0, t.copies.center?.[1] ?? 0, t.copies.center?.[2] ?? 0);
        Quat.set(q, poses[k * 8 + 4], poses[k * 8 + 5], poses[k * 8 + 6], poses[k * 8 + 7]);
        Vec3.transformQuat(center, center, q);
        copyPoses[o + 8] = poses[k * 8] + center[0];
        copyPoses[o + 9] = poses[k * 8 + 1] + center[1];
        copyPoses[o + 10] = poses[k * 8 + 2] + center[2];
        copyPoses[o + 11] = spacing;
        copyPoses[o + 12] = radius;
        copyPoses[o + 13] = blockScale;
        blockScales.set(t.pages, blockScale);
        copyPoses[o + 14] = 0;
        copyPoses[o + 15] = group + k;
        for (let i = 0; i < t.pages.length; ++i) {
          const p = t.pages[i];
          // the frame's instances at most: every copy the template whole (a level's fewer), resident or not (a copy
          // drawn at a coarse level while its leaf block is not)
          if (level === 0) instances += src.records[p * recordWords + 6];
          // a leaf page not written too: tested by the wants, not drawn (`paged-cull-select`)
          entries.push((p | tag) >>> 0, c + k);
        }
      }
      c += n;
    }
    // the copies' centers' box: its sphere (the select's order, `setOrderRange`)
    const lo = Vec3.create(Infinity, Infinity, Infinity);
    const hi = Vec3.create(-Infinity, -Infinity, -Infinity);
    for (let k = 0; k < copies; ++k) {
      const o = k * CopyVec4s * 4 + 8;
      for (let i = 0; i < 3; ++i) {
        lo[i] = Math.min(lo[i], copyPoses[o + i]);
        hi[i] = Math.max(hi[i], copyPoses[o + i]);
      }
    }
    if (copies > 0) {
      Vec3.scale(copyDepthSphere.center, Vec3.add(copyDepthSphere.center, lo, hi), 0.5);
      copyDepthSphere.radius = Vec3.distance(lo, hi) / 2;
    } else {
      Vec3.set(copyDepthSphere.center, 0, 0, 0);
      copyDepthSphere.radius = 0;
    }
    activeEntries = Uint32Array.from(entries.length ? entries : [0, NoCopy]);
    activeCount = entries.length / 2;
    return { copies, instances };
  };
  /**
   * The frame's slots sized: what a cull needed (`frameWanted`, `PagedFrameInitialSlots` at first), the copies'
   * instances at most, within `PagedFrameMaxSlots` and what the device's limit leaves after the persistent slots (a
   * zoom into the yeast cell's copies: its template's 36 M slots, when a template's particles had slots of their own,
   * and a frame's 32 M past 2 GiB); the poses' buffer
   * grown for them (the frame's slots not kept: a cull writes them)
   */
  const sizeFrame = () => {
    const fit = Math.max(1, limitSlots() - slots.capacity);
    const room = Math.min(frameInstances, frameWanted, PagedFrameMaxSlots, fit);
    // a room below the copies' instances: a share of it the fallback's (`PagedFallbackShare`)
    const fallback = coarseCopies && room < frameInstances ? Math.floor(room * PagedFallbackShare) : 0;
    ValueCell.updateIfChanged(uniforms.uFrameRoom, room - fallback);
    ValueCell.updateIfChanged(uniforms.uFallbackRoom, fallback);
    if (room === frameRoom) return;
    frameRoom = room;
    const front = posesFront.byteLength;
    sizePoses(slots.capacity, front, front);
  };

  /**
   * The cursors of the cull just run read back (its last phase's `counts`: the regions' and the frame's, the first's
   * and the second's slots; in the second phase the first's region counts too, in `bases`), one read at a time, where
   * the frame's slots or a region's are fewer than could be needed. Past a room instances were not drawn: the view
   * culled again (`frameVersion`) and a draw asked for (a settled promise tracked as pending work); the room grown to
   * half again what it needed (`grow`), not where a pyramid is to come (`PagedCullView.pyramidToCome`: an app's first
   * frames, nothing occluded, every instance in view needed; the next cull's pyramid hides most of them: the isg cell
   * on a 4x4x4 grid 29.4 M without, 0.56 M with), nor where a cull of a pyramid (`tested`) selected entries past the
   * occlusion bits' room (not tested against it: `readSelected` grows it; a cull without one tests none)
   */
  const readCursors = (phaseCounts: StorageBuffer, grow: boolean, second: boolean, tested: boolean) => {
    const frameBudgeted = frameRoom > 0 && frameRoom < frameInstances;
    if (cursorsRead || !(frameBudgeted || budgeted)) return;
    cursorsRead = true;
    const n = regionCount;
    const rooms = regionRooms;
    const bits = bitsRoom;
    const phasesRoom = uniforms.uFrameRoom.ref.value;
    compute.submit();
    const reads = [phaseCounts.readAsync(0, n * 4 + 4), selectArgs.ref.value.readAsync(0, 4)];
    if (second && budgeted) reads.push(bases.ref.value.readAsync(n * 4, n * 4));
    Promise.all(reads)
      .then(([data, selectedData, priorData]) => {
        cursorsRead = false;
        const c = new Uint32Array(data);
        const prior = priorData ? new Uint32Array(priorData) : undefined;
        // every selected entry tested against the pyramid, else what it needed not what the view will
        if (tested && new Uint32Array(selectedData)[0] > bits) grow = false;
        let overflowed = false;
        // a room at its cap (`sizeFrame`) not grown by a cull again: overflowing it asks for another while each needs
        // fewer than the one before (its pyramid this view's: a zoom's first cull of the last view's pyramid, the isg
        // cell's 16x16x16 grid 87.6 M of 33.5 M, the rest not drawn till another cull), none past that (a still
        // view's every frame culled else)
        const frameCapped = frameRoom >= Math.min(PagedFrameMaxSlots, Math.max(1, limitSlots() - slots.capacity));
        readBack.reads += 1;
        readBack.frameNeed = c[n];
        readBack.phasesRoom = phasesRoom;
        readBack.frameRoom = frameRoom;
        readBack.frameCapped = frameCapped;
        readBack.regionsOver = 0;
        readBack.regionExcess = 0;
        for (let r = 0; r < n; ++r) {
          const excess = c[r] + (prior ? prior[r] : 0) - rooms[r];
          if (excess <= 0) continue;
          readBack.regionsOver += 1;
          readBack.regionExcess = Math.max(readBack.regionExcess, excess);
        }
        // the phases' room (the fallback's share after it)
        if (frameBudgeted && c[n] > phasesRoom) {
          if (!frameCapped) {
            overflowed = true;
            if (grow) frameWanted = Math.min(PagedFrameMaxSlots, Math.ceil(c[n] * 1.5));
          } else if (c[n] < cappedNeed) {
            overflowed = true;
          }
          cappedNeed = frameCapped ? c[n] : Number.POSITIVE_INFINITY;
        } else {
          cappedNeed = Number.POSITIVE_INFINITY;
        }
        // the regions as laid out when culled (another layout since: its rooms of its own)
        if (budgeted && rooms === regionRooms) {
          let grown = false;
          for (let r = 0; r < n; ++r) {
            const needed = c[r] + (prior ? prior[r] : 0);
            if (needed <= rooms[r]) continue;
            overflowed = true;
            if (!grow) continue;
            regionWanted[r] = Math.max(regionWanted[r], Math.ceil(needed * 1.5));
            grown = true;
          }
          if (grown) roomsVersion += 1;
        }
        if (!overflowed) {
          // the first to fit after one that did not: its pyramid of a frame with instances missing (holes, less
          // occluded); culled once more with one of a whole frame (a still view's lists otherwise kept)
          if (!lastOverflowed) return;
          lastOverflowed = false;
        } else {
          lastOverflowed = true;
          readBack.overflows += 1;
        }
        frameVersion += 1;
        void resources.pending.track(Promise.resolve(), 'paged cull: slots overflowed');
      })
      .catch(() => {
        cursorsRead = false;
      });
  };

  /**
   * The extras' section grown to `capacity` (the tables' last: no other section moved, none written again): the
   * tables grown by a copy of what was there, submitted before what is written now (a write goes ahead of the commands
   * not submitted: the copy would overwrite it); the extras the store added written as its changes
   */
  const growExtras = (capacity: number) => {
    const extrasBase = uniforms.uExtrasBase.ref.value;
    const kept = (extrasBase + Math.ceil((tablesExtras * 2) / 4)) * 16;
    tablesExtras = capacity;
    if (ensure(tables, (extrasBase + Math.ceil((capacity * 2) / 4)) * 16, kept)) compute.submit();
  };
  const writeActive = () => {
    tables.ref.value.write(activeEntries, uniforms.uActiveBase.ref.value * 16);
    if (copyPoses.length > 0) tables.ref.value.write(copyPoses, uniforms.uCopiesBase.ref.value * 16);
  };

  /**
   * The pages pending composed (their values as written by store index, kept): a pass over their particles, at a cull
   * alone (the page list a queue write: one compose between submits, a frame's)
   */
  const compose = () => {
    if (pendingPages.size === 0 || !source?.format) return;
    composeKernel ??= compute.createKernel(PagedCullComposeKernel, PagedCullComposeSchema, composeValues);
    const list: number[] = [];
    // a page without pose slots too (a template's, one past the device's limit): its group boxes; a page written (its
    // value block: its values on the CPU maybe given up, `releaseValues`)
    for (const p of pendingPages) {
      if (p < pageBlocks.length && pageBlocks[p] !== NoSlots)
        list.push(p, p < pageFirsts.length ? pageFirsts[p] : NoSlots);
    }
    pendingPages.clear();
    if (list.length === 0) return;
    const count = list.length / 2;
    ensure(composePages, list.length * 4);
    composePages.ref.value.write(Uint32Array.from(list), 0);
    const groups = count * PagedCullGroupsPerPage;
    const { x, y } = getDispatchSize(groups, 65535);
    ValueCell.updateIfChanged(composeValues.uComposeCount, count);
    ValueCell.updateIfChanged(composeValues.uRowThreads, x * PagedCullWorkgroupKeys);
    compute.beginComputePass({ label: 'paged cull compose' });
    composeKernel.dispatch(x, y);
    compute.endComputePass();
  };

  /**
   * The regions' rooms and firsts, one after another in region order: a level's region (`member * levels + level`), a
   * raster region's (from `levels * members`, level-major; its first its end, its slots down from it), the ghost's
   * after them. A region's room what a cull needed (`regionWanted`), at first its member's instances times
   * `regionRoomShare` (at least `PagedRegionMinRoom`; the shares within `PagedRegionInitialSlots` in all), at most
   * its member's instances (the ghost's too). Laid
   * out where the store's members' instances, the regions or a room changed
   */
  const layoutBases = (levels: number, raster: boolean, ghost: boolean) => {
    const n = memberCount;
    const shape = `${n},${levels},${raster},${ghost}`;
    if (shape !== regionShape) {
      regionShape = shape;
      regionWanted = new Float64Array(regionCount);
    }
    const key = `${syncVersion},${shape},${roomsVersion}`;
    if (key === basesKey) return;
    basesKey = key;
    const count = Math.max(1, regionCount);
    const rasterBase = n * levels;
    const ghostBase = rasterBase + (raster ? levels * n : 0);
    const rooms = new Uint32Array(count);
    const firsts = new Uint32Array(count);
    budgeted = false;
    rasterRoom = 0;
    const memberOf = (r: number) =>
      r < rasterBase ? Math.floor(r / levels) : r < ghostBase ? (r - rasterBase) % n : r - ghostBase;
    // the shares scaled down where their sum passes the first rooms' budget (not where every instance has room)
    let shares = 0;
    for (let r = 0; r < regionCount; ++r) shares += memberCounts[memberOf(r)] * regionRoomShare;
    const share =
      regionRoomShare < 1 && shares > regionInitialSlots
        ? (regionRoomShare * regionInitialSlots) / shares
        : regionRoomShare;
    let sum = 0;
    for (let r = 0; r < regionCount; ++r) {
      const m = memberOf(r);
      const all = memberCounts[m];
      const start = Math.min(all, Math.max(PagedRegionMinRoom, Math.ceil(all * share)));
      const room = Math.min(all, Math.max(start, regionWanted[r]));
      rooms[r] = room;
      if (room < all) budgeted = true;
      const isRaster = r >= rasterBase && r < ghostBase;
      if (isRaster) rasterRoom += room;
      // a first past 2^32 slots: its region none written (the indices' room far below)
      firsts[r] = Math.min(0xffffffff, isRaster ? sum + room : sum);
      sum += room;
    }
    regionRooms = rooms;
    indexTotal = sum;
    // the first phase's counts after the firsts (the second's), the rooms after them
    ensure(bases, count * 12);
    bases.ref.value.write(firsts, 0);
    bases.ref.value.write(rooms, count * 8);
  };

  /** the phase's pass: its cursors cleared, its particles culled into the regions */
  let dispatchSize = [1, 1];
  let culledWithPyramid = false;
  const runPass = (label: string, frameCursor: boolean) => {
    // the first phase's frame cursors from 0 (the phases' and the fallback's); the second's the first's (copied there:
    // its frame slots past them)
    counts.ref.value.clear(0, regionCount * 4 + (frameCursor ? 8 : 0));
    cullGroups(label);
  };
  /**
   * The selected entries' groups tested whole (`PagedCullGroupsKernel`, a thread a group: the survivors), then the cull
   * a workgroup a survivor, by the uniforms as set
   */
  const cullGroups = (label: string) => {
    groupsKernel ??= compute.createKernel(PagedCullGroupsKernel, PagedCullGroupsSchema, groupsValues);
    // the count and the buckets' cursors cleared; counted, their firsts and the cull's dispatch, placed
    bases.ref.value.clear(uniforms.uSurvivorsBase.ref.value * 4, (1 + PagedCullOrderBuckets) * 4);
    compute.beginComputePass({ label: `${label}, groups` });
    ValueCell.updateIfChanged(uniforms.uGroupsPass, 0);
    groupsKernel.dispatchIndirect(selectArgs.ref.value, 4);
    ValueCell.updateIfChanged(selectArgsValues.uArgsMode, 1);
    selectArgsKernel!.dispatch(1);
    ValueCell.updateIfChanged(selectArgsValues.uArgsMode, 0);
    ValueCell.updateIfChanged(uniforms.uGroupsPass, 1);
    groupsKernel.dispatchIndirect(selectArgs.ref.value, 4);
    compute.endComputePass();
    compute.beginComputePass({ label });
    cullKernel!.dispatchIndirect(selectArgs.ref.value, (PagedCullGroupArgs + 1) * 4);
    compute.endComputePass();
  };
  /**
   * The selected count of the cull just run read back (one read at a time), where the occlusion bits hold fewer than
   * the active entries: past their room the bits grown to half again the count, the view culled again
   * (`frameVersion`: its entries past the room were not tested against the pyramid, their groups past the survivors'
   * room not culled)
   */
  const readSelected = (withPyramid: boolean) => {
    const bitsShort = withPyramid && bitsRoom < activeCount;
    if (selectedRead || (!bitsShort && survivorsEntries >= activeCount)) return;
    selectedRead = true;
    const room = bitsRoom;
    const entries = survivorsEntries;
    compute.submit();
    selectArgs.ref.value
      .readAsync(0, 4)
      .then((data) => {
        selectedRead = false;
        const count = new Uint32Array(data)[0];
        readBack.selected = count;
        readBack.bitsRoom = room;
        readBack.survivorsEntries = entries;
        let grown = false;
        if (bitsShort && count > room) {
          bitsWanted = Math.max(bitsWanted, Math.ceil(count * 1.5));
          grown = true;
        }
        if (count > entries) {
          survivorsWanted = Math.max(survivorsWanted, Math.ceil(count * 1.5));
          grown = true;
        }
        if (!grown) return;
        frameVersion += 1;
        void resources.pending.track(Promise.resolve(), 'paged cull: occlusion bits or survivors grown');
      })
      .catch(() => {
        selectedRead = false;
      });
  };
  /** the entries the phases cull (`PagedCullSelectKernel`), by the view's uniforms as set */
  const select = (append = false) => {
    selectKernel ??= compute.createKernel(PagedCullSelectKernel, PagedCullSelectSchema, selectValues);
    selectArgsKernel ??= compute.createKernel(PagedCullSelectArgsKernel, PagedCullSelectArgsSchema, selectArgsValues);
    const [x, y] = dispatchSize;
    selectArgs.ref.value.clear(16, PagedCullOrderBuckets * 4);
    ValueCell.updateIfChanged(selectValues.uSelectAppend, append ? 1 : 0);
    compute.beginComputePass({ label: append ? 'paged cull select, second phase' : 'paged cull select' });
    ValueCell.updateIfChanged(selectValues.uSelectPass, 0);
    selectKernel.dispatch(x, y);
    selectArgsKernel.dispatch(1);
    ValueCell.updateIfChanged(selectValues.uSelectPass, 1);
    selectKernel.dispatch(x, y);
    compute.endComputePass();
  };
  /**
   * the copies' centers' sphere (`sync`), what the select's order spans in a view (`PagedCullOrderBuckets`): its
   * depth's range, moved by the object's transient transform
   */
  const copyDepthSphere = { center: Vec3(), radius: 0 };
  const _orderCenter = Vec3();
  const setOrderRange = (view: PagedCullView, objectTransform: Mat4) => {
    const center = Vec3.transformMat4(_orderCenter, copyDepthSphere.center, objectTransform);
    const radius = copyDepthSphere.radius * Mat4.getMaxScaleOnAxis(objectTransform);
    const d = Vec3.dot(view.cameraPlane.normal, center) + view.cameraPlane.constant;
    ValueCell.updateIfChanged(uniforms.uOrderNear, d - radius);
    ValueCell.updateIfChanged(uniforms.uOrderScale, radius > 0 ? PagedCullOrderBuckets / (2 * radius) : 0);
  };

  /**
   * The fallback after the last phase (`PagedFallbackShare`, where the frame's room has a share for it): the dropped
   * copy groups' entries a level coarser selected, culled into the last phase's counts (its regions' slots after its)
   * without the pyramid, the frame's slots from the fallback's. A select and an indirect dispatch a cull: none dropped
   * none culled. Not by the pyramid: one primed from an occluder grid (`HiZPass.setOccluder`, the particles' densities)
   * hides what is behind where the copies would be, the dropped ones' own instances too (the isg cell's 16x16x16 grid a
   * zoom step's first frame, every copy a level coarser: 11.1 M of its 14.3 M hidden, the holes kept); a copy's sphere
   * by it hid few
   */
  const fallback = () => {
    if (uniforms.uFallbackRoom.ref.value === 0 || groupCount === 0) return;
    ValueCell.updateIfChanged(uniforms.uFallback, 1);
    select();
    const hiZ = uniforms.uHiZ.ref.value;
    ValueCell.updateIfChanged(uniforms.uHiZ, false);
    cullGroups('paged cull fallback');
    ValueCell.updateIfChanged(uniforms.uHiZ, hiZ);
    ValueCell.updateIfChanged(uniforms.uFallback, 0);
  };

  /** the leaf pages the view wants (`PagedCullWantsKernel`) into `wanted` and `lacking`, cleared first */
  const wants = () => {
    const pageWords = Math.max(1, Math.ceil(Math.max(1, source?.pageCapacity ?? 1) / 32));
    if (ensure(wanted, pageWords * 4)) ValueCell.update(wantsValues.wanted, wanted.ref.value);
    ensure(lacking, Math.max(1, groupCount) * 4);
    wanted.ref.value.clear(0, pageWords * 4);
    lacking.ref.value.clear(0, Math.max(1, groupCount) * 4);
    if (!source?.templates?.some(isLeafBlock)) return;
    wantsKernel ??= compute.createKernel(PagedCullWantsKernel, PagedCullWantsSchema, wantsValues);
    const [x, y] = dispatchSize;
    compute.beginComputePass({ label: 'paged cull wants' });
    wantsKernel.dispatch(x, y);
    compute.endComputePass();
  };
  /**
   * The wants of the cull just run read back (one read at a time): a wanted page's time and read (`wantedAt`,
   * `wantedRead`: what its release goes by), the unwritten ones written by the next culls (`pendingWants`: the view
   * culled again)
   */
  const readWants = () => {
    if (wantsRead || !source?.templates?.some(isLeafBlock)) return;
    // every cull while wanted pages are not written yet, else a few times a second (the written ones' wants kept up
    // for their release, a page wanted anew found): not a submit and a read a frame
    const now = performance.now();
    if (pendingWants.size === 0 && now - wantsReadAt < (options.wantsReadMs ?? PagedWantsReadMs)) return;
    wantsReadAt = now;
    wantsRead = true;
    const pageWords = Math.max(1, Math.ceil(Math.max(1, source.pageCapacity) / 32));
    const groups = Math.max(1, groupCount);
    compute.submit();
    Promise.all([wanted.ref.value.readAsync(0, pageWords * 4), lacking.ref.value.readAsync(0, groups * 4)])
      .then(([w, l]) => {
        wantsRead = false;
        const bits = new Uint32Array(w);
        const at = performance.now();
        let added = false;
        for (let k = 0; k < bits.length; ++k) {
          let b = bits[k];
          while (b !== 0) {
            const low = b & -b;
            const p = k * 32 + 31 - Math.clz32(low);
            b ^= low;
            wantedAt.set(p, at);
            wantedRead.set(p, wantsReads + 1);
            if (unwrittenPages.has(p) && !pendingWants.has(p)) {
              pendingWants.add(p);
              added = true;
            }
          }
        }
        wantsReads += 1;
        const lack = new Uint32Array(l);
        lackingCount = 0;
        for (let k = 0; k < lack.length; ++k) if (lack[k] !== 0) lackingCount += 1;
        if (added) {
          frameVersion += 1;
          void resources.pending.track(Promise.resolve(), 'paged cull: leaf pages wanted');
        }
      })
      .catch(() => {
        wantsRead = false;
      });
  };
  const leafUploadBytes = options.leafUploadBytes ?? PagedLeafUploadBytes;
  /** a page's bytes in the values' buffer: its positions, then its rotations */
  const valueBytes = () => PagedCullPageSize * (positionWords + rotationWords) * 4;
  /** page `p`'s values as stored written into its block */
  /** page `p`'s value block, made where it has none (the buffer grown by half again, what it held copied) */
  const blockOf = (p: number) => {
    if (pageBlocks.length <= p) {
      const blocks = new Uint32Array(Math.max(p + 1, pageBlocks.length * 2)).fill(NoSlots);
      blocks.set(pageBlocks);
      pageBlocks = blocks;
    }
    if (pageBlocks[p] !== NoSlots) return pageBlocks[p];
    const b = freeBlocks.pop() ?? blockCount++;
    pageBlocks[p] = b;
    // a copy of what was there submitted before what is written now (a write goes ahead of the commands not submitted)
    if (ensure(positions, (b + 1) * valueBytes(), b * valueBytes())) compute.submit();
    return b;
  };
  /** room for `n` more value blocks (past the ones let go) in the values' buffer, what it held copied */
  const reserveBlocks = (n: number) => {
    const blocks = blockCount + Math.max(0, n - freeBlocks.length);
    if (ensure(positions, blocks * valueBytes(), blockCount * valueBytes())) compute.submit();
  };
  /** page `p`'s value block let go (reused by another page) */
  const freeBlock = (p: number) => {
    if (p >= pageBlocks.length || pageBlocks[p] === NoSlots) return;
    freeBlocks.push(pageBlocks[p]);
    pageBlocks[p] = NoSlots;
  };
  /** page `p`'s values written into its block (copied: the source's given up then, `releaseValues`) */
  const writeValues = (p: number, values: NonNullable<PagedCullSource['values'][number]>) => {
    const at = blockOf(p) * valueBytes();
    positions.ref.value.write(values.positions, at);
    if (rotationWords > 0 && values.rotations) {
      positions.ref.value.write(values.rotations, at + PagedCullPageSize * positionWords * 4);
    }
    source?.releaseValues?.(p);
  };
  /** the value blocks there is room for (`PagedValueBytes`, the device's limit) */
  const blocksBudget = () =>
    Math.floor(Math.min(options.valueBytes ?? PagedValueBytes, compute.maxStorageBufferSize) / valueBytes());
  /** written leaf page `p` let go (its value block, its slots): wanted again, written again */
  const releaseLeaf = (p: number) => {
    freeBlock(p);
    freeSlots(p);
    pendingPages.delete(p);
    unwrittenPages.add(p);
  };
  /**
   * up to `n` written leaf pages the last read of the wants did not want let go for room (the least recently wanted
   * first): how many
   */
  const releaseForRoom = (src: PagedCullSource, n: number) => {
    const unwanted: number[] = [];
    for (const t of src.templates ?? []) {
      if (!isLeafBlock(t)) continue;
      for (let i = 0; i < t.pages.length; ++i) {
        const p = t.pages[i];
        if (unwrittenPages.has(p) || p >= pageBlocks.length || pageBlocks[p] === NoSlots) continue;
        if ((wantedRead.get(p) ?? 0) < wantsReads) unwanted.push(p);
      }
    }
    const when = (p: number) => wantedAt.get(p) ?? Number.NEGATIVE_INFINITY;
    unwanted.sort((a, b) => when(a) - when(b));
    const count = Math.min(n, unwanted.length);
    for (let i = 0; i < count; ++i) releaseLeaf(unwanted[i]);
    return count;
  };
  /** page `p`'s pose slots: its particles', none for a template's (its copies' composed by the cull) */
  const poseCount = (src: PagedCullSource, p: number) => {
    const o = p * RecordVec4s * 4;
    return src.records[o + PagedFlagsWord] & PagedTemplatePage ? 0 : src.records[o + 6];
  };
  /**
   * The leaf pages the cull wanted written (`pendingWants`, within `leafUploadBytes` a cull: the view culled again till
   * all are, `frameVersion`), those no cull wanted for `PagedLeafReleaseMs` let go (their value blocks); whether pages
   * were written or let go (their value blocks' table synced again)
   */
  /** the pages whose values are being read (`requestValues`): the view culled again as each arrives */
  const readingPages = new Set<number>();
  /**
   * whether wanted page `p`'s values are read (a deferred list's page: `PagedCullSource.requestValues`, begun once):
   * kept wanted, written when they arrive (the view culled again then)
   */
  function readValues(p: number) {
    if (readingPages.has(p)) return true;
    const read = source?.requestValues?.(p);
    if (!read) return false;
    readingPages.add(p);
    void resources.pending.track(
      read.then(
        () => {
          readingPages.delete(p);
          frameVersion += 1;
        },
        () => readingPages.delete(p),
      ),
      'paged cull: a page read',
    );
    return true;
  }
  function updateResidency() {
    const src = source;
    if (!src) return false;
    const at = performance.now();
    let changed = false;
    const pageBytes = valueBytes();
    const toWrite: number[] = [];
    // a template drawn whole (no coarse level: not a leaf block) of a deferred list: its pages not read wanted
    if (src.requestValues) {
      for (const t of src.templates ?? []) {
        if (isLeafBlock(t) || t.copies.poses.length === 0) continue;
        for (let i = 0; i < t.pages.length; ++i) {
          const p = t.pages[i];
          if (unwrittenPages.has(p) && !src.values[p]) pendingWants.add(p);
        }
      }
    }
    for (const p of pendingWants) {
      if (unwrittenPages.has(p) && src.values[p]) toWrite.push(p);
      else if (unwrittenPages.has(p) && readValues(p)) {
        // its values kept (`releaseValues`) set at once
        if (src.values[p]) toWrite.push(p);
      } else pendingWants.delete(p);
    }
    if (toWrite.length > 0) {
      const upload = Math.min(toWrite.length, Math.max(1, Math.floor(leafUploadBytes / pageBytes)));
      // within the budget: the pages not wanted lately let go for room, a wanted one past it left unwritten (its copies
      // coarse; not culled again for it: tried again as the wants change)
      let room = blocksBudget() - (blockCount - freeBlocks.length);
      if (room < upload) {
        const released = releaseForRoom(src, upload - room);
        room += released;
        if (released > 0) changed = true;
      }
      const count = Math.max(0, Math.min(upload, room));
      // every wanted page's block within the budget made room for at once (one growth, not one a cull)
      reserveBlocks(Math.min(toWrite.length, Math.max(0, room)));
      for (let i = 0; i < count; ++i) {
        const p = toWrite[i];
        writeValues(p, src.values[p]!);
        unwrittenPages.delete(p);
        pendingWants.delete(p);
        allocateSlots(p, poseCount(src, p));
        pendingPages.add(p);
      }
      if (count > 0) changed = true;
      // the occlusion bits grown ahead by the entries they bring (a page a copy lacking pages): every selected entry
      // tested against the pyramid in its first cull. Why: past the bits' room drawn untested, the yeast cell's first
      // near cull wanted 38.9 M frame slots (1.5 M after), its rooms and the raster's grown to them
      bitsWanted += Math.ceil(count * Math.max(1, lackingCount) * 1.5);
      if (count === upload && count < toWrite.length) {
        // pages left to write: the view culled again (a still view's too) and a draw asked for
        frameVersion += 1;
        void resources.pending.track(Promise.resolve(), 'paged cull: leaf pages written');
      }
    }
    // the leaf pages no cull wanted for a while let go
    for (const t of src.templates ?? []) {
      if (!isLeafBlock(t)) continue;
      for (let i = 0; i < t.pages.length; ++i) {
        const p = t.pages[i];
        if (unwrittenPages.has(p) || at - (wantedAt.get(p) ?? Number.NEGATIVE_INFINITY) <= PagedLeafReleaseMs) continue;
        if (wantsReads - (wantedRead.get(p) ?? 0) < 2) continue;
        releaseLeaf(p);
        changed = true;
      }
    }
    return changed;
  }
  const NoChanges: PagedCullChanges = { written: [], released: [], runs: [], extras: [], grown: false };

  return {
    sync(src, changes, reaches) {
      source = src;
      const format = src.format;
      if (!format) return;
      positionWords = format.positions === 'grid' ? 2 : 3;
      rotationWords = format.rotations === 'snorm16' ? 2 : format.rotations === 'float32' ? 4 : 0;
      ValueCell.updateIfChanged(uniforms.uPositions, format.positions === 'grid' ? 0 : 1);
      ValueCell.updateIfChanged(uniforms.uRotations, rotationWords === 4 ? 1 : rotationWords === 2 ? 2 : 0);
      // nothing the cull reads changed (a marker's change: the object updated): none of it again. Why: the active
      // entries and the copies' poses collected again were most of a hover's frame with many copies (the isg cell on
      // a 16x16x16 grid, 12 k copies: 7.2 ms)
      if (!forceSync && everSynced && src === syncedSource && isUnchanged(changes) && sameReaches(src, reaches)) return;
      forceSync = false;
      syncedSource = src;
      syncedReaches = Array.from({ length: src.pool.count }, (_, m) => reaches[m] ?? 0);

      // the tables: written whole where a section's capacity changed, else by what changed
      // a pose's segment its member in 16 bits (`PagedSegmentMask`: the copy of a pick above it)
      if (src.pool.count > PagedSegmentMask + 1) {
        throw new Error(`paged cull: ${src.pool.count} members, at most ${PagedSegmentMask + 1}`);
      }
      const members = Math.max(src.pool.count, 1);
      const memberCapacity = Math.max(64, 2 ** Math.ceil(Math.log2(members)));
      const { copies, instances } = collectActive(src);
      const activeCapacity = Math.max(64, 2 ** Math.ceil(Math.log2(Math.max(1, activeCount))));
      const copyCapacity = copies > 0 ? Math.max(16, 2 ** Math.ceil(Math.log2(copies))) : 0;
      const whole =
        src.pageCapacity !== tablesPages ||
        src.runCapacity !== tablesRuns ||
        memberCapacity !== tablesMembers ||
        activeCapacity !== tablesActive ||
        copyCapacity !== tablesCopies;
      // the frame's slots: what a cull needs (`sizeFrame`)
      frameInstances = instances;
      sizeFrame();
      if (whole) {
        tablesPages = src.pageCapacity;
        tablesRuns = src.runCapacity;
        tablesMembers = memberCapacity;
        tablesExtras = src.extraCapacity;
        tablesActive = activeCapacity;
        tablesCopies = copyCapacity;
        const runsBase = tablesPages * RecordVec4s;
        const membersBase = runsBase + tablesRuns;
        const activeBase = membersBase + tablesMembers;
        ValueCell.updateIfChanged(uniforms.uRunsBase, runsBase);
        ValueCell.updateIfChanged(uniforms.uMembersBase, membersBase);
        ValueCell.updateIfChanged(uniforms.uActiveBase, activeBase);
        const poseFirstsBase = activeBase + Math.ceil((tablesActive * 2) / 4);
        ValueCell.updateIfChanged(uniforms.uPoseFirstsBase, poseFirstsBase);
        const valueBlocksBase = poseFirstsBase + Math.ceil(tablesPages / 4);
        ValueCell.updateIfChanged(uniforms.uValueBlocksBase, valueBlocksBase);
        const clipBase = valueBlocksBase + Math.ceil(tablesPages / 4);
        ValueCell.updateIfChanged(uniforms.uClipBase, clipBase);
        const rasterRowsBase = clipBase + Clip.MaxObjectCount * ClipStorageVec4s;
        ValueCell.updateIfChanged(uniforms.uRasterRowsBase, rasterRowsBase);
        const selectedBase = rasterRowsBase + tablesMembers;
        ValueCell.updateIfChanged(uniforms.uSelectedBase, selectedBase);
        const copiesBase = selectedBase + 1 + Math.ceil(tablesActive / 4);
        ValueCell.updateIfChanged(uniforms.uCopiesBase, copiesBase);
        // the extras last: grown in place, no other section moved (`growExtras`)
        const extrasBase = copiesBase + tablesCopies * CopyVec4s;
        ValueCell.updateIfChanged(uniforms.uExtrasBase, extrasBase);
        const vec4s = extrasBase + Math.ceil((tablesExtras * 2) / 4);
        ensure(tables, vec4s * 16);
        writeClip();
        writeRasterRows();
        tables.ref.value.write(src.records.subarray(0, src.pageCapacity * RecordVec4s * 4), 0);
        tables.ref.value.write(src.runs.subarray(0, src.runCapacity * 4), runsBase * 16);
        if (tablesExtras > 0) tables.ref.value.write(src.extras.subarray(0, tablesExtras * 2), extrasBase * 16);
        membersWhole = true;
      } else {
        if (src.extraCapacity !== tablesExtras) growExtras(src.extraCapacity);
        const recordWords = RecordVec4s * 4;
        for (const p of [...changes.written, ...changes.released]) {
          tables.ref.value.write(src.records.subarray(p * recordWords, (p + 1) * recordWords), p * recordWords * 4);
        }
        const r = changes.runs;
        const runsBase = uniforms.uRunsBase.ref.value;
        for (let k = 0; k < r.length; k += 2) {
          tables.ref.value.write(src.runs.subarray(r[k] * 4, r[k + 1] * 4), (runsBase + r[k]) * 16);
        }
        const e = changes.extras;
        const extrasBase = uniforms.uExtrasBase.ref.value * 16;
        for (let k = 0; k < e.length; k += 2) {
          tables.ref.value.write(src.extras.subarray(e[k] * 2, e[k + 1] * 2), extrasBase + e[k] * 8);
        }
      }
      writeMembers(src, reaches, whole ? undefined : changes.members);
      writeActive();
      // a page's padding: the largest reach times the largest scale of the runs and extras so far (only grown: a scale
      // gone kept, conservative), by what this sync uploaded (all where written whole)
      const runScales = new Float32Array(src.runs.buffer, src.runs.byteOffset, src.runs.length);
      const extraScales = new Float32Array(src.extras.buffer, src.extras.byteOffset, src.extras.length);
      const scanScales = (scales: Float32Array, words: number, from: number, to: number) => {
        const end = Math.min(to * words, scales.length);
        for (let k = from * words + words - 1; k < end; k += words) {
          const v = scales[k];
          if (v > maxScale) maxScale = v;
        }
      };
      if (whole || !scannedScales) {
        scanScales(runScales, 4, 0, src.runCapacity);
        scanScales(extraScales, 2, 0, src.extraCapacity);
        scannedScales = true;
      } else {
        for (let k = 0; k < changes.runs.length; k += 2) scanScales(runScales, 4, changes.runs[k], changes.runs[k + 1]);
        for (let k = 0; k < changes.extras.length; k += 2) {
          scanScales(extraScales, 2, changes.extras[k], changes.extras[k + 1]);
        }
      }
      let maxReach = 0;
      for (let m = 0; m < src.pool.count; ++m) maxReach = Math.max(maxReach, reaches[m] ?? 0);
      ValueCell.updateIfChanged(uniforms.uMaxReach, maxReach);
      memberCount = src.pool.count;

      // the values by page (their offsets kept as the store grows): the buffers grown by a copy of what was there (the
      // pages before this sync's), submitted before what is written now
      // (a write goes ahead of the commands not submitted: the copy would overwrite it)
      const kept = keptPages;
      const groupBytes = PagedCullGroupsPerPage * 32;
      const copied = ensure(groupBoxes, src.pageCapacity * groupBytes, kept * groupBytes);
      if (copied && kept > 0) compute.submit();
      keptPages = src.pageCount;
      // the pages released their slots let go, the written given slots (a culler's first sync every page: the changes
      // taken before, another culler's of the store, one made again)
      const first = !everSynced;
      everSynced = true;
      for (const p of changes.released) {
        freeSlots(p);
        freeBlock(p);
        pendingPages.delete(p);
        unwrittenPages.delete(p);
      }
      const written = first ? Array.from({ length: src.pageCount }, (_, p) => p) : changes.written;
      // the pages of leaf blocks: written where a cull wants them (`updateResidency`), not here
      const deferred = new Set<number>();
      for (const t of src.templates ?? []) {
        if (!isLeafBlock(t)) continue;
        for (let i = 0; i < t.pages.length; ++i) deferred.add(t.pages[i]);
      }
      // the value blocks of the pages written made room for at once (one growth, one copy)
      let newBlocks = 0;
      for (const p of written) {
        if (src.values[p] && !deferred.has(p) && (p >= pageBlocks.length || pageBlocks[p] === NoSlots)) newBlocks += 1;
      }
      reserveBlocks(newBlocks);
      for (const p of written) {
        const values = src.values[p];
        // a leaf page (its values read or not: a deferred list's read when wanted)
        if (deferred.has(p)) {
          freeSlots(p);
          freeBlock(p);
          pendingPages.delete(p);
          unwrittenPages.add(p);
          continue;
        }
        if (!values) {
          // a deferred list's page not read (a template's not in a leaf block: its copies not set yet, or no coarse
          // level): wanted where its template is drawn whole (`updateResidency`)
          if (src.requestValues) unwrittenPages.add(p);
          continue;
        }
        unwrittenPages.delete(p);
        writeValues(p, values);
        allocateSlots(p, poseCount(src, p));
        pendingPages.add(p);
      }
      writePoseFirsts();
      instanceTotal = 0;
      if (memberCounts.length < src.pool.count) memberCounts = new Int32Array(src.pool.counts.length);
      for (let m = 0; m < src.pool.count; ++m) {
        memberCounts[m] = src.pool.counts[m];
        instanceTotal += src.pool.counts[m];
      }
      syncVersion += 1;
    },

    cull(view, pyramid) {
      // the leaf blocks this view draws written (within a budget a cull: their copies drawn coarse till all are)
      // a cull's counts growing the rooms where its pyramid is of the same pages: not the first after a leaf block was
      // made resident (its pyramid of the copies drawn coarse, a sample's gaps hiding little: the yeast's first near
      // cull 27.6 M frame slots, 1.5 M after, its rooms and the raster's chunks grown to it), culled again instead
      residencyChanged = updateResidency();
      if (residencyChanged && source) {
        forceSync = true;
        this.sync(source, NoChanges, syncedReaches);
      }
      culled = activeCount > 0 && memberCount > 0;
      culledWithPyramid = culled && !!pyramid && view.secondPhase !== false;
      pyramidToCome = !pyramid && !!view.pyramidToCome;
      if (!culled) return;
      // the frame's slots grown where a cull read back needed more (`readCursors`)
      sizeFrame();
      compose();
      cullKernel ??= compute.createKernel(PagedCullKernel, PagedCullSchema, cullValues);
      const planes = [
        uniforms.uPlane0,
        uniforms.uPlane1,
        uniforms.uPlane2,
        uniforms.uPlane3,
        uniforms.uPlane4,
        uniforms.uPlane5,
      ];
      for (let i = 0; i < 6; ++i) {
        const { normal, constant } = view.frustum[i as 0];
        ValueCell.update(planes[i], Vec4.set(planes[i].ref.value, normal[0], normal[1], normal[2], constant));
      }
      const { normal, constant } = view.cameraPlane;
      ValueCell.update(
        uniforms.uCameraPlane,
        Vec4.set(uniforms.uCameraPlane.ref.value, normal[0], normal[1], normal[2], constant),
      );
      ValueCell.update(uniforms.uLevels, Mat4.copy(uniforms.uLevels.ref.value, view.levels));
      const objectTransform = view.objectTransform ?? IdentityTransform;
      ValueCell.update(uniforms.uObjectTransform, Mat4.copy(uniforms.uObjectTransform.ref.value, objectTransform));
      ValueCell.updateIfChanged(uniforms.uObjectScale, Mat4.getMaxScaleOnAxis(objectTransform));
      setOrderRange(view, objectTransform);

      ValueCell.updateIfChanged(uniforms.uLevelCount, view.levelCount);
      const levels = Math.max(1, view.levelCount);
      const ghost = !!view.ghost;
      const raster = !!view.rasterRegions;
      // the regions: the levels', the raster's (level-major), the ghost's
      const rasterRegions = raster ? levels * memberCount : 0;
      regionCount = memberCount * levels + rasterRegions + (ghost ? memberCount : 0);
      ValueCell.updateIfChanged(uniforms.uGhost, ghost ? 1 : 0);
      ValueCell.updateIfChanged(uniforms.uGhostBase, memberCount * levels + rasterRegions);
      ValueCell.updateIfChanged(uniforms.uMemberCount, memberCount);
      ValueCell.updateIfChanged(uniforms.uRasterBase, memberCount * levels);
      const split = raster ? view.raster : undefined;
      ValueCell.updateIfChanged(uniforms.uRasterRadius, split?.radius ?? 0);
      ValueCell.updateIfChanged(uniforms.uRasterScale, split?.scale ?? 0);
      ValueCell.updateIfChanged(uniforms.uRasterNear, split?.near ?? 0);
      ValueCell.updateIfChanged(uniforms.uRasterOrtho, split?.ortho ? 1 : 0);
      ValueCell.updateIfChanged(uniforms.uCoarseScale, view.coarse?.scale ?? 0);
      ValueCell.updateIfChanged(uniforms.uCoarseOrtho, view.coarse?.ortho ? 1 : 0);
      ValueCell.updateIfChanged(uniforms.uCoarsePixels, view.coarse?.pixels ?? 0);
      // the phases' cursors (the second's bound as the first's prior)
      ValueCell.update(counts, firstCounts);
      ValueCell.update(prior, secondCounts);
      // the frame's cursors after the regions' (the phases', the fallback's)
      ensure(counts, regionCount * 4 + 8);
      ensure(prior, regionCount * 4 + 8);
      ValueCell.updateIfChanged(uniforms.uFrameCursor, regionCount);
      firstCounts = counts.ref.value;
      secondCounts = prior.ref.value;
      ValueCell.updateIfChanged(uniforms.uRegionCount, regionCount);
      layoutBases(levels, raster, ghost);
      // the indices: the regions' rooms (`layoutBases`); within the device's limit (the slots past it not written,
      // their regions' counts cut)
      const maxSlots = Math.floor(compute.maxStorageBufferSize / 4);
      indexRoom = Math.min(Math.max(1, indexTotal), maxSlots);
      ValueCell.updateIfChanged(uniforms.uIndexRoom, indexRoom);
      ensure(indices, indexRoom * 4);

      // the select's: a thread an active entry; the cull's rows of `PagedCullMaxRow` workgroups (its dispatch's x where
      // in 2d: `paged-cull-select-args`)
      const { x, y } = getDispatchSize(Math.ceil(Math.max(1, activeCount) / 64), PagedCullMaxRow);
      ValueCell.updateIfChanged(uniforms.uActiveCount, activeCount);
      ValueCell.updateIfChanged(selectValues.uRowThreads, x * 64);
      ValueCell.updateIfChanged(uniforms.uRowThreads, PagedCullMaxRow * PagedCullWorkgroupKeys);

      dispatchSize = [x, y];
      // the first phase: the last frame's pyramid, the occluded marked for the second
      setHiZValues(uniforms, pyramid, emptyHiZ);
      ValueCell.updateIfChanged(uniforms.uPhase, 0);
      // a bit a selected entry's particle, then a bit its group's (occluded whole in the first phase): by selected
      // entry, not store index (a template's particles drawn once a copy); the entries `bitsRoom` holds, within the
      // active ones and the device's limit what the words after them leave (the buffer bound whole)
      // the dropped copy groups after them (a word each: the fallback's), cleared every cull
      const droppedWords = Math.max(1, groupCount);
      // the copies occluded whole by the pyramid after them (a word each: the second phase's select)
      const copyCount = Math.max(1, copyPoses.length / (CopyVec4s * 4));
      // two words for the sections' rounding up
      const bitsWords = Math.floor(compute.maxStorageBufferSize / 4) - droppedWords - copyCount - 2;
      const bitsLimit = Math.floor((bitsWords * 32) / (PagedCullPageSize + PagedCullGroupsPerPage));
      bitsRoom = Math.max(1, Math.min(bitsWanted, activeCount, bitsLimit));
      const particleWords = Math.ceil(bitsRoom * (PagedCullPageSize / 32));
      const droppedBase = particleWords + Math.ceil((bitsRoom * PagedCullGroupsPerPage) / 32);
      ensure(occludedBits, (droppedBase + droppedWords + copyCount) * 4);
      ValueCell.updateIfChanged(selectValues.uCopyMarksBase, droppedBase + droppedWords);
      ValueCell.updateIfChanged(uniforms.uPageBitsBase, particleWords);
      ValueCell.updateIfChanged(uniforms.uBitsRoom, bitsRoom);
      // the survivors after the regions' firsts, counts and rooms (their count, then a word a group): the groups of the
      // entries a cull selected (read back, grown: `readSelected`), not of every active entry (the isg 16x16x16's
      // 2.7 M: 340 MB); past it a frame's groups not culled, culled again
      const survivorsBase = Math.max(1, regionCount) * 3;
      survivorsEntries = Math.max(1, Math.min(survivorsWanted, activeCount));
      const survivorsRoom = survivorsEntries * PagedCullGroupsPerPage;
      ensure(bases, (survivorsBase + 1 + PagedCullOrderBuckets + survivorsRoom) * 4, survivorsBase * 4);
      ValueCell.updateIfChanged(uniforms.uSurvivorsBase, survivorsBase);
      ValueCell.updateIfChanged(uniforms.uSurvivorsRoom, survivorsRoom);
      ValueCell.updateIfChanged(uniforms.uDroppedBase, droppedBase);
      if (pyramid) occludedBits.ref.value.clear();
      else occludedBits.ref.value.clear(droppedBase * 4, droppedWords * 4);
      wants();
      select();
      readWants();
      runPass('paged cull', true);
      // no second phase to follow: this phase's slots the frame's (read before the fallback's), its selected entries
      // the survivors' room's (`PagedCullGroupsKernel`) and with a pyramid the bits'; else read where the phases end
      // (one read at a time: the second's not made while the first's is on its way)
      if (!culledWithPyramid) {
        readSelected(!!pyramid);
        readCursors(firstCounts, !pyramidToCome && !residencyChanged, false, !!pyramid);
        fallback();
      }
    },

    endFirst() {
      if (!culledWithPyramid) return;
      culledWithPyramid = false;
      readSelected(true);
      readCursors(firstCounts, !residencyChanged, false, true);
    },

    cullSecond(pyramid) {
      if (!culledWithPyramid || !cullKernel) return false;
      culledWithPyramid = false;
      // the second's cursors in `counts`, the first's counts after the regions' firsts (its slots past them)
      firstCounts.copyTo(bases.ref.value, 0, regionCount * 4, regionCount * 4);
      // the frame's cursors the first's: the second's frame slots past them (the fallback's none yet)
      firstCounts.copyTo(secondCounts, regionCount * 4, regionCount * 4, 8);
      ValueCell.update(counts, secondCounts);
      setHiZValues(uniforms, pyramid, emptyHiZ);
      ValueCell.updateIfChanged(uniforms.uPhase, 1);
      // the copies the first left out that the frame's pyramid shows, after the first's entries
      select(true);
      runPass('paged cull, second phase', false);
      // the first's entries and those the second added: the survivors' room for both, the bits' for the first's (the
      // count of both: grown for a few more)
      readSelected(true);
      readCursors(secondCounts, !residencyChanged, true, true);
      fallback();
      ValueCell.updateIfChanged(uniforms.uPhase, 0);
      ValueCell.update(counts, firstCounts);
      return true;
    },
    setRasterRows(rows) {
      rasterRows = rows;
      writeRasterRows();
    },
    setPosesFront(front) {
      const from = posesFront.byteLength;
      posesFront = front;
      ValueCell.updateIfChanged(uniforms.uPoseBase, front.length / 4);
      // another count of vec4s before the poses: the poses moved after it (a copy), within the device's limit
      if (front.byteLength !== from || poses.ref.value.byteLength < front.byteLength + slots.capacity * poseBytes) {
        const capacity = Math.min(slots.capacity, maxPoseSlots());
        sizePoses(capacity, from, front.byteLength);
        ValueCell.updateIfChanged(uniforms.uPoseRoom, capacity);
      } else {
        writeFront();
      }
    },
    setClip(c) {
      clip = c;
      ValueCell.updateIfChanged(uniforms.uClipObjectCount, Math.min(c.count, Clip.MaxObjectCount));
      ValueCell.updateIfChanged(uniforms.uClipInstance, c.instance && c.count > 0 ? 1 : 0);
      writeClip();
    },
    get counts() {
      return firstCounts;
    },
    readBack,
    get secondCounts() {
      return secondCounts;
    },
    get bases() {
      return bases.ref.value;
    },
    get indices() {
      return indices.ref.value;
    },
    get maxScale() {
      return maxScale;
    },
    get indexRoom() {
      return indexRoom;
    },
    get poses() {
      return poses.ref.value;
    },
    get instanceTotal() {
      return instanceTotal;
    },
    get rasterRoom() {
      return rasterRoom;
    },
    get regionRooms() {
      return regionRooms;
    },
    get regionCount() {
      return regionCount;
    },
    get memberCount() {
      return memberCount;
    },
    get frameVersion() {
      return frameVersion;
    },
    async readSelectedEntries() {
      const base = uniforms.uSelectedBase.ref.value * 16;
      compute.submit();
      const count = new Uint32Array(await tables.ref.value.readAsync(base, 16))[0];
      if (count === 0) return [];
      const list = new Uint32Array(await tables.ref.value.readAsync(base + 16, count * 4));
      return Array.from(list, (a) => ({
        page: activeEntries[a * 2] & PagedEntryPageMask,
        copy: activeEntries[a * 2 + 1],
      }));
    },
    get valueBlocks() {
      return blockCount - freeBlocks.length;
    },
    get deferredPages() {
      return unwrittenPages.size;
    },
    get culled() {
      return culled;
    },
    destroy() {
      composeKernel?.destroy();
      selectKernel?.destroy();
      selectArgsKernel?.destroy();
      cullKernel?.destroy();
      wantsKernel?.destroy();
      groupsKernel?.destroy();
      for (const c of [
        tables,
        positions,
        composePages,
        bases,
        poses,
        groupBoxes,
        indices,
        occludedBits,
        wanted,
        lacking,
      ]) {
        c.ref.value.destroy();
      }
      firstCounts.destroy();
      secondCounts.destroy();
      selectArgs.ref.value.destroy();
      emptyHiZ.destroy();
    },
  };
}
