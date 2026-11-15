/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
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

// The culling of a page store's particles (`ParticlePageStore`): a thread a particle, nothing of
// the instances laid out on the CPU. A particle's pose composed on the GPU when its page is written (`compose`, at
// the next cull), at a slot of its page's range (a range a page, its count's: a pose a particle),
// not each frame; each frame a pass a phase: a particle's run and position as stored tested (8 B a particle, not its
// pose's 32: the yeast's pass 2.6-3.0 ms against 3.0-3.3), a
// visible instance's pose slot written into its region (a member and level, the raster's, the ghost's) at its cursor;
// the regions laid out on the CPU as the merged culler's (a block a level, a member's room in it its instances: a
// member's instances all at one level fit; a level's raster region the same room from its end), so no count and no
// prefix sum before. The draws read a pose by the slot (`CulledData`, `CulledPagedFlag`). The store's tables in one
// buffer (the records, the runs, the members, the active pages, the pages' pose slots: the storage buffers a stage has
// are few), its positions and rotations in one each (by store index), its extras in one.

/** a store page's particles as bits (`StorePageBits`): 8k, a tile its own pages (a 16³ grid's yeast tile 11-20k) */
export const PagedCullPageBits = 13;
/** particles a store page */
export const PagedCullPageSize = 1 << PagedCullPageBits;
/** a workgroup's invocations and its table of regions (the instances of a region a workgroup counted at once) */
export const PagedCullWorkgroupKeys = 256;
/** words of a page's record (`PageRecordWords`), and its runs' count's flag of its members a particle */
export const PagedRecordWords = 16;
export const PagedPerParticleMembers = 0x80000000;
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

/** a workgroup a group of a page's particles (the cull's): its centers' box reduced (`groupBoxes`) */
export const PagedCullComposeKernel = KernelCode('paged-cull-compose', [PagedCullWorkgroupKeys]);
export const PagedCullKernel = KernelCode('paged-cull', [PagedCullWorkgroupKeys]);
export const PagedCullArgsKernel = KernelCode('paged-cull-args', [64]);

const PagedCullUniforms = {
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
  /** the active pages: the dispatch's pages, a page's particles a thread each */
  uActiveCount: UniformSpec('i'),
  /** threads a row of workgroups (the dispatch in 2d past the per-dimension limit) */
  uRowThreads: UniformSpec('i'),
  /** the tables' sections, in vec4s: the runs, the members (center, reach), the active pages and the pages' first pose slots (four a vec4) */
  uRunsBase: UniformSpec('i'),
  uMembersBase: UniformSpec('i'),
  uActiveBase: UniformSpec('i'),
  uPoseFirstsBase: UniformSpec('i'),
  /** positions: 0 a grid's words (two a particle), 1 float32 (three) */
  uPositions: UniformSpec('i'),
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
  uRasterBase: UniformSpec('i'),
  uRasterRowsBase: UniformSpec('i'),
  /**
   * a page's sphere padded by this: the largest member's reach times the largest scale (a particle's sphere's
   * largest radius), what a page tested whole against the pyramid takes (none hidden whose particles show)
   */
  uMaxReach: UniformSpec('f'),
  /** the groups' bits in `occludedBits` from this word: a group occluded whole in the first phase */
  uPageBitsBase: UniformSpec('i'),
  /** the regions (`regions`: their firsts, then the first phase's counts) */
  uRegionCount: UniformSpec('i'),
  ...HiZSchema,
};

/** The poses of the pages in `composePages` (a page and its first pose slot each), from their values as stored */
export const PagedCullComposeSchema = {
  uRowThreads: UniformSpec('i'),
  uComposeCount: UniformSpec('i'),
  uRunsBase: UniformSpec('i'),
  uMembersBase: UniformSpec('i'),
  uPoseBase: UniformSpec('i'),
  uPoseRoom: UniformSpec('i'),
  /** positions: 0 a grid's words (two a particle), 1 float32 (three) */
  uPositions: UniformSpec('i'),
  /** rotations: 0 none, 1 float32 (four words a particle), 2 snorm16 (two) */
  uRotations: UniformSpec('i'),
  tables: StorageBufferSpec('vec4<u32>', 'read'),
  positions: StorageBufferSpec('u32', 'read'),
  rotations: StorageBufferSpec('u32', 'read'),
  extrasWords: StorageBufferSpec('u32', 'read'),
  composePages: StorageBufferSpec('u32', 'read'),
  poses: StorageBufferSpec('vec4<u32>', 'read-write'),
  /** per group (`PagedCullGroupsPerPage` a page) its instances' centers' box, two vec4s (min, max) */
  groupBoxes: StorageBufferSpec('vec4<f32>', 'read-write'),
};

export const PagedCullSchema = {
  ...PagedCullUniforms,
  tables: StorageBufferSpec('vec4<u32>', 'read'),
  positions: StorageBufferSpec('u32', 'read'),
  extrasWords: StorageBufferSpec('u32', 'read'),
  /** the regions' cursors (the phase's) */
  counts: StorageBufferSpec('atomic<u32>', 'read-write'),
  /**
   * the regions' firsts (a raster region's its end, laid out on the CPU), then the first phase's counts (copied there
   * for the second, its slots past them): one binding, the stage's eight
   */
  regions: StorageBufferSpec('u32', 'read'),
  /** the groups' boxes (the compose's) */
  groupBoxes: StorageBufferSpec('vec4<f32>', 'read'),
  /** a bit a store index: occluded in the first phase (`uPhase` 0), what the second tests */
  occludedBits: StorageBufferSpec('atomic<u32>', 'read-write'),
  /** the visible instances' pose slots, by region */
  indices: StorageBufferSpec('u32', 'read-write'),
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
  counts: StorageBufferSpec('u32', 'read'),
  bases: StorageBufferSpec('u32', 'read'),
  prior: StorageBufferSpec('u32', 'read'),
  args: StorageBufferSpec('u32', 'read-write'),
};

export const PagedCullKernels: readonly KernelRequest[] = [
  { code: PagedCullComposeKernel, schema: PagedCullComposeSchema },
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
  /** a page's values as stored (none: free) */
  readonly values: readonly (
    | { readonly positions: ArrayBufferView; readonly rotations?: ArrayBufferView }
    | undefined
  )[];
  readonly pageCount: number;
  /** the members: their centers (three a member) and their count */
  readonly pool: { readonly centers: Float64Array; readonly count: number; readonly counts: Int32Array };
}

/** What changed in the source since the last sync (`PageStoreChanges`) */
export interface PagedCullChanges {
  readonly written: readonly number[];
  readonly released: readonly number[];
  readonly runs: readonly number[];
  readonly extras: readonly number[];
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
  readonly raster?: { readonly radius: number; readonly scale: number; readonly near: number };
}

/** The clip objects a cull is of (`packStorageClipObjects`'s, `ClipStorageVec4s` a clip object) */
export interface PagedCullClip {
  readonly data: Float32Array;
  readonly count: number;
  /** instance clipping: an instance's center tested too */
  readonly instance: boolean;
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
  readonly counts: StorageBuffer;
  readonly secondCounts: StorageBuffer;
  /** the regions' firsts in `indices` (a raster region's its end) */
  readonly bases: StorageBuffer;
  /** the visible instances' pose slots by region (a block of regions a level and its raster's, the ghost's) */
  readonly indices: StorageBuffer;
  /** the slots of `indices` written (within the device's buffer limit), as of the last cull */
  readonly indexRoom: number;
  /** the composed poses, `PagedPoseWords` u32s a slot, after the front (`setPosesFront`) */
  readonly poses: StorageBuffer;
  /** the store's instances: a block of regions' slots (a level's), as of the last sync */
  readonly instanceTotal: number;
  /** regions: members times levels (at least one), the raster's, and a member the ghost's after, as of the last cull */
  readonly regionCount: number;
  /** the members of the regions, as of the last cull */
  readonly memberCount: number;
  /** whether the last `cull` culled (a store with pages and members) */
  readonly culled: boolean;
  destroy: () => void;
}

export function createPagedCuller(compute: ComputeContext, resources: GraphicsResources): PagedCuller {
  // a storage buffer's size a multiple of 4: of 16 (a vec4's) where that is within the device's binding limit, else of
  // 4 (the limit 2^31 - 4: rounded up to 16, 2^31, refused); a growth's 1.5 times capped at the limit by its caller
  const createBuffer = (bytes: number, label: string) => {
    const aligned = Math.ceil(bytes / 16) * 16;
    const size = aligned <= compute.maxStorageBufferSize ? aligned : Math.ceil(bytes / 4) * 4;
    return compute.createBuffer(Math.max(16, size), label);
  };
  const cell = (label: string) => ValueCell.create(createBuffer(16, label));
  const tables = cell('paged cull tables');
  const positions = cell('paged cull positions');
  const rotations = cell('paged cull rotations');
  const extras = cell('paged cull extras');
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

  const uniforms = {
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
    uPositions: ValueCell.create(0),
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
    uRasterBase: ValueCell.create(0),
    uRasterRowsBase: ValueCell.create(0),
    uMaxReach: ValueCell.create(0),
    uPageBitsBase: ValueCell.create(0),
    uRegionCount: ValueCell.create(0),
    ...createHiZValues(emptyHiZ),
  };
  const composeValues = {
    uRowThreads: ValueCell.create(0),
    uComposeCount: ValueCell.create(0),
    uRunsBase: uniforms.uRunsBase,
    uMembersBase: uniforms.uMembersBase,
    uPoseBase: uniforms.uPoseBase,
    uPoseRoom: uniforms.uPoseRoom,
    uPositions: uniforms.uPositions,
    uRotations: ValueCell.create(0),
    tables,
    positions,
    rotations,
    extrasWords: extras,
    composePages,
    poses,
    groupBoxes,
  };
  const cullValues = {
    ...uniforms,
    tables,
    positions,
    extrasWords: extras,
    counts,
    regions: bases,
    groupBoxes,
    occludedBits,
    indices,
  };
  let composeKernel: ComputeKernel | undefined;
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
  let positionWords = 2;
  let rotationWords = 0;
  let activeCount = 0;
  let memberCount = 0;
  let instanceTotal = 0;
  let regionCount = 1;
  /** the pages there were at the last sync: what a grown buffer of positions keeps */
  let keptPages = 0;
  let everSynced = false;
  /** the extras there were at the last sync: what a grown buffer keeps */
  let syncedExtras = 0;
  let culled = false;
  let source: PagedCullSource | undefined;
  /** the pages written and not composed yet */
  const pendingPages = new Set<number>();
  /** the members' instances as of the last sync, and that sync's: the regions laid out again where another */
  let memberCounts = new Int32Array(0);
  let syncVersion = 0;
  let basesKey = '';

  // the pose slots: a range a page (its count's), its first in the tables' section (`NoSlots`: none)
  const poseBytes = PagedPoseWords * 4;
  const maxPoseSlots = () =>
    Math.max(1, Math.floor((compute.maxStorageBufferSize - posesFront.byteLength) / poseBytes));
  const slots = createRangeAllocator(0, 1);
  let pageFirsts = new Uint32Array(0);
  let pageSlotCounts = new Int32Array(0);
  /** the poses' buffer for `capacity` slots after the front, what was there copied (the front's poses moved to `toFront`) */
  const sizePoses = (capacity: number, fromFront: number, toFront: number) => {
    const bytes = toFront + capacity * poseBytes;
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
  };

  /** the members' rows as written (center, reach), and whether the next are written whole (the tables made again) */
  let writtenMembers = new Float32Array(0);
  let membersWhole = true;
  /** the largest scale of the runs and extras synced (1: an instance's scale where none) */
  let maxScale = 1;
  let scannedScales = false;
  const writeMembers = (src: PagedCullSource, reaches: ArrayLike<number>) => {
    const n = src.pool.count;
    const data = new Float32Array(Math.max(1, n) * 4);
    const { centers } = src.pool;
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

  const writeActive = (src: PagedCullSource) => {
    const active: number[] = [];
    const recordWords = RecordVec4s * 4;
    for (let p = 0; p < src.pageCount; ++p) if (src.records[p * recordWords + 6] > 0) active.push(p);
    activeCount = active.length;
    tables.ref.value.write(Uint32Array.from(active.length ? active : [0]), uniforms.uActiveBase.ref.value * 16);
  };

  /**
   * The pages pending composed (their values as written by store index, kept): a pass over their particles, at a cull
   * alone (the page list a queue write: one compose between submits, a frame's)
   */
  const compose = () => {
    if (pendingPages.size === 0 || !source?.format) return;
    composeKernel ??= compute.createKernel(PagedCullComposeKernel, PagedCullComposeSchema, composeValues);
    const list: number[] = [];
    for (const p of pendingPages) {
      if (source.values[p] && p < pageFirsts.length && pageFirsts[p] !== NoSlots) list.push(p, pageFirsts[p]);
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
   * The regions' firsts: a level's region (block, member) at the block's first (`level * instanceTotal`) plus the
   * member's (its instances before it), member-major (`member * levels + level`); the level's raster region (from
   * `levels * members`, level-major) its end (its slots down from it); the ghost's (a block after the levels') after
   * them; laid out where the store's members' instances or the blocks changed
   */
  const layoutBases = (levels: number, raster: boolean, ghost: boolean) => {
    const key = `${syncVersion},${levels},${raster},${ghost}`;
    if (key === basesKey) return;
    basesKey = key;
    const n = memberCount;
    const total = instanceTotal;
    const first = new Float64Array(n);
    for (let m = 0, sum = 0; m < n; ++m) {
      first[m] = sum;
      sum += memberCounts[m];
    }
    const out = new Uint32Array(Math.max(1, regionCount));
    const rasterBase = n * levels;
    const ghostBase = rasterBase + (raster ? levels * n : 0);
    // a first past 2^32 slots at it: its region none written (the indices' room far below)
    const at = (block: number, m: number) => Math.min(0xffffffff, block * total + first[m]);
    for (let m = 0; m < n; ++m) {
      for (let l = 0; l < levels; ++l) {
        out[m * levels + l] = at(l, m);
        if (raster) out[rasterBase + l * n + m] = Math.min(0xffffffff, at(l, m) + memberCounts[m]);
      }
      if (ghost) out[ghostBase + m] = at(levels, m);
    }
    // room for the first phase's counts after them (the second's)
    ensure(bases, out.byteLength * 2);
    bases.ref.value.write(out, 0);
  };

  /** the phase's pass: its cursors cleared, its particles culled into the regions */
  let dispatchSize = [1, 1];
  let culledWithPyramid = false;
  const runPass = (label: string) => {
    const [x, y] = dispatchSize;
    counts.ref.value.clear(0, regionCount * 4);
    compute.beginComputePass({ label });
    cullKernel!.dispatch(x, y);
    compute.endComputePass();
  };

  return {
    sync(src, changes, reaches) {
      source = src;
      const format = src.format;
      if (!format) return;
      positionWords = format.positions === 'grid' ? 2 : 3;
      rotationWords = format.rotations === 'snorm16' ? 2 : format.rotations === 'float32' ? 4 : 0;
      ValueCell.updateIfChanged(uniforms.uPositions, format.positions === 'grid' ? 0 : 1);
      ValueCell.updateIfChanged(composeValues.uRotations, rotationWords === 4 ? 1 : rotationWords === 2 ? 2 : 0);

      // the tables: written whole where a section's capacity changed, else by what changed
      const members = Math.max(src.pool.count, 1);
      const memberCapacity = Math.max(64, 2 ** Math.ceil(Math.log2(members)));
      const whole =
        src.pageCapacity !== tablesPages || src.runCapacity !== tablesRuns || memberCapacity !== tablesMembers;
      if (whole) {
        tablesPages = src.pageCapacity;
        tablesRuns = src.runCapacity;
        tablesMembers = memberCapacity;
        const runsBase = tablesPages * RecordVec4s;
        const membersBase = runsBase + tablesRuns;
        const activeBase = membersBase + tablesMembers;
        ValueCell.updateIfChanged(uniforms.uRunsBase, runsBase);
        ValueCell.updateIfChanged(uniforms.uMembersBase, membersBase);
        ValueCell.updateIfChanged(uniforms.uActiveBase, activeBase);
        const poseFirstsBase = activeBase + Math.ceil(tablesPages / 4);
        ValueCell.updateIfChanged(uniforms.uPoseFirstsBase, poseFirstsBase);
        const clipBase = poseFirstsBase + Math.ceil(tablesPages / 4);
        ValueCell.updateIfChanged(uniforms.uClipBase, clipBase);
        const rasterRowsBase = clipBase + Clip.MaxObjectCount * ClipStorageVec4s;
        ValueCell.updateIfChanged(uniforms.uRasterRowsBase, rasterRowsBase);
        const vec4s = rasterRowsBase + tablesMembers;
        ensure(tables, vec4s * 16);
        writeClip();
        writeRasterRows();
        tables.ref.value.write(src.records.subarray(0, src.pageCapacity * RecordVec4s * 4), 0);
        tables.ref.value.write(src.runs.subarray(0, src.runCapacity * 4), runsBase * 16);
        membersWhole = true;
      } else {
        const recordWords = RecordVec4s * 4;
        for (const p of [...changes.written, ...changes.released]) {
          tables.ref.value.write(src.records.subarray(p * recordWords, (p + 1) * recordWords), p * recordWords * 4);
        }
        const r = changes.runs;
        const runsBase = uniforms.uRunsBase.ref.value;
        for (let k = 0; k < r.length; k += 2) {
          tables.ref.value.write(src.runs.subarray(r[k] * 4, r[k + 1] * 4), (runsBase + r[k]) * 16);
        }
      }
      writeMembers(src, reaches);
      writeActive(src);
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
      ValueCell.updateIfChanged(uniforms.uMaxReach, maxReach * maxScale);
      memberCount = src.pool.count;

      // the positions by store index and the extras (their offsets kept as the store grows, its pool copied): the
      // buffers grown by a copy of what was there (the pages before this sync's), submitted before what is written now
      // (a write goes ahead of the commands not submitted: the copy would overwrite it)
      const kept = keptPages;
      let copied = ensure(
        positions,
        src.pageCapacity * PagedCullPageSize * positionWords * 4,
        kept * PagedCullPageSize * positionWords * 4,
      );
      if (rotationWords > 0) {
        copied =
          ensure(
            rotations,
            src.pageCapacity * PagedCullPageSize * rotationWords * 4,
            kept * PagedCullPageSize * rotationWords * 4,
          ) || copied;
      }
      const groupBytes = PagedCullGroupsPerPage * 32;
      copied = ensure(groupBoxes, src.pageCapacity * groupBytes, kept * groupBytes) || copied;
      const keptExtras = everSynced ? syncedExtras : 0;
      copied = ensure(extras, src.extraCapacity * 8, keptExtras * 8) || copied;
      if (copied && (kept > 0 || keptExtras > 0)) compute.submit();
      keptPages = src.pageCount;
      syncedExtras = src.extraCapacity;
      // the pages released their slots let go, the written given slots (a culler's first sync every page: the changes
      // taken before, another culler's of the store, one made again)
      const first = !everSynced;
      everSynced = true;
      for (const p of changes.released) {
        freeSlots(p);
        pendingPages.delete(p);
      }
      const written = first ? Array.from({ length: src.pageCount }, (_, p) => p) : changes.written;
      const recordWords = RecordVec4s * 4;
      for (const p of written) {
        const values = src.values[p];
        if (!values) continue;
        positions.ref.value.write(values.positions, p * PagedCullPageSize * positionWords * 4);
        if (rotationWords > 0 && values.rotations) {
          rotations.ref.value.write(values.rotations, p * PagedCullPageSize * rotationWords * 4);
        }
        allocateSlots(p, src.records[p * recordWords + 6]);
        pendingPages.add(p);
      }
      if (first && src.extraCapacity > 0) {
        extras.ref.value.write(src.extras.subarray(0, src.extraCapacity * 2), 0);
      } else {
        const e = changes.extras;
        for (let k = 0; k < e.length; k += 2) {
          extras.ref.value.write(src.extras.subarray(e[k] * 2, e[k + 1] * 2), e[k] * 8);
        }
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
      culled = activeCount > 0 && memberCount > 0;
      culledWithPyramid = culled && !!pyramid;
      if (!culled) return;
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
      // the phases' cursors (the second's bound as the first's prior)
      ValueCell.update(counts, firstCounts);
      ValueCell.update(prior, secondCounts);
      ensure(counts, regionCount * 4);
      ensure(prior, regionCount * 4);
      firstCounts = counts.ref.value;
      secondCounts = prior.ref.value;
      ValueCell.updateIfChanged(uniforms.uRegionCount, regionCount);
      layoutBases(levels, raster, ghost);
      // the indices: a block of the store's instances a level (its raster's in it too), the ghost's; within the
      // device's limit (the slots past it not written, their regions' counts cut)
      const blocks = levels + (ghost ? 1 : 0);
      const maxSlots = Math.floor(compute.maxStorageBufferSize / 4);
      indexRoom = Math.min(Math.max(1, blocks * instanceTotal), maxSlots);
      ValueCell.updateIfChanged(uniforms.uIndexRoom, indexRoom);
      ensure(indices, indexRoom * 4);

      const groups = activeCount * (PagedCullPageSize / PagedCullWorkgroupKeys);
      const { x, y } = getDispatchSize(groups, 65535);
      ValueCell.updateIfChanged(uniforms.uActiveCount, activeCount);
      ValueCell.updateIfChanged(uniforms.uRowThreads, x * PagedCullWorkgroupKeys);

      dispatchSize = [x, y];
      // the first phase: the last frame's pyramid, the occluded marked for the second
      setHiZValues(uniforms, pyramid, emptyHiZ);
      ValueCell.updateIfChanged(uniforms.uPhase, 0);
      // a bit a particle, then a bit a group (occluded whole in the first phase)
      const pageCapacity = source?.pageCapacity ?? 1;
      const particleWords = Math.ceil(pageCapacity * (PagedCullPageSize / 32));
      ensure(occludedBits, (particleWords + Math.ceil((pageCapacity * PagedCullGroupsPerPage) / 32)) * 4);
      ValueCell.updateIfChanged(uniforms.uPageBitsBase, particleWords);
      if (pyramid) occludedBits.ref.value.clear();
      runPass('paged cull');
    },

    cullSecond(pyramid) {
      if (!culledWithPyramid || !cullKernel) return false;
      culledWithPyramid = false;
      // the second's cursors in `counts`, the first's counts after the regions' firsts (its slots past them)
      firstCounts.copyTo(bases.ref.value, 0, regionCount * 4, regionCount * 4);
      ValueCell.update(counts, secondCounts);
      setHiZValues(uniforms, pyramid, emptyHiZ);
      ValueCell.updateIfChanged(uniforms.uPhase, 1);
      runPass('paged cull, second phase');
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
    get secondCounts() {
      return secondCounts;
    },
    get bases() {
      return bases.ref.value;
    },
    get indices() {
      return indices.ref.value;
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
    get regionCount() {
      return regionCount;
    },
    get memberCount() {
      return memberCount;
    },
    get culled() {
      return culled;
    },
    destroy() {
      composeKernel?.destroy();
      cullKernel?.destroy();
      for (const c of [
        tables,
        positions,
        rotations,
        extras,
        composePages,
        bases,
        poses,
        groupBoxes,
        indices,
        occludedBits,
      ]) {
        c.ref.value.destroy();
      }
      firstCounts.destroy();
      secondCounts.destroy();
      emptyHiZ.destroy();
    },
  };
}
