/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UserTiming } from '@molstar/core/task/util/user-timing';
import { ValueCell } from '@molstar/core/util';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import { Quat } from '@molstar/core/math/linear-algebra/3d/quat';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import type { Frustum3D } from '@molstar/core/math/geometry/primitives/frustum3d';
import type { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import type { ComputeContext, ComputeKernel, StorageBuffer } from '../compute.js';
import type { GraphicsResources } from '../context.js';
import { MaxCullLevels } from '../compute/cell-cull.js';
import {
  createPagedCuller,
  PagedCullArgsKernel,
  PagedCullArgsSchema,
  type PagedCullChanges,
  PagedCullKernels,
  PagedPerParticleMembers,
  PagedRecordWords,
} from '../compute/paged-cull.js';
import { IndirectRecordBytes, type MultiDrawBaseData } from '../render-item.js';
import type { ItemIndexPattern } from '../quad-indices.js';
import { type LodLevelsValue, LodView } from '../renderable.js';
import { CulledPagedFlag, CulledSegmentVec4s, getCulledSegmentData } from '../culled-data.js';
import { ClipStorageVec4s, getCullClipCount, getInstanceClipCount, packStorageClipObjects } from '../clip-storage.js';
import type { DataArray } from './util.js';
import type { CulledValues, Merged, MergeRecord } from './merged.js';
import { isRasterSize, type MergedCuller, type MergedRasterInput, type MergedRasterParams } from './merged-culler.js';

// The culler of a merged object of a page store (`MergedInstanceSource.paged`): its instances'
// poses composed on the GPU as their pages are written, culled each frame into regions of their store indices by
// member and level (`paged-cull.ts`), a draw record a level and segment whose instance count and first instance a
// kernel writes from the regions (`paged-cull-args`). The poses a page store's `CulledData` (`CulledPagedFlag`: its
// header and segments, then a pose a store index), what the culled programs (`dCulledPaged`, a pose by the index at
// `instance_index`) and the sphere raster read, as the merged culler's. The records' blocks: the levels',
// the raster's (the levels' again, with its split), the ghost's. Nothing of the instances on the CPU: a view of its
// own culled on the GPU too (`Renderable.prepareCullView`), into the same buffers (the next frame's view culled again).

const EmptyInt = new Int32Array(0);
const EmptyUint = new Uint32Array(0);

/**
 * a template's spacing in CSS pixels below which a copy draws it at a coarse level (a level coarser each halving: its
 * sample's spacing twice its level's before), `PagedCullView.coarse`
 */
export const PagedCoarsePixels = 4;

/** the view's pixels a unit for a copy's coarse level (`LodView`: perspective by depth, orthographic everywhere) */
function getCoarse(lod: LodView) {
  if (lod.depthScale > 0) return { scale: 1 / lod.depthScale, ortho: false, pixels: PagedCoarsePixels };
  if (lod.pixelSize > 0) return { scale: 1 / lod.pixelSize, ortho: true, pixels: PagedCoarsePixels };
  return undefined;
}

const dispatchItems = (kernel: ComputeKernel, n: number) => {
  kernel.dispatch(Math.max(1, Math.ceil(n / 64)));
};

export function createPagedMergedCuller(
  compute: ComputeContext,
  resources: GraphicsResources,
  merged: Merged,
  members: readonly MergeRecord[],
  indexed: boolean,
  items?: ItemIndexPattern,
): MergedCuller {
  const { values, segments, source } = merged;
  const paged = source!.paged!;
  const culler = createPagedCuller(compute, resources);
  const stride = indexed ? IndirectRecordBytes.indexed : IndirectRecordBytes.vertex;
  const words = stride / 4;
  // a storage buffer's size a multiple of 4: of 16 (a vec4's) where that is within the device's binding limit, else of
  // 4 (the limit 2^31 - 4: rounded up to 16, 2^31, refused); a growth's 1.5 times capped at the limit by its caller
  const createBuffer = (bytes: number, label: string) => {
    const aligned = Math.ceil(bytes / 16) * 16;
    const size = aligned <= compute.maxStorageBufferSize ? aligned : Math.ceil(bytes / 4) * 4;
    return compute.createBuffer(Math.max(16, size), label);
  };

  const args = ValueCell.create(createBuffer(16, 'paged draw records'));
  /** the second phase's records (`dispatchSecond`): its own regions, drawn after the first's */
  let secondArgs = createBuffer(16, 'paged draw records, second phase');
  let secondPending = false;
  let secondDispatched = false;
  /** the poses, a page store's `CulledData` (its segments in front), and the regions' store indices */
  const poses = ValueCell.create(culler.poses);
  const indices = ValueCell.create(culler.indices);
  /** its copies' markers (`MergedInstanceSource.paged.copyMarkers`) as last uploaded */
  const copyMarkers = ValueCell.create(createBuffer(16, 'paged copy markers'));
  let copyMarkersVersion = -1;
  /** its copies' own layers (`MergedInstanceSource.paged.copyLayers`) as last uploaded */
  const copyLayers = ValueCell.create(createBuffer(16, 'paged copy layers'));
  let copyLayersVersion = -1;
  const culledValues: CulledValues = {
    sCulledIndices: indices,
    sCulledData: poses,
    dCulled: ValueCell.create(true),
    dCulledPaged: ValueCell.create(true),
    sCopyMarkers: copyMarkers,
    dCopyMarkers: ValueCell.create(false),
    sCopyLayers: copyLayers,
    dCopyLayers: ValueCell.create(false),
  };
  /** its copies' own layers uploaded whole where they changed (a state's), the buffer grown to half again */
  function syncCopyLayers() {
    const l = paged.copyLayers;
    if (!l || l.version === copyLayersVersion) return;
    copyLayersVersion = l.version;
    ValueCell.updateIfChanged(culledValues.dCopyLayers, l.buffer.length > 0);
    if (l.buffer.length === 0) return;
    if (copyLayers.ref.value.byteLength < l.buffer.byteLength) {
      copyLayers.ref.value.destroy();
      ValueCell.update(copyLayers, createBuffer(Math.ceil(l.buffer.byteLength * 1.5), 'paged copy layers'));
    }
    copyLayers.ref.value.write(l.buffer);
  }
  /**
   * the copies' markers uploaded where they changed (a mark between frames: before the frame's draws, its cull), the
   * buffer grown to half again; read where the source has copies
   */
  function syncCopyMarkers() {
    const m = paged.copyMarkers;
    if (!m || m.version === copyMarkersVersion) return;
    const first = copyMarkersVersion < 0;
    copyMarkersVersion = m.version;
    // the blocks' bytes changed (all where the buffer is new, the first sync's)
    let ranges = m.takeChanged();
    ValueCell.updateIfChanged(culledValues.dCopyMarkers, m.words.length > 0);
    if (m.words.length === 0) return;
    if (copyMarkers.ref.value.byteLength < m.byteLength) {
      copyMarkers.ref.value.destroy();
      ValueCell.update(copyMarkers, createBuffer(Math.ceil(m.byteLength * 1.5), 'paged copy markers'));
      ranges = undefined;
    }
    if (first) ranges = undefined;
    const buffer = copyMarkers.ref.value;
    buffer.write(m.words);
    // whole words (the blocks' bytes a multiple of four)
    const blocks = Math.min(Math.ceil((m.byteLength - m.words.byteLength) / 4) * 4, m.bytes.length);
    if (blocks === 0) return;
    if (!ranges) {
      buffer.write(new Uint8Array(m.bytes.buffer, 0, blocks), m.words.byteLength);
      return;
    }
    for (let k = 0; k < ranges.length; k += 2) {
      const start = Math.floor(ranges[k] / 4) * 4;
      const end = Math.min(Math.ceil(ranges[k + 1] / 4) * 4, blocks);
      if (end > start) buffer.write(new Uint8Array(m.bytes.buffer, start, end - start), m.words.byteLength + start);
    }
  }
  const argsValues = {
    uRecordCount: ValueCell.create(0),
    uRecordWords: ValueCell.create(words),
    uSegmentCount: ValueCell.create(0),
    uLevels: ValueCell.create(1),
    uMemberCount: ValueCell.create(0),
    uSecond: ValueCell.create(0),
    uIndexRoom: ValueCell.create(0),
    uRasterLevels: ValueCell.create(0),
    uRoomsBase: ValueCell.create(0),
    counts: ValueCell.create(culler.counts),
    bases: ValueCell.create(culler.bases),
    prior: ValueCell.create(culler.secondCounts),
    args,
  };
  let argsKernel: ComputeKernel | undefined;
  let preparing = false;
  const kernelsReady = () => {
    if (compute.hasKernels(PagedCullKernels)) return true;
    if (!preparing) {
      preparing = true;
      compute.prepareKernels(PagedCullKernels).finally(() => (preparing = false));
    }
    return false;
  };

  let dataVersion = 0;
  /** the culler's frame version (`PagedCuller.frameVersion`) as the last dispatch saw it */
  let seenFrameVersion = 0;
  let dispatchVersion = 0;
  let levelCount = 0;
  let segmentCount = 0;
  /** the records' words as written; per segment its member's changes and its offset as read (`layoutRecords`) */
  let writtenRecords = new Uint32Array(0);
  /**
   * the records' template written (`layoutRecords`), and as the second phase's buffer last had it: written there only
   * where it changed (the kernel writes a record's instance count and first alone), not a phase (a frame while moving)
   */
  let templateVersion = 0;
  let secondTemplate = -1;
  /** the members whose records changed since the second phase's buffer was last written (all: -1 to `segmentCount`) */
  let secondLo = 0;
  let secondHi = -1;
  let secondAll = true;
  let seenChanges = new Uint32Array(0);
  let seenOffsets = new Float64Array(0);
  let seenBlocks = -1;
  /** the view of the last dispatch and the data's version: a still view culls nothing again (its records stand) */
  // the frustum, the camera plane, the level ranges, the rest, the raster's split, the transient transform, the view's
  // pixels a unit, the budget
  const dispatchedKey = new Float64Array(6 * 4 + 4 + 16 + 4 + 3 + 16 + 2 + 1).fill(NaN);
  const dispatchKey = new Float64Array(dispatchedKey.length);
  let segmentDataVersion = -1;
  let frontSegments = -1;

  /**
   * whether the records have the raster's blocks (from the first dispatch with its split) and the ghost's (from the
   * first with the ghost on), kept after
   */
  let rasterRegion = false;
  let rasterDispatched = false;
  let ghostRegion = false;
  let ghostDispatched = false;
  /** the blocks of records (a segment's each): the levels', the raster's, the ghost's */
  const levelBlocks = () => Math.max(1, levelCount);
  const ghostBlock = () => levelBlocks() * (rasterRegion ? 2 : 1);
  const recordBlocks = () => ghostBlock() + (ghostRegion ? 1 : 0);
  /** the level of block `k` (a raster block its level's, the ghost's the coarsest) */
  const blockLevel = (k: number) => {
    const levels = levelBlocks();
    return k < levels ? k : k < ghostBlock() ? k - levels : levels - 1;
  };

  function getLevelCount() {
    const lodLevels = members[0].lodLevels?.ref.value as LodLevelsValue | undefined;
    return lodLevels ? Math.min(lodLevels.length, MaxCullLevels) : 0;
  }

  /**
   * The draw records' words but the instance count and first instance (the kernel's): a level's draw count (a
   * member's levels' counts, its draws without), its first index or vertex; written where they changed
   */
  function layoutRecords() {
    const blocks = recordBlocks();
    const size = blocks * segmentCount * words;
    const all = size !== writtenRecords.length || seenBlocks !== blocks * 16 + levelCount;
    const records = all ? new Uint32Array(size) : writtenRecords;
    const changes = merged.memberChanges;
    let changed = all;
    // the changed members' span: a block's records of it written, not all (a step changes a few members of thousands)
    let lo = segmentCount;
    let hi = -1;
    for (let s = 0; s < segmentCount; ++s) {
      // a member's row read again where its cells changed (`Merged.memberChanges`) or its offsets moved: thousands of
      // members an update, a values object's reads scattered and slow
      const offset = indexed ? segments.elementOffsets[s] : segments.drawBases[s];
      if (!all && seenChanges[s] === changes[s] && seenOffsets[s] === offset) continue;
      seenChanges[s] = changes[s];
      seenOffsets[s] = offset;
      changed = true;
      if (s < lo) lo = s;
      hi = s;
      const member = members[s];
      const memberLevels = levelCount > 0 ? (member.lodLevels?.ref.value as LodLevelsValue | undefined) : undefined;
      const draws = member.drawCount.ref.value as number;
      for (let k = 0; k < blocks; ++k) {
        const o = (k * segmentCount + s) * words;
        records[o] = memberLevels ? memberLevels[blockLevel(k)][2] : draws;
        if (indexed && items) {
          records[o + 3] = (segments.elementOffsets[s] / 4 / items.pattern.length) * items.verticesPerItem;
        } else if (indexed) {
          records[o + 2] = segments.elementOffsets[s] / 4;
        } else {
          records[o + 2] = segments.drawBases[s];
        }
      }
    }
    if (!changed) return false;
    seenBlocks = blocks * 16 + levelCount;
    let grown = false;
    if (args.ref.value.byteLength < records.byteLength) {
      args.ref.value.destroy();
      ValueCell.update(args, createBuffer(records.byteLength, 'paged draw records'));
      grown = true;
    }
    if (all || grown) args.ref.value.write(records);
    else writeSpan(args.ref.value, records, blocks, lo, hi);
    // the members whose records changed, for the chunk capacity: all, or a span (`setChunkCapacity`)
    laidOut = all ? undefined : [lo, hi];
    writtenRecords = records;
    templateVersion += 1;
    if (all) secondAll = true;
    else if (!secondAll) {
      secondLo = Math.min(secondLo, lo);
      secondHi = Math.max(secondHi, hi);
    }
    return true;
  }

  /** `records`' members `lo` to `hi` of each of its `blocks` (level-major: a block's members together) written */
  function writeSpan(
    buffer: { write(data: Uint32Array, offset: number): void },
    records: Uint32Array,
    blocks: number,
    lo: number,
    hi: number,
  ) {
    if (hi < lo) return;
    for (let k = 0; k < blocks; ++k) {
      const from = (k * segmentCount + lo) * words;
      const to = (k * segmentCount + hi + 1) * words;
      buffer.write(records.subarray(from, to), from * 4);
    }
  }

  /**
   * The poses' front (`CulledData` of a page store): its header (no instance arrays, the slots' first after the
   * segments, `CulledPagedFlag`) and the segments' data
   */
  function writePosesFront() {
    const tSegment = (values as { tSegment?: ValueCell<DataArray<Float32Array>> }).tSegment;
    const version = tSegment?.ref.version ?? -1;
    if (version === segmentDataVersion && frontSegments === segmentCount) return false;
    segmentDataVersion = version;
    frontSegments = segmentCount;
    const data = getCulledSegmentData(segments, tSegment?.ref.value);
    const front = new Uint32Array(4 + segmentCount * CulledSegmentVec4s * 4);
    front.set([0, segmentCount, 1 + segmentCount * CulledSegmentVec4s, CulledPagedFlag]);
    front.set(new Uint32Array(data.buffer, data.byteOffset, Math.min(data.length, front.length - 4)), 4);
    culler.setPosesFront(front);
    ValueCell.updateIfChanged(poses, culler.poses);
    return true;
  }

  /**
   * per member its raster row (`rasterized`): its largest sphere (`padding`) at each level's scale, read again where
   * the member changed (`Merged.memberChanges`: thousands of members, a values object's reads slow, 2.7 ms the yeast's
   * all), written where a row changed
   */
  let rasterRows = new Float32Array(0);
  let rowChanges = new Uint32Array(0);
  function writeRasterRows() {
    const changes = merged.memberChanges;
    const all = rasterRows.length !== Math.max(1, segmentCount) * 4 || changes.length < segmentCount;
    const rows = all ? new Float32Array(Math.max(1, segmentCount) * 4) : rasterRows;
    if (all) rowChanges = new Uint32Array(segmentCount);
    let changed = all;
    for (let s = 0; s < segmentCount; ++s) {
      if (!all && rowChanges[s] === changes[s]) continue;
      rowChanges[s] = changes[s];
      const member = members[s] as MergeRecord & { padding?: ValueCell<number> };
      const padding = member.padding?.ref.value ?? 0;
      const memberLevels = member.lodLevels?.ref.value as LodLevelsValue | undefined;
      for (let j = 0; j < 4; ++j) {
        const r = padding * (memberLevels?.[j]?.[3] ?? 1);
        if (rows[s * 4 + j] !== r) changed = true;
        rows[s * 4 + j] = r;
      }
    }
    if (!changed) return false;
    rasterRows = rows;
    culler.setRasterRows(rows);
    return true;
  }

  /** chunks of 64 spheres of all raster regions full (each member's instances at each level, `MergedRasterInput`) */
  let chunkCapacity = 0;
  /** the last `layoutRecords`' members whose records it wrote: a span, or all (undefined) */
  let laidOut: [number, number] | undefined;
  /** per member its chunks (`setChunkCapacity`): the capacity kept by the members that changed */
  let memberChunks = new Float64Array(0);
  /**
   * the raster's chunks at most (per member and level its instances times its chunks); `changed`: the members whose
   * counts changed alone (a step's: a few of thousands), else all (the records laid out again, other members)
   */
  function setChunkCapacity(changed?: Iterable<number>) {
    const counts = paged.source.pool.counts;
    const n = Math.min(segmentCount, paged.source.pool.count);
    const levels = levelBlocks();
    const chunksOf = (s: number) => {
      let c = 0;
      for (let j = 0; j < levels; ++j)
        c += Math.ceil(writtenRecords[(j * segmentCount + s) * words] / 6 / 64) * counts[s];
      return c;
    };
    if (!changed || memberChunks.length !== n) {
      memberChunks = new Float64Array(n);
      chunkCapacity = 0;
      for (let s = 0; s < n; ++s) chunkCapacity += memberChunks[s] = chunksOf(s);
      return;
    }
    for (const s of changed) {
      if (s < 0 || s >= n) continue;
      const c = chunksOf(s);
      chunkCapacity += c - memberChunks[s];
      memberChunks[s] = c;
    }
  }

  /**
   * per member its sphere (center, radius: four a member): the box of the pages its runs are in, padded by its reach;
   * none (radius -1) without instances. What a dispatch predicts its records by: a draw a level and member costs
   * whether it draws or not, thousands of members, most out of view. Kept by what changed: per page its members (as
   * last read), per member its box and its pages' count; a page written grows its members' boxes, a page let go leaves
   * them as they are (larger than they need be: a draw more predicted, never one less) until many are, then all boxes
   * again from the pages' (a pass over the pages' members). Not a pass over every page (of a page of a member a
   * particle, its particles) each sync: 1.4 s of a Chrome load
   */
  let memberSpheres = new Float64Array(0);
  let memberBoxes = new Float64Array(0);
  let memberPageCounts = new Int32Array(0);
  let pageMembers: (Int32Array | undefined)[] = [];
  let spheresSource: unknown;
  /** members whose box is larger than their pages' since the last pass over all */
  let staleBoxes = 0;
  /** the pages' members, a pair each: what a pass over all costs */
  let pagePairs = 0;
  /**
   * per member the read of a page that last listed it (`stampCount`): a page's members once each without a set a page
   * (the yeast drawn whole: tens of thousands of pages at once)
   */
  let memberStamps = new Int32Array(0);
  let stampCount = 0;
  let pageList = new Int32Array(0);
  /** the members of page `p` as its record and runs (or extras) say, once each; none for a free page */
  function getPageMembers(p: number): Int32Array | undefined {
    const src = paged.source;
    const r = src.records;
    const o = p * PagedRecordWords;
    const count = r[o + 6];
    if (count === 0) return undefined;
    if (memberStamps.length < src.pool.count) {
      const stamps = new Int32Array(Math.max(src.pool.count, memberStamps.length * 2));
      stamps.set(memberStamps);
      memberStamps = stamps;
    }
    if (stampCount === 0x7fffffff) {
      memberStamps.fill(0);
      stampCount = 0;
    }
    const stamps = memberStamps;
    const stamp = ++stampCount;
    let n = 0;
    const add = (m: number) => {
      if (m >= 0 && m < stamps.length) {
        if (stamps[m] === stamp) return;
        stamps[m] = stamp;
      }
      if (n === pageList.length) {
        const grown = new Int32Array(Math.max(64, n * 2));
        grown.set(pageList);
        pageList = grown;
      }
      pageList[n++] = m;
    };
    if (r[o + 5] === PagedPerParticleMembers) {
      const e = r[o + 11];
      const extras = src.extras;
      let last = -1;
      for (let j = 0; j < count; ++j) {
        const m = extras[(e + j) * 2];
        if (m !== last) add(m);
        last = m;
      }
    } else {
      // inline, not `add` a run: a tile-ordered template's 2.3 M runs
      const runs = src.runs;
      const k0 = r[o + 4];
      const kn = r[o + 5];
      if (pageList.length < n + kn) {
        const grown = new Int32Array(Math.max(64, (n + kn) * 2));
        grown.set(pageList.subarray(0, n));
        pageList = grown;
      }
      for (let k = 0; k < kn; ++k) {
        const m = runs[(k0 + k) * 4 + 2];
        if (m >= 0 && m < stamps.length) {
          if (stamps[m] === stamp) continue;
          stamps[m] = stamp;
        }
        pageList[n++] = m;
      }
    }
    return pageList.slice(0, n);
  }
  /** the records as floats (a page's box), made again where the records are another array */
  let recordFloats: Float32Array = new Float32Array(0);
  function getRecordFloats() {
    const r = paged.source.records;
    if (
      recordFloats.buffer !== r.buffer ||
      recordFloats.byteOffset !== r.byteOffset ||
      recordFloats.length !== r.length
    ) {
      recordFloats = new Float32Array(r.buffer, r.byteOffset, r.length);
    }
    return recordFloats;
  }
  /**
   * per template page (a nested list's) the box its copies place it in: its template's pages' box (in the template's
   * frame, as its records hold it) placed by each copy, their union. Why: a member's sphere of a template's page's box
   * is where no copy need be, a near copy's records not predicted (its spheres through the near plane, not the
   * raster's, not drawn at all)
   */
  let placedBoxes = new Map<number, Float64Array>();
  /** the template pages `placedBoxes` was made of */
  let placedPages = new Set<number>();
  /**
   * a placed box's id (never one again), and per member the id of the placed box its box was last grown by (-1 none,
   * reset with its box): a template's pages share one box, a member of many of them grown by it once, not once a page
   */
  const placedIds = new Map<Float64Array, number>();
  let placedIdCount = 0;
  let memberPlaced = new Int32Array(0);
  function setPlacedBoxes() {
    placedBoxes = new Map();
    placedPages = new Set();
    placedIds.clear();
    const src = paged.source;
    const r = src.records;
    const f = getRecordFloats();
    const q = Quat();
    const v = Vec3();
    for (const t of src.templates ?? []) {
      for (let i = 0; i < t.pages.length; ++i) placedPages.add(t.pages[i]);
      const lo = [Infinity, Infinity, Infinity];
      const hi = [-Infinity, -Infinity, -Infinity];
      for (let i = 0; i < t.pages.length; ++i) {
        const o = t.pages[i] * PagedRecordWords;
        if (r[o + 6] === 0) continue;
        for (let a = 0; a < 3; ++a) {
          if (f[o + 8 + a] < lo[a]) lo[a] = f[o + 8 + a];
          if (f[o + 12 + a] > hi[a]) hi[a] = f[o + 12 + a];
        }
      }
      if (!(hi[0] >= lo[0])) continue;
      const box = new Float64Array([Infinity, Infinity, Infinity, -Infinity, -Infinity, -Infinity]);
      const poses = t.copies.poses;
      // as the kernel places a point (`copyPoint`: t + R v), a corner at a time
      for (let k = 0, n = poses.length / 8; k < n; ++k) {
        const o = k * 8;
        Quat.set(q, poses[o + 4], poses[o + 5], poses[o + 6], poses[o + 7]);
        for (let c = 0; c < 8; ++c) {
          Vec3.set(v, c & 1 ? hi[0] : lo[0], c & 2 ? hi[1] : lo[1], c & 4 ? hi[2] : lo[2]);
          Vec3.transformQuat(v, v, q);
          for (let a = 0; a < 3; ++a) {
            const x = v[a] + poses[o + a];
            if (x < box[a]) box[a] = x;
            if (x > box[3 + a]) box[3 + a] = x;
          }
        }
      }
      for (let i = 0; i < t.pages.length; ++i) placedBoxes.set(t.pages[i], box);
      placedIds.set(box, placedIdCount++);
    }
  }
  /** whether a sync's changes move a template page's placed box (its copies, or a template page written or let go) */
  function changesPlacement(changes: PagedCullChanges) {
    if (changes.copies && changes.copies.length > 0) return true;
    if (placedPages.size === 0 && !paged.source.templates?.length) return false;
    const templatePages = new Set<number>();
    for (const t of paged.source.templates ?? [])
      for (let i = 0; i < t.pages.length; ++i) templatePages.add(t.pages[i]);
    for (const p of changes.written) if (templatePages.has(p) || placedPages.has(p)) return true;
    for (const p of changes.released) if (templatePages.has(p) || placedPages.has(p)) return true;
    return false;
  }
  /** the boxes of page `p`'s members (those in range) grown by its box (a template page's as placed), read once */
  function growBoxes(members: Int32Array, p: number) {
    const f = getRecordFloats();
    const o = p * PagedRecordWords;
    const placed = placedBoxes.get(p);
    const x0 = placed ? placed[0] : f[o + 8];
    const y0 = placed ? placed[1] : f[o + 9];
    const z0 = placed ? placed[2] : f[o + 10];
    const x1 = placed ? placed[3] : f[o + 12];
    const y1 = placed ? placed[4] : f[o + 13];
    const z1 = placed ? placed[5] : f[o + 14];
    const box = memberBoxes;
    const id = placed ? placedIds.get(placed)! : -1;
    for (let i = 0; i < members.length; ++i) {
      const m = members[i];
      if (m < 0 || m >= segmentCount) continue;
      if (placed) {
        if (memberPlaced[m] === id) continue;
        memberPlaced[m] = id;
      }
      const b = m * 6;
      if (x0 < box[b]) box[b] = x0;
      if (y0 < box[b + 1]) box[b + 1] = y0;
      if (z0 < box[b + 2]) box[b + 2] = z0;
      if (x1 > box[b + 3]) box[b + 3] = x1;
      if (y1 > box[b + 4]) box[b + 4] = y1;
      if (z1 > box[b + 5]) box[b + 5] = z1;
    }
  }
  /** the members whose sphere is to be set again, once each (a flag a member and their list: a pair a page's member) */
  let affectedFlags = new Uint8Array(0);
  let affectedList = new Int32Array(0);
  let affectedCount = 0;
  function markAffected(m: number) {
    if (affectedFlags[m]) return;
    affectedFlags[m] = 1;
    affectedList[affectedCount++] = m;
  }
  function emptyBox(m: number) {
    const b = m * 6;
    memberBoxes.fill(Infinity, b, b + 3);
    memberBoxes.fill(-Infinity, b + 3, b + 6);
    memberPlaced[m] = -1;
  }
  /** member `m`'s sphere from its box and reach */
  function setMemberSphere(m: number) {
    const o4 = m * 4;
    const b = m * 6;
    const box = memberBoxes;
    if (memberPageCounts[m] === 0 || !(box[b + 3] >= box[b])) {
      memberSpheres.fill(-1, o4, o4 + 4);
      return;
    }
    memberSpheres[o4] = (box[b] + box[b + 3]) / 2;
    memberSpheres[o4 + 1] = (box[b + 1] + box[b + 4]) / 2;
    memberSpheres[o4 + 2] = (box[b + 2] + box[b + 5]) / 2;
    memberSpheres[o4 + 3] =
      Math.hypot(box[b + 3] - box[b], box[b + 4] - box[b + 1], box[b + 5] - box[b + 2]) / 2 + (paged.reaches[m] ?? 0);
  }
  function setMemberSpheres(changes: PagedCullChanges) {
    const src = paged.source;
    if (affectedFlags.length !== segmentCount) {
      affectedFlags = new Uint8Array(segmentCount);
      affectedList = new Int32Array(segmentCount);
    }
    affectedCount = 0;
    const inRange = (m: number) => m >= 0 && m < segmentCount;
    const setPage = (p: number) => {
      const before = pageMembers[p];
      if (before) {
        for (let i = 0; i < before.length; ++i) {
          const m = before[i];
          if (!inRange(m)) continue;
          memberPageCounts[m] -= 1;
          pagePairs -= 1;
          if (memberPageCounts[m] === 0) emptyBox(m);
          else staleBoxes += 1;
          markAffected(m);
        }
      }
      const after = getPageMembers(p);
      pageMembers[p] = after;
      if (after) {
        growBoxes(after, p);
        for (let i = 0; i < after.length; ++i) {
          const m = after[i];
          if (!inRange(m)) continue;
          memberPageCounts[m] += 1;
          pagePairs += 1;
          markAffected(m);
        }
      }
    };
    // a template page's placed box another: all again (its members' boxes grown by the old one)
    const placement = changesPlacement(changes);
    if (spheresSource !== src || memberSpheres.length !== segmentCount * 4 || placement) {
      spheresSource = src;
      setPlacedBoxes();
      memberSpheres = new Float64Array(segmentCount * 4).fill(-1);
      memberBoxes = new Float64Array(segmentCount * 6);
      memberPlaced = new Int32Array(segmentCount);
      memberPageCounts = new Int32Array(segmentCount);
      for (let m = 0; m < segmentCount; ++m) emptyBox(m);
      pageMembers = [];
      staleBoxes = 0;
      pagePairs = 0;
      for (let p = 0; p < src.pageCount; ++p) setPage(p);
    } else {
      for (const p of changes.released) setPage(p);
      for (const p of changes.written) setPage(p);
      if (changes.members) for (const m of changes.members) if (inRange(m)) markAffected(m);
    }
    // as many boxes larger than their pages' as there are pairs: all again from the pages' members (their records'
    // boxes), its cost a pair a removal
    if (staleBoxes > Math.max(1024, pagePairs)) {
      UserTiming.span('cull: spheres again', () => {
        staleBoxes = 0;
        for (let m = 0; m < segmentCount; ++m) emptyBox(m);
        for (let p = 0; p < pageMembers.length; ++p) {
          const members = pageMembers[p];
          if (members) growBoxes(members, p);
        }
        for (let m = 0; m < segmentCount; ++m) setMemberSphere(m);
      });
    } else {
      for (let i = 0; i < affectedCount; ++i) setMemberSphere(affectedList[i]);
    }
    for (let i = 0; i < affectedCount; ++i) affectedFlags[affectedList[i]] = 0;
    affectedCount = 0;
  }

  /**
   * per block and segment whether it may draw: its member's sphere in the frustum and the level's range (the raster's
   * blocks as their levels' where the split is on, the ghost's in the frustum); a level's impostors none where the
   * raster takes all of it (`raster`: its member's largest sphere, by the largest scale uploaded, raster-sized at the
   * level's nearest distance, as the kernel's `rasterized` decides per instance: conservative). Why: inside a cell
   * every member's sphere is in view and the near levels' ranges, its records all predicted and nearly all empty
   */
  let predicted = new Uint8Array(0);
  function predictRecords(
    cameraPlane: Plane3D,
    frustum: Frustum3D,
    raster: MergedRasterParams | undefined,
    objectTransform: Mat4,
  ) {
    const levels = levelBlocks();
    const blocks = recordBlocks();
    if (predicted.length !== blocks * segmentCount) predicted = new Uint8Array(blocks * segmentCount);
    predicted.fill(0);
    const sp = memberSpheres;
    // as the kernel tests them: moved by the transient transform
    const t = objectTransform;
    const objectScale = Mat4.getMaxScaleOnAxis(t);
    for (let s = 0; s < segmentCount; ++s) {
      const o = s * 4;
      if (!(sp[o + 3] >= 0)) continue;
      const r = sp[o + 3] * objectScale;
      const x = t[0] * sp[o] + t[4] * sp[o + 1] + t[8] * sp[o + 2] + t[12];
      const y = t[1] * sp[o] + t[5] * sp[o + 1] + t[9] * sp[o + 2] + t[13];
      const z = t[2] * sp[o] + t[6] * sp[o + 1] + t[10] * sp[o + 2] + t[14];

      let inside = true;
      for (let i = 0; i < 6 && inside; ++i) {
        const { normal: n, constant } = frustum[i as 0];
        inside = n[0] * x + n[1] * y + n[2] * z + constant >= -r;
      }
      if (!inside) continue;
      if (ghostRegion) predicted[ghostBlock() * segmentCount + s] = 1;
      const n = cameraPlane.normal;
      const d = n[0] * x + n[1] * y + n[2] * z + cameraPlane.constant;
      for (let j = 0; j < levels; ++j) {
        if (levelCount > 0 && (d + r < levelRanges[j * 4] || d - r >= levelRanges[j * 4 + 1])) continue;
        if (raster) {
          predicted[(levels + j) * segmentCount + s] = 1;
          // its instances' nearest distance at the level, their largest sphere
          const nearest = Math.max(d - r, levelCount > 0 ? levelRanges[j * 4] : -Infinity);
          const largest = (rasterRows[s * 4 + Math.min(j, 3)] ?? 0) * culler.maxScale * objectScale * 1.001;
          if (largest > 0 && nearest - largest > raster.near && isRasterSize(raster, largest, nearest)) continue;
        }
        predicted[j * segmentCount + s] = 1;
      }
    }
  }
  /** the lists' records: the predicted ones (`counts`: those the last read of the records counted instances in) */
  function setListRecords(counts: Uint32Array | undefined, of: MultiDrawBaseData[], firstBlock: number) {
    if (predicted.length !== recordBlocks() * segmentCount) return;
    recordsVersion += 1;
    for (let i = 0; i < of.length; ++i) {
      const l = of[i];
      const j = firstBlock + i;
      const records = new Uint32Array(segmentCount);
      let n = 0;
      const base = (j * segmentCount * stride) / 4;
      for (let s = 0; s < segmentCount; ++s) {
        if (!predicted[j * segmentCount + s]) continue;
        if (counts && counts[base + (s * stride) / 4 + 1] === 0) continue;
        records[n++] = s;
      }
      l.count = n;
      l.indirect = { ...l.indirect!, records: records.subarray(0, n), recordsVersion };
    }
  }
  /** every list's records (`first`, `second`: the reads of the records, the second phase's) */
  function setAllRecords(first?: Uint32Array, second?: Uint32Array) {
    const levels = levelBlocks();
    setListRecords(first, lists, 0);
    setListRecords(first, rasterLists, levels);
    setListRecords(first, ghostLists, ghostBlock());
    if (!first || second) {
      setListRecords(second, secondLists, 0);
      setListRecords(second, secondRasterLists, levels);
    }
  }
  /** a still view's records read once: the lists then of those with instances alone */
  let recordsGeneration = 0;
  let recordsRead = -1;
  let recordsCompact = false;
  function compactRecords() {
    if (recordsRead === recordsGeneration || recordsCompact || lists.length === 0) return;
    recordsRead = recordsGeneration;
    const generation = recordsGeneration;
    args.ref.value
      .readAsync(0, recordBlocks() * segmentCount * stride)
      .then(async (data) => {
        if (generation !== recordsGeneration) return;
        const second =
          secondDispatched && secondLists.length > 0 ? await secondArgs.readAsync(0, data.byteLength) : undefined;
        if (generation !== recordsGeneration) return;
        setAllRecords(new Uint32Array(data), second ? new Uint32Array(second) : undefined);
        recordsCompact = true;
      })
      // the predicted lists stay: the same image
      .catch(() => {});
  }

  /** the merged values' cull version as `sync` last read it */
  let syncedCullVersion = -1;
  function sync() {
    UserTiming.span('cull: sync', syncPaged);
    syncCopyMarkers();
    syncCopyLayers();
  }
  function syncPaged() {
    const changes = paged.takeChanges();
    UserTiming.span('cull: upload', () => culler.sync(paged.source, changes, paged.reaches));
    const counts = segmentCount !== segments.count || levelCount !== getLevelCount();
    segmentCount = segments.count;
    levelCount = getLevelCount();
    if (seenChanges.length !== segmentCount) {
      seenChanges = new Uint32Array(segmentCount);
      seenOffsets = new Float64Array(segmentCount).fill(NaN);
      writtenRecords = new Uint32Array(0);
    }
    const front = writePosesFront();
    const rows = writeRasterRows();
    const records = layoutRecords();
    if (records) recordsVersion += 1;
    // the members whose counts changed, and those whose records did (a span where not all)
    let changedChunks: Iterable<number> | undefined = counts ? undefined : changes.members;
    if (changedChunks && records) {
      if (!laidOut) changedChunks = undefined;
      else {
        const span = new Set(changedChunks);
        for (let m = laidOut[0]; m <= laidOut[1]; ++m) span.add(m);
        changedChunks = span;
      }
    }
    setChunkCapacity(changedChunks);
    UserTiming.span('cull: spheres', () => setMemberSpheres(changes));
    const members = changes.members;
    const cull = merged.cullVersion.ref.value !== syncedCullVersion;
    syncedCullVersion = merged.cullVersion.ref.value;
    // what the cull reads changed: its records cull again (`dataVersion`); not where nothing did (an update of the
    // markers alone: the frame's cull stands, the marking pass drawing its lists, `merged.ts`'s update)
    const stored =
      changes.written.length > 0 ||
      changes.released.length > 0 ||
      changes.runs.length > 0 ||
      changes.extras.length > 0 ||
      changes.grown ||
      (members?.size ?? 0) > 0;
    if (stored || counts || front || rows || records || cull) dataVersion += 1;
  }

  /** the levels' distance ranges for `lod` (`LodView.getRange`), per level a column */
  const levelRanges = Mat4.zero();
  const range: number[] = [0, 0, 0];
  function setLevelRanges(lod: LodView) {
    Mat4.setZero(levelRanges);
    const lodLevels = members[0].lodLevels?.ref.value as LodLevelsValue | undefined;
    for (let j = 0; j < levelCount; ++j) {
      LodView.getRange(lodLevels!, j, lod, range);
      levelRanges[j * 4] = range[0];
      levelRanges[j * 4 + 1] = range[1];
    }
  }

  /** the lists of `levels`' uniforms over blocks from `firstBlock` of `buffer` */
  const blockLists = (levels: readonly MultiDrawBaseData[], buffer: StorageBuffer, firstBlock: number) =>
    levels.slice(0, levelBlocks()).map(
      (level, j): MultiDrawBaseData => ({
        firsts: EmptyInt,
        counts: EmptyInt,
        offsets: EmptyInt,
        instanceCounts: EmptyInt,
        baseVertices: EmptyInt,
        baseInstances: EmptyUint,
        count: segmentCount,
        uniforms: level.uniforms,
        indirect: { buffer, byteOffset: (firstBlock + j) * segmentCount * stride, indexed },
      }),
    );

  /** the levels' lists and the raster's (the impostors drawing its instances where the raster does not) */
  let lists: MultiDrawBaseData[] = [];
  let rasterLists: MultiDrawBaseData[] = [];
  let listsLevels: readonly MultiDrawBaseData[] | undefined;
  let listsKey = '';
  function buildLists(levels: readonly MultiDrawBaseData[]) {
    const key = `${segmentCount},${levelCount},${args.ref.value.id},${rasterRegion}`;
    if (levels === listsLevels && key === listsKey) return;
    listsLevels = levels;
    listsKey = key;
    recordsVersion += 1;
    lists = blockLists(levels, args.ref.value, 0);
    rasterLists = rasterRegion ? blockLists(levels, args.ref.value, levelBlocks()) : [];
    setListRecords(undefined, lists, 0);
    setListRecords(undefined, rasterLists, levelBlocks());
    combined.clear();
  }
  /** the second phase's lists: over its records (its instances after the first's in the same regions) */
  let secondLists: MultiDrawBaseData[] = [];
  let secondRasterLists: MultiDrawBaseData[] = [];
  let secondKey = '';
  function buildSecondLists(levels: readonly MultiDrawBaseData[]) {
    const key = `${secondArgs.id},${listsKey}`;
    if (key === secondKey && secondLists.length > 0) return;
    secondKey = key;
    secondLists = blockLists(levels, secondArgs, 0);
    secondRasterLists = rasterRegion ? blockLists(levels, secondArgs, levelBlocks()) : [];
    setListRecords(undefined, secondLists, 0);
    setListRecords(undefined, secondRasterLists, levelBlocks());
    combined.clear();
  }
  /** the lists drawn, by what they are of (the same array while none is made again) */
  const combined = new Map<string, MultiDrawBaseData[]>();
  function getCombined(key: string, parts: () => MultiDrawBaseData[][]) {
    let c = combined.get(key);
    if (!c) {
      c = parts().flat();
      combined.set(key, c);
    }
    return c;
  }

  /** the ghost's list: its block at the coarsest level, at every distance (the coarsest level list's uniforms, its range) */
  let ghostLists: MultiDrawBaseData[] = [];
  let ghostKey = '';
  let ghostLevels: readonly MultiDrawBaseData[] | undefined;
  const ghostLod = ValueCell.create(Vec4());
  const ghostLodLevel = ValueCell.create(0);
  function getGhostListsOf(levels: readonly MultiDrawBaseData[]) {
    const key = `${segmentCount},${levelCount},${args.ref.value.id},${ghostRegion},${rasterRegion}`;
    if (key === ghostKey && levels === ghostLevels) return ghostLists;
    ghostKey = key;
    ghostLevels = levels;
    const j = levelBlocks() - 1;
    const level = levels[Math.min(j, levels.length - 1)];
    const lodLevels = members[0].lodLevels?.ref.value as LodLevelsValue | undefined;
    let uniforms: MultiDrawBaseData['uniforms'] = level?.uniforms ?? [];
    if (levelCount > 0 && lodLevels && level?.uniforms.length) {
      ValueCell.update(ghostLod, Vec4.set(ghostLod.ref.value, 0, Number.MAX_VALUE, 0, lodLevels[j][3]));
      ValueCell.updateIfChanged(ghostLodLevel, Math.min(j, 3));
      uniforms = level.uniforms.map(([name, cell]) =>
        name === 'uLod' ? ['uLod', ghostLod] : name === 'uLodLevel' ? ['uLodLevel', ghostLodLevel] : [name, cell],
      );
    }
    ghostLists = [
      {
        firsts: EmptyInt,
        counts: EmptyInt,
        offsets: EmptyInt,
        instanceCounts: EmptyInt,
        baseVertices: EmptyInt,
        baseInstances: EmptyUint,
        count: segmentCount,
        uniforms,
        indirect: { buffer: args.ref.value, byteOffset: ghostBlock() * segmentCount * stride, indexed },
      },
    ];
    setListRecords(undefined, ghostLists, ghostBlock());
    return ghostLists;
  }

  /**
   * The clip objects cull the instances they clip whole (`getCullClipCount`), their centers with instance clipping
   * (`getInstanceClipCount`); where they changed set, true
   */
  const clipVersions = [-1, -1, -1, -1];
  function syncClip() {
    const count = getCullClipCount(values);
    const frame = values.uClipObjectFrame.ref.version;
    const shape = values.uClipObjectShape.ref.version;
    const instance = getInstanceClipCount(values) > 0 ? 1 : 0;
    const v = clipVersions;
    if (v[0] === count && v[1] === frame && v[2] === shape && v[3] === instance) return false;
    v[0] = count;
    v[1] = frame;
    v[2] = shape;
    v[3] = instance;
    const data = new Float32Array(Math.max(1, count) * ClipStorageVec4s * 4);
    if (count > 0) packStorageClipObjects(values, data);
    culler.setClip({ data, count, instance: instance > 0 });
    return true;
  }

  /** a block added (the raster's, the ghost's): the records laid out again, the lists made again */
  function addBlocks(add: () => void) {
    add();
    layoutRecords();
    setChunkCapacity();
    recordsVersion += 1;
    listsKey = '';
    ghostKey = '';
    secondKey = '';
  }

  let recordsVersion = 0;

  // laid out at once, as the merged culler (`update` syncs it after)
  sync();

  const rasterInput = (): MergedRasterInput => {
    const levels = levelBlocks();
    return {
      indices: indices.ref.value,
      data: poses.ref.value,
      records: args.ref.value,
      // read when the second phase is drawn: its dispatch grows them (the buffer of the first's let go)
      get secondRecords() {
        return secondArgs;
      },
      get secondDispatched() {
        return secondDispatched;
      },
      recordWords: words,
      levels,
      levelCount,
      lodLevels: members[0].lodLevels?.ref.value as LodLevelsValue | undefined,
      // the raster regions' slots (their rooms: the chunks' first capacity)
      instanceCount: culler.rasterRoom,
      segmentCount,
      chunkCapacity,
      dispatchVersion,
      recordBase: levels * segmentCount,
      slotBase: 0,
    };
  };

  return {
    paged: true,
    values: culledValues,
    identityBase: 0,
    setIdentityOrder: () => {},
    sync,
    get dataVersion() {
      return dataVersion;
    },
    areKernelsReady: kernelsReady,
    dispatch(
      cameraPlane,
      frustum,
      lod,
      pyramid,
      raster?: MergedRasterParams,
      _batch?,
      ghost = false,
      twoPhase = false,
    ) {
      syncCopyMarkers();
      syncCopyLayers();
      if (segmentCount === 0 || !kernelsReady()) return false;
      argsKernel ??= compute.createKernel(PagedCullArgsKernel, PagedCullArgsSchema, argsValues);
      // the raster's and the ghost's blocks laid out where first asked for (kept after)
      if (raster && !rasterRegion) addBlocks(() => (rasterRegion = true));
      if (ghost && !ghostRegion) addBlocks(() => (ghostRegion = true));
      // clipped instances culled: a change of the clip objects culls again
      if (syncClip()) dataVersion += 1;
      // a cull's instances past the frame's slots (not drawn): culled again (the room grown, or with a pyramid)
      if (culler.frameVersion !== seenFrameVersion) {
        seenFrameVersion = culler.frameVersion;
        dataVersion += 1;
      }
      setLevelRanges(lod);
      // the same view of the same data: the last dispatch's records stand
      const key = dispatchKey;
      for (let i = 0; i < 6; ++i) {
        const { normal, constant } = frustum[i as 0];
        key[i * 4] = normal[0];
        key[i * 4 + 1] = normal[1];
        key[i * 4 + 2] = normal[2];
        key[i * 4 + 3] = constant;
      }
      key[24] = cameraPlane.normal[0];
      key[25] = cameraPlane.normal[1];
      key[26] = cameraPlane.normal[2];
      key[27] = cameraPlane.constant;
      key.set(levelRanges, 28);
      key[44] = levelCount;
      key[45] = dataVersion;
      // a view culled without a pyramid culled again with one (an app's first frames: none built yet)
      key[46] = pyramid ? 1 : 0;
      key[47] = ghost ? 1 : 0;
      key[48] = raster?.radius ?? 0;
      key[49] = raster?.scale ?? 0;
      key[50] = raster?.near ?? 0;
      // the object's transient transform (a drag's preview): what is drawn moves, culled again
      const objectTransform = values.uObjectTransform.ref.value;
      key.set(objectTransform, 51);
      // a copy's coarse level by the view's pixels a unit (`LodView`)
      key[67] = lod.depthScale;
      key[68] = lod.pixelSize;
      // the copies' levels by a budget (`PagedCullView.maxInstances`): a new one culls again
      const maxInstances = paged.maxInstances ?? 0;
      key[69] = maxInstances;
      let same = true;
      for (let i = 0; same && i < key.length; ++i) same = key[i] === dispatchedKey[i];
      if (same) {
        compactRecords();
        return true;
      }
      dispatchedKey.set(key);
      recordsGeneration += 1;
      recordsCompact = false;
      dispatchVersion += 1;
      predictRecords(cameraPlane, frustum, rasterRegion ? raster : undefined, objectTransform);
      setAllRecords();
      // the last cull's second phase not run (another view first): its cursors read, its rooms grown
      if (secondPending) culler.endFirst();
      culler.cull(
        {
          frustum,
          cameraPlane,
          levels: levelRanges,
          levelCount,
          ghost,
          rasterRegions: rasterRegion,
          raster,
          objectTransform,
          pyramidToCome: twoPhase && !pyramid,
          secondPhase: twoPhase,
          coarse: getCoarse(lod),
          maxInstances,
        },
        pyramid,
      );
      ghostDispatched = ghost;
      rasterDispatched = !!raster;
      // its buffers as they are now (grown by the cull)
      ValueCell.updateIfChanged(poses, culler.poses);
      ValueCell.updateIfChanged(indices, culler.indices);
      ValueCell.updateIfChanged(argsValues.counts, culler.counts);
      ValueCell.updateIfChanged(argsValues.bases, culler.bases);
      ValueCell.updateIfChanged(argsValues.prior, culler.secondCounts);
      const levels = levelBlocks();
      const recordCount = recordBlocks() * segmentCount;
      ValueCell.updateIfChanged(argsValues.uRecordCount, recordCount);
      ValueCell.updateIfChanged(argsValues.uSegmentCount, segmentCount);
      ValueCell.updateIfChanged(argsValues.uLevels, levels);
      ValueCell.updateIfChanged(argsValues.uMemberCount, culler.culled ? culler.memberCount : 0);
      ValueCell.updateIfChanged(argsValues.uRasterLevels, rasterRegion ? levels : 0);
      ValueCell.updateIfChanged(argsValues.uIndexRoom, culler.indexRoom);
      ValueCell.updateIfChanged(argsValues.uRoomsBase, culler.regionCount * 2);
      compute.beginComputePass({ label: 'paged draw records' });
      dispatchItems(argsKernel, recordCount);
      compute.endComputePass();
      secondPending = culler.culled && !!pyramid;
      secondDispatched = false;
      combined.clear();
      return true;
    },
    get ghostDispatched() {
      return ghostDispatched;
    },
    getGhostLists: (levels) => getGhostListsOf(levels),
    get secondPending() {
      return secondPending;
    },
    dispatchSecond(pyramid) {
      if (!secondPending) return false;
      secondPending = false;
      if (!pyramid || !argsKernel || !culler.cullSecond(pyramid)) {
        // none to run: the first's cursors read (its rooms grown)
        culler.endFirst();
        return false;
      }
      // the second's records: its regions' counts and bases into records of the same draws
      const bytes = writtenRecords.byteLength;
      if (secondArgs.byteLength < bytes) {
        secondArgs.destroy();
        secondArgs = createBuffer(bytes, 'paged draw records, second phase');
        secondKey = '';
        secondTemplate = -1;
      }
      if (secondTemplate !== templateVersion) {
        if (secondTemplate < 0 || secondAll) secondArgs.write(writtenRecords);
        else writeSpan(secondArgs, writtenRecords, recordBlocks(), secondLo, secondHi);
        secondTemplate = templateVersion;
        secondAll = false;
        secondLo = segmentCount;
        secondHi = -1;
      }
      // the second's records: its cursors, its regions' firsts after the first's instances (`prior`)
      const first = { args: args.ref.value };
      ValueCell.update(argsValues.counts, culler.secondCounts);
      ValueCell.update(argsValues.prior, culler.counts);
      ValueCell.update(args, secondArgs);
      ValueCell.updateIfChanged(argsValues.uSecond, 1);
      compute.beginComputePass({ label: 'paged draw records, second phase' });
      dispatchItems(argsKernel, argsValues.uRecordCount.ref.value);
      compute.endComputePass();
      ValueCell.updateIfChanged(argsValues.uSecond, 0);
      ValueCell.update(argsValues.counts, culler.counts);
      ValueCell.update(argsValues.prior, culler.secondCounts);
      ValueCell.update(args, first.args);
      secondDispatched = true;
      combined.clear();
      return true;
    },
    get secondDispatched() {
      return secondDispatched;
    },
    getSecondLists(levels, withRaster = true) {
      buildLists(levels);
      buildSecondLists(levels);
      const raster = withRaster && rasterRegion;
      return raster ? getCombined('second+raster', () => [secondLists, secondRasterLists]) : secondLists;
    },
    getLists(levels, withRaster = true) {
      buildLists(levels);
      const raster = withRaster && rasterRegion;
      // after a second phase both (a still view's next frames draw what it found too)
      if (!secondDispatched) return raster ? getCombined('raster', () => [lists, rasterLists]) : lists;
      buildSecondLists(levels);
      return raster
        ? getCombined('both+raster', () => [lists, rasterLists, secondLists, secondRasterLists])
        : getCombined('both', () => [lists, secondLists]);
    },
    getRasterInput: () => (rasterRegion && rasterDispatched ? rasterInput() : undefined),
    get recordsCompact() {
      return recordsCompact;
    },
    get recordsVersion() {
      return recordsVersion;
    },
    destroy() {
      copyMarkers.ref.value.destroy();
      copyLayers.ref.value.destroy();
      argsKernel?.destroy();
      culler.destroy();
      args.ref.value.destroy();
      secondArgs.destroy();
    },
  };
}
