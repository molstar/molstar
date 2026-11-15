/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
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
import type { CulledValues, Merged, MergeableValues } from './merged.js';
import type { MergedCuller, MergedRasterInput, MergedRasterParams } from './merged-culler.js';

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

const dispatchItems = (kernel: ComputeKernel, n: number) => {
  kernel.dispatch(Math.max(1, Math.ceil(n / 64)));
};

export function createPagedMergedCuller(
  compute: ComputeContext,
  resources: GraphicsResources,
  merged: Merged,
  members: readonly MergeableValues[],
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
  const culledValues: CulledValues = {
    sCulledIndices: indices,
    sCulledData: poses,
    dCulled: ValueCell.create(true),
    dCulledPaged: ValueCell.create(true),
  };
  const argsValues = {
    uRecordCount: ValueCell.create(0),
    uRecordWords: ValueCell.create(words),
    uSegmentCount: ValueCell.create(0),
    uLevels: ValueCell.create(1),
    uMemberCount: ValueCell.create(0),
    uSecond: ValueCell.create(0),
    uIndexRoom: ValueCell.create(0),
    uRasterLevels: ValueCell.create(0),
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
  let dispatchVersion = 0;
  let levelCount = 0;
  let segmentCount = 0;
  /** the records' words as written; per segment its member's changes and its offset as read (`layoutRecords`) */
  let writtenRecords = new Uint32Array(0);
  let seenChanges = new Uint32Array(0);
  let seenOffsets = new Float64Array(0);
  let seenBlocks = -1;
  /** the view of the last dispatch and the data's version: a still view culls nothing again (its records stand) */
  // the frustum, the camera plane, the level ranges, the rest, the raster's split, the transient transform
  const dispatchedKey = new Float64Array(6 * 4 + 4 + 16 + 4 + 3 + 16).fill(NaN);
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
    for (let s = 0; s < segmentCount; ++s) {
      // a member's row read again where its cells changed (`Merged.memberChanges`) or its offsets moved: thousands of
      // members an update, a values object's reads scattered and slow
      const offset = indexed ? segments.elementOffsets[s] : segments.drawBases[s];
      if (!all && seenChanges[s] === changes[s] && seenOffsets[s] === offset) continue;
      seenChanges[s] = changes[s];
      seenOffsets[s] = offset;
      changed = true;
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
    if (args.ref.value.byteLength < records.byteLength) {
      args.ref.value.destroy();
      ValueCell.update(args, createBuffer(records.byteLength, 'paged draw records'));
    }
    args.ref.value.write(records);
    writtenRecords = records;
    return true;
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
      const member = members[s] as MergeableValues & { padding?: ValueCell<number> };
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
  function setChunkCapacity() {
    const counts = paged.source.pool.counts;
    const n = Math.min(segmentCount, paged.source.pool.count);
    chunkCapacity = 0;
    for (let s = 0; s < n; ++s) {
      for (let j = 0; j < levelBlocks(); ++j) {
        chunkCapacity += Math.ceil(writtenRecords[(j * segmentCount + s) * words] / 6 / 64) * counts[s];
      }
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
      const runs = src.runs;
      for (let k = 0, k0 = r[o + 4], kn = r[o + 5]; k < kn; ++k) add(runs[(k0 + k) * 4 + 2]);
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
  /** the boxes of page `p`'s members (those in range) grown by its box, read once */
  function growBoxes(members: Int32Array, p: number) {
    const f = getRecordFloats();
    const o = p * PagedRecordWords;
    const x0 = f[o + 8];
    const y0 = f[o + 9];
    const z0 = f[o + 10];
    const x1 = f[o + 12];
    const y1 = f[o + 13];
    const z1 = f[o + 14];
    const box = memberBoxes;
    for (let i = 0; i < members.length; ++i) {
      const m = members[i];
      if (m < 0 || m >= segmentCount) continue;
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
    if (spheresSource !== src || memberSpheres.length !== segmentCount * 4) {
      spheresSource = src;
      memberSpheres = new Float64Array(segmentCount * 4).fill(-1);
      memberBoxes = new Float64Array(segmentCount * 6);
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
      staleBoxes = 0;
      for (let m = 0; m < segmentCount; ++m) emptyBox(m);
      for (let p = 0; p < pageMembers.length; ++p) {
        const members = pageMembers[p];
        if (members) growBoxes(members, p);
      }
      for (let m = 0; m < segmentCount; ++m) setMemberSphere(m);
    } else {
      for (let i = 0; i < affectedCount; ++i) setMemberSphere(affectedList[i]);
    }
    for (let i = 0; i < affectedCount; ++i) affectedFlags[affectedList[i]] = 0;
    affectedCount = 0;
  }

  /**
   * per block and segment whether it may draw: its member's sphere in the frustum and the level's range (the raster's
   * blocks as their levels' where the split is on, the ghost's in the frustum)
   */
  let predicted = new Uint8Array(0);
  function predictRecords(cameraPlane: Plane3D, frustum: Frustum3D, raster: boolean, objectTransform: Mat4) {
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
        predicted[j * segmentCount + s] = 1;
        if (raster) predicted[(levels + j) * segmentCount + s] = 1;
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
    const changes = paged.takeChanges();
    culler.sync(paged.source, changes, paged.reaches);
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
    setChunkCapacity();
    setMemberSpheres(changes);
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
      secondRecords: secondArgs,
      get secondDispatched() {
        return secondDispatched;
      },
      recordWords: words,
      levels,
      levelCount,
      lodLevels: members[0].lodLevels?.ref.value as LodLevelsValue | undefined,
      // a block of regions a level the store's instances (its raster region's slots in it, down from a member's end)
      instanceCount: culler.instanceTotal,
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
    dispatch(cameraPlane, frustum, lod, pyramid, raster?: MergedRasterParams, _batch?, ghost = false) {
      if (segmentCount === 0 || !kernelsReady()) return false;
      argsKernel ??= compute.createKernel(PagedCullArgsKernel, PagedCullArgsSchema, argsValues);
      // the raster's and the ghost's blocks laid out where first asked for (kept after)
      if (raster && !rasterRegion) addBlocks(() => (rasterRegion = true));
      if (ghost && !ghostRegion) addBlocks(() => (ghostRegion = true));
      // clipped instances culled: a change of the clip objects culls again
      if (syncClip()) dataVersion += 1;
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
      predictRecords(cameraPlane, frustum, rasterRegion && !!raster, objectTransform);
      setAllRecords();
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
      if (!pyramid || !argsKernel || !culler.cullSecond(pyramid)) return false;
      // the second's records: its regions' counts and bases into records of the same draws
      const bytes = writtenRecords.byteLength;
      if (secondArgs.byteLength < bytes) {
        secondArgs.destroy();
        secondArgs = createBuffer(bytes, 'paged draw records, second phase');
        secondKey = '';
      }
      secondArgs.write(writtenRecords);
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
      argsKernel?.destroy();
      culler.destroy();
      args.ref.value.destroy();
      secondArgs.destroy();
    },
  };
}
