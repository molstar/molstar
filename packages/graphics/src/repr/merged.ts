/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type GraphicsRenderObject,
  type MergedGraphicsRenderObject,
  canMergeRenderObjects,
  createMergedRenderObject,
} from '@molstar/graphics/gl/render-object';
import {
  type GpuCullingValue,
  type RenderMergedValue,
  resolveGpuCulling,
  resolveRenderMerged,
} from '@molstar/graphics/geo/geometry/base';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { type Loci, EmptyLoci } from '@molstar/model/model/loci';
import type { ValueCell } from '@molstar/core/util/value-cell';

/**
 * What `canMergeRenderObjects` reads of a member (its defines, the color/size/layer types, instance granularity and
 * LOD), as the cells and their versions when it last passed. A group of the same members with none of them changed
 * was mergeable then and is now: the check (a sorted key string per member) is skipped (perf-audit-2 07/P17).
 * `uLod` and `lodLevels` are kept by value, as `getMergeKey` reads them: `BaseGeometry.updateValues` sets `uLod` on
 * every update, and a structure update `lodLevels`, with the same values (perf-audit-3 07/P17).
 */
type MergeFingerprint = { keyCount: number; cells: ValueCell<any>[]; versions: number[]; lod: string };
const mergeFingerprints = new WeakMap<GraphicsRenderObject, MergeFingerprint>();

/** The `*Type` cells `canMergeValues` reads */
const MergeTypeKeys = new Set([
  'colorType',
  'sizeType',
  'overpaintType',
  'transparencyType',
  'emissiveType',
  'substanceType',
  'clippingType',
  'wiggleType',
]);

function isMergeRelevantKey(k: string) {
  return (
    (k.length > 1 && k[0] === 'd' && k[1] >= 'A' && k[1] <= 'Z') || MergeTypeKeys.has(k) || k === 'uInstanceGranularity'
  );
}

/** `uLod` and `lodLevels` as `getMergeKey` writes them */
function lodOf(values: { [k: string]: ValueCell<any> }) {
  const lod = values.uLod?.ref.value as ArrayLike<number> | undefined;
  const levels = values.lodLevels?.ref.value as ArrayLike<number>[] | undefined;
  let s = lod ? `${lod[0]},${lod[1]},${lod[2]},${lod[3]}` : '';
  if (levels) for (const l of levels) s += `;${l[0]},${l[1]},${l[2]}`;
  return s;
}

function recordMergeFingerprint(ro: GraphicsRenderObject) {
  const values = ro.values as unknown as { [k: string]: ValueCell<any> };
  const keys = Object.keys(values);
  const cells: ValueCell<any>[] = [];
  const versions: number[] = [];
  for (const k of keys) {
    if (!isMergeRelevantKey(k) || !values[k]?.ref) continue;
    cells.push(values[k]);
    versions.push(values[k].ref.version);
  }
  mergeFingerprints.set(ro, { keyCount: keys.length, cells, versions, lod: lodOf(values) });
}

function isMergeFingerprintUnchanged(ro: GraphicsRenderObject) {
  const f = mergeFingerprints.get(ro);
  if (!f) return false;
  const values = ro.values as unknown as { [k: string]: ValueCell<any> };
  if (Object.keys(values).length !== f.keyCount) return false;
  for (let i = 0, il = f.cells.length; i < il; ++i) {
    if (f.cells[i].ref.version !== f.versions[i]) return false;
  }
  return lodOf(values) === f.lod;
}

/** One visual and its renderObject, ready to be grouped/merged by `mergeRenderObjectsByType`. */
export type MergeableEntry<V> = { visual: V; renderObject: GraphicsRenderObject };

/** A currently-merged group of render objects, keyed by their shared geometry type. */
export type MergedEntry<V> = { object: MergedGraphicsRenderObject; visuals: V[]; key: string };

/**
 * Groups `entries` by `renderObject.type` and merges each group of >= 2 mergeable
 * render objects into one shared render item, when `renderMerged` resolves to true
 * for the group's size. Reuses `prevMergedByType`'s entry for a group whose member
 * ids are unchanged, instead of rebuilding a merged render object every update.
 *
 * Clears and refills `renderObjects` with the merged and any left-over unmerged
 * render objects. Returns the new by-type map, to pass as `prevMergedByType` next time.
 */
export function mergeRenderObjectsByType<V>(
  entries: readonly MergeableEntry<V>[],
  renderMerged: RenderMergedValue,
  prevMergedByType: ReadonlyMap<string, MergedEntry<V>>,
  renderObjects: GraphicsRenderObject[],
  /** culled on the GPU (`BaseGeometry.Params.gpuCulling`): a merged object of the other kind is a new one */
  gpuCulling: GpuCullingValue = false,
  /** whether the device culls on the GPU (`canCullOnGpu`): else no group is merged for it */
  gpuSupported = false,
): Map<string, MergedEntry<V>> {
  renderObjects.length = 0;
  const nextMerged = new Map<string, MergedEntry<V>>();
  const consumed = new Set<GraphicsRenderObject>();

  if (renderMerged !== false) {
    const byType = new Map<string, MergeableEntry<V>[]>();
    for (const e of entries) {
      const group = byType.get(e.renderObject.type);
      if (group) group.push(e);
      else byType.set(e.renderObject.type, [e]);
    }
    byType.forEach((group, type) => {
      // culled on the GPU where many instances per draw (a segment and level): merged for it, also one object alone
      let instances = 0;
      for (const e of group) instances += e.renderObject.values.instanceCount.ref.value;
      const levels = (group[0].renderObject.values as { lodLevels?: { ref: { value: unknown[] } } }).lodLevels;
      const draws = group.length * Math.max(1, levels?.ref.value.length ?? 0);
      const gpu = gpuSupported && resolveGpuCulling(gpuCulling, instances, draws);
      const merge = gpu || (group.length >= 2 && resolveRenderMerged(renderMerged, group.length));
      if (!merge) return;

      const key = `${group.map((e) => e.renderObject.id).join(',')}|${gpu}`;
      const prev = prevMergedByType.get(type);
      const same = !!prev && prev.key === key;
      // the members checked again only where one of them changed what merging reads (the key: the same members, the
      // same culling and so the same least count)
      if (!same || !group.every((e) => isMergeFingerprintUnchanged(e.renderObject))) {
        if (
          !canMergeRenderObjects(
            group.map((e) => e.renderObject),
            gpu ? 1 : 2,
          )
        )
          return;
        for (const e of group) recordMergeFingerprint(e.renderObject);
      }

      const entry = same
        ? prev!
        : {
            object: createMergedRenderObject(
              group.map((e) => e.renderObject),
              { gpuCulling: gpu },
            ),
            visuals: group.map((e) => e.visual),
            key,
          };
      nextMerged.set(type, entry);
      renderObjects.push(entry.object);
      for (const e of group) consumed.add(e.renderObject);
    });
  }

  for (const e of entries) {
    if (!consumed.has(e.renderObject)) renderObjects.push(e.renderObject);
  }

  return nextMerged;
}

/**
 * Resolves a `pickingId` that may target a merged render object: finds the
 * `MergedEntry` whose object id matches, translates the global instance id to
 * the owning member's index and local instance id, and delegates to `resolve`.
 * Returns undefined when `pickingId` does not target any currently merged
 * object, letting the caller fall back to per-visual picking.
 */
export function pickMergedLoci<V>(
  mergedByType: ReadonlyMap<string, MergedEntry<V>>,
  pickingId: PickingId,
  resolve: (visual: V, index: number, localInstanceId: number, object: MergedGraphicsRenderObject) => Loci,
): Loci | undefined {
  let loci: Loci | undefined;
  mergedByType.forEach(({ object, visuals }) => {
    if (loci !== undefined || pickingId.objectId !== object.id) return;
    const { segments } = object.merged;
    const { instanceId } = pickingId;
    for (let i = 0, il = segments.count; i < il; ++i) {
      const base = segments.instanceBases[i];
      if (instanceId >= base && instanceId < base + segments.instanceCounts[i]) {
        loci = resolve(visuals[i], i, instanceId - base, object);
        return;
      }
    }
    loci = EmptyLoci;
  });
  return loci;
}
