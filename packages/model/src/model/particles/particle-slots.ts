/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomProperties } from '../custom-property.js';
import type { ModelFormat } from '../../formats/format.js';
import type { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import type { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import {
  concatAttributes,
  concatPositions,
  concatRotations,
  concatScales,
  concatValues,
  type Segment,
} from './particle-concat.js';
import {
  Particle,
  type ParticleCompartmentInfo,
  type ParticleEntityInfo,
  type ParticleList,
  type ParticleTarget,
  type ParticleTargetInfo,
  ParticleValues,
} from './particle-list.js';

// Slot lists: a particle's index its slot, kept while it lives, the slots in pages of any size, a page a range of them
// with a flat list of its own (a tile's: no copy); the slots no page holds gaps. What reads every particle walks the
// pages (`forEachPage`), each a flat list read as one; a lookup a slot only where one is wanted (picking, a label). A
// version says what changed since the one before (`pagesFrom`), what consumers update by.

/** A page of a slot list: the slots `[start, start + list.count)`, its particles `list`'s */
export interface ParticlePage {
  readonly start: number;
  /** a flat list (no pages of its own) */
  readonly list: ParticleList;
}

export interface ParticlePages {
  /** by `start`, none overlapping */
  readonly pages: readonly ParticlePage[];
  /** the pages' starts, ascending: a slot's page found by a search */
  readonly starts: Int32Array;
  /** the slots the pages hold */
  readonly liveCount: number;
  /** the page found last (`findPage`): a slot of the same page, as a loop over slots in order finds, without a search */
  last: number;
  /**
   * the sphere the particles of all slots are within, whatever pages are live (a tile set's): its consumers' bounds
   * from the first version on, not grown page by page
   */
  readonly sphere?: Sphere3D;
}

/**
 * What changed from the version whose pages are `from` to the version that holds this: the pages added and removed (by
 * identity), and those posed otherwise (a page at the same slots of the same particles, `Particle.isSamePoses`: a frame
 * step; the new page). By the earlier version's pages, not its list: a list would hold its own `pagesFrom`, every
 * version since the first reachable from the last (their removed pages' particles never collected)
 */
export interface ParticlePagesChange {
  readonly from: ParticlePages;
  readonly added: readonly ParticlePage[];
  readonly removed: readonly ParticlePage[];
  readonly posed: readonly ParticlePage[];
}

/** What a slot list shares across its pages (a tile set's catalog: one object each for all) */
export interface ParticleSlotsInfo {
  readonly label?: string;
  readonly entryId?: string;
  readonly cell?: Cell;
  readonly targetInfo: ReadonlyMap<number, ParticleTargetInfo>;
  readonly targetMapping?: ReadonlyMap<number, ParticleTarget>;
  readonly entityInfo?: ReadonlyMap<number, ParticleEntityInfo>;
  readonly compartmentInfo?: ReadonlyMap<number, ParticleCompartmentInfo>;
  readonly sourceData: ModelFormat;
}

/** A slot and its particle: its page's list and its index in it */
export interface ParticleSlot {
  list: ParticleList;
  index: number;
}

export namespace ParticleSlots {
  export function isPaged(list: ParticleList): boolean {
    return !!list.pages;
  }

  /** The particles: a slot list's its pages hold, a flat list's all */
  export function getLiveCount(list: ParticleList): number {
    return list.pages ? list.pages.liveCount : list.count;
  }

  /** `f` of each page, by slot: a flat list one page of all its particles */
  export function forEachPage(list: ParticleList, f: (page: ParticlePage, index: number) => void) {
    if (!list.pages) {
      f({ start: 0, list }, 0);
      return;
    }
    const { pages } = list.pages;
    for (let i = 0, il = pages.length; i < il; ++i) f(pages[i], i);
  }

  /** The page holding `slot`, -1 for a gap (a flat list: 0 where `slot` is one of its particles) */
  export function findPage(list: ParticleList, slot: number): number {
    if (!list.pages) return slot >= 0 && slot < list.count ? 0 : -1;
    const cache = list.pages;
    const { pages, starts } = cache;
    const last = pages[cache.last];
    if (last && slot >= last.start && slot < last.start + last.list.count) return cache.last;
    let lo = 0;
    let hi = starts.length - 1;
    if (hi < 0 || slot < starts[0]) return -1;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (starts[mid] <= slot) lo = mid;
      else hi = mid - 1;
    }
    if (slot >= pages[lo].start + pages[lo].list.count) return -1;
    cache.last = lo;
    return lo;
  }

  /** `slot`'s particle into `out` (its page's list, its index there), undefined for a gap */
  export function locate(list: ParticleList, slot: number, out: ParticleSlot = { list, index: 0 }) {
    if (!list.pages) {
      if (slot < 0 || slot >= list.count) return undefined;
      out.list = list;
      out.index = slot;
      return out;
    }
    const p = findPage(list, slot);
    if (p < 0) return undefined;
    const page = list.pages.pages[p];
    out.list = page.list;
    out.index = slot - page.start;
    return out;
  }

  /**
   * A slot list of `count` slots, its particles `pages`' (by start; none overlapping, each within the slots, each a
   * flat list), what they share `info`'s. Given `from`, a slot list of the same slots this is the next version of:
   * the pages added and removed between them (`getChanges`), by identity. `sphere`: what all slots' particles are
   * within (`ParticlePages.sphere`), `from`'s by default
   */
  export function create(
    info: ParticleSlotsInfo,
    count: number,
    pages: readonly ParticlePage[],
    from?: ParticleList,
    sphere: Sphere3D | undefined = from?.pages?.sphere,
  ): ParticleList {
    const sorted = [...pages].sort((a, b) => a.start - b.start);
    let liveCount = 0;
    for (let i = 0; i < sorted.length; ++i) {
      const { start, list } = sorted[i];
      if (list.pages) throw new Error('particle slots: a page of pages');
      if (start < 0 || start + list.count > count) throw new Error(`particle slots: a page past the slots (${start})`);
      if (i > 0 && start < sorted[i - 1].start + sorted[i - 1].list.count) {
        throw new Error(`particle slots: pages overlapping (${start})`);
      }
      liveCount += list.count;
    }
    const starts = Int32Array.from(sorted, (p) => p.start);
    let pagesFrom: ParticlePagesChange | undefined;
    if (from?.pages) {
      const before = new Set(from.pages.pages);
      const now = new Set(sorted);
      const added = sorted.filter((p) => !before.has(p));
      const removed = from.pages.pages.filter((p) => !now.has(p));
      // a page replaced by one of the same particles posed otherwise: posed, neither added nor removed
      const removedAt = new Map(removed.map((p) => [p.start, p]));
      const posed: ParticlePage[] = [];
      const newPages: ParticlePage[] = [];
      for (const p of added) {
        const old = removedAt.get(p.start);
        if (old && old.list.count === p.list.count && Particle.isSamePoses(old.list, p.list)) {
          posed.push(p);
          removedAt.delete(p.start);
        } else {
          newPages.push(p);
        }
      }
      pagesFrom = { from: from.pages, added: newPages, removed: [...removedAt.values()], posed };
    }
    const slot: ParticleSlot = { list: undefined as unknown as ParticleList, index: 0 };
    const list: ParticleList = {
      ...info,
      count,
      targets: new Int32Array(0),
      coordinates: new Float32Array(0),
      getParticleLabel: (index) => {
        const s = locate(list, index, slot);
        return s ? s.list.getParticleLabel(s.index) : `Slot ${index + 1} (empty)`;
      },
      pages: { pages: sorted, starts, liveCount, last: 0, sphere },
      pagesFrom,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    return list;
  }

  /** The next version of slot list `list`: `add` pages added, those of `remove` (by start) removed */
  export function update(
    list: ParticleList,
    add: readonly ParticlePage[],
    remove: Iterable<number> = [],
  ): ParticleList {
    if (!list.pages) throw new Error('particle slots: not a slot list');
    const gone = new Set(remove);
    const kept = list.pages.pages.filter((p) => !gone.has(p.start));
    return create(list, list.count, [...kept, ...add], list);
  }

  /** What changed from `from` to `to` where `to` is its next version (`pagesFrom`), else undefined: built again */
  export function getChanges(from: ParticleList, to: ParticleList): ParticlePagesChange | undefined {
    return from.pages !== undefined && to.pagesFrom?.from === from.pages ? to.pagesFrom : undefined;
  }

  /**
   * A flat list of a slot list's particles, by slot (copied; a particle's key its slot, its values where every page
   * has them): for what reads flat lists only; a flat list as it is
   */
  export function flatten(list: ParticleList): ParticleList {
    if (!list.pages) return list;
    const segments: Segment[] = [];
    const firstSlots: number[] = [];
    let count = 0;
    for (const page of list.pages.pages) {
      const n = page.list.count;
      segments.push({
        tile: page.start,
        list: page.list,
        indices: Int32Array.from({ length: n }, (_, i) => i),
        start: 0,
        end: n,
      });
      firstSlots.push(count);
      count += n;
    }
    const keys = new Int32Array(count);
    const targets = new Int32Array(count);
    let o = 0;
    for (const page of list.pages.pages) {
      const n = page.list.count;
      for (let i = 0; i < n; ++i, ++o) {
        keys[o] = page.start + i;
        targets[o] = ParticleValues.get(page.list.targets, i);
      }
    }
    const all = (f: (l: ParticleList) => unknown) =>
      segments.length > 0 && segments.every((s) => f(s.list) !== undefined);
    const pages = list.pages.pages;
    return {
      entryId: list.entryId,
      label: list.label,
      count,
      cell: list.cell,
      keys,
      targets: ParticleValues.compact(targets),
      targetInfo: list.targetInfo,
      targetMapping: list.targetMapping,
      entities: all((l) => l.entities) ? concatValues(segments, count, (l) => l.entities!, Int32Array) : undefined,
      entityInfo: list.entityInfo,
      compartments: all((l) => l.compartments)
        ? concatValues(segments, count, (l) => l.compartments!, Int32Array)
        : undefined,
      compartmentInfo: list.compartmentInfo,
      coordinates: concatPositions(segments, count),
      rotations: concatRotations(segments, count),
      radii: all((l) => l.radii) ? concatValues(segments, count, (l) => l.radii!, Float32Array) : undefined,
      // its pages' coarse nodes drawn larger (`ParticleList.scale`, `scales`), as the page store draws them
      scales: concatScales(segments, count),
      attributes: concatAttributes(segments, count),
      getParticleLabel: (index) => {
        let p = 0;
        while (p + 1 < firstSlots.length && firstSlots[p + 1] <= index) p += 1;
        return pages[p].list.getParticleLabel(index - firstSlots[p]);
      },
      sourceData: list.sourceData,
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
  }
}
