/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomProperties, CustomPropertyDescriptor } from '../custom-property.js';
import {
  concatAttributes,
  concatPositions,
  concatRotations,
  concatScales,
  concatValues,
  findSegment,
  type Segment,
} from './particle-concat.js';
import {
  Particle,
  type ParticleList,
  type ParticleTarget,
  type ParticleTargetInfo,
  ParticleValues,
} from './particle-list.js';
import { type ParticlePage, ParticleSlots } from './particle-slots.js';
import { type ParticleTileSet, ParticleTiles } from './particle-tiles.js';

// A view of a tile set: one `ParticleList` of some tiles at a frame, what the representations draw as any list. Its
// particles by target (a target one run: its instances one range), within a target by tile (ascending), within a tile
// in the tile's order; a particle's key its index in the whole (its tile's base and its index in the tile: kept across
// views of other tiles). What a target is made of (its composition: the tiles and their ranges) says which targets'
// instances a view of other tiles changes; a frame step the same composition, the same particles posed otherwise.

/** A target's part of a view: the tiles of its particles (ascending) and where each one's are, `[start, end)` */
export interface ParticleTileViewPart {
  readonly tiles: Int32Array;
  /** per tile its first particle in the view, the end last (tiles + 1 values) */
  readonly offsets: Int32Array;
}

export interface ParticleTileView {
  readonly set: ParticleTileSet;
  readonly frame: number;
  /** the tiles, ascending */
  readonly tiles: Int32Array;
  readonly list: ParticleList;
  /** per target of the view (ascending) its part; none for a slot view (its changes by page instead) */
  readonly composition: ReadonlyMap<number, ParticleTileViewPart>;
  /** a slot view (`createSlots`): its list a slot list, a tile a page at its base */
  readonly slots?: boolean;
  /** a slot view's: its source's `targetsVersion` its mapping is of */
  readonly targetsVersion?: number;
}

/** What gives a view its tiles' lists and targets: a loader's ready ones (`ParticleTileLoader`) */
export interface ParticleTileViewSource {
  getTile(id: number, frame: number): ParticleList | undefined;
  getTarget(id: number): ParticleTarget | undefined;
  /** the targets loaded (a loader's): a slot view's mapping, else each tile's looked up */
  readonly targets?: ReadonlyMap<number, ParticleTarget>;
  /** advanced as `targets` changes: a slot view of the same version keeps its mapping (`createSlots`) */
  readonly targetsVersion?: number;
}

const ViewDescriptor = CustomPropertyDescriptor({ name: 'particle-tile-view' });

/** a tile's particles by target: its indices grouped by target (each in the tile's order), a target's range in them */
interface TileTargets {
  readonly indices: Int32Array;
  readonly ranges: ReadonlyMap<number, readonly [number, number]>;
}

/** per tile list (the first of its line of poses: its frames' the same) its particles by target */
const tileTargets = new WeakMap<ParticleList, TileTargets>();
function getTileTargets(list: ParticleList): TileTargets {
  const key = list.posesFrom ?? list;
  let t = tileTargets.get(key);
  if (t) return t;
  // a counting sort by target: a pack's tile a run a target, so mostly ranges copied
  const counts = new Map<number, number>();
  ParticleValues.forEachRun(list.targets, list.count, (start, end, target) => {
    counts.set(target, (counts.get(target) ?? 0) + end - start);
  });
  const ranges = new Map<number, [number, number]>();
  const fill = new Map<number, number>();
  let o = 0;
  for (const [target, n] of counts) {
    ranges.set(target, [o, o + n]);
    fill.set(target, o);
    o += n;
  }
  const indices = new Int32Array(list.count);
  ParticleValues.forEachRun(list.targets, list.count, (start, end, target) => {
    let k = fill.get(target)!;
    for (let i = start; i < end; ++i) indices[k++] = i;
    fill.set(target, k);
  });
  t = { indices, ranges };
  tileTargets.set(key, t);
  return t;
}

export namespace ParticleTileView {
  /** The view a list is (`create`), undefined for another list */
  export function of(list: ParticleList): ParticleTileView | undefined {
    // a stub list (param definitions made without data) has none
    return list._propertyData?.[ViewDescriptor.name];
  }

  /**
   * The view of `tiles` of `set` at `frame`, their lists and targets from `source` (each must be there: a loader's
   * ready tiles). Given `previous`, a view of the same set and tiles at another frame: the same particles posed
   * otherwise (`Particle.withPoses` of its list: `Particle.isSamePoses` holds, what was built of their identity kept)
   */
  export function create(
    set: ParticleTileSet,
    source: ParticleTileViewSource,
    tiles: Iterable<number>,
    frame: number,
    previous?: ParticleTileView,
  ): ParticleTileView {
    const ids = Int32Array.from(new Set(tiles)).sort();
    const lists = new Map<number, ParticleList>();
    for (const id of ids) {
      const list = source.getTile(id, frame);
      if (!list) throw new Error(`particle tile view: tile ${id} at frame ${frame} not ready`);
      lists.set(id, list);
    }

    // per target (ascending) its tiles' particles
    const byTarget = new Map<number, Segment[]>();
    for (const id of ids) {
      const list = lists.get(id)!;
      const { indices, ranges } = getTileTargets(list);
      for (const [target, [start, end]] of ranges) {
        let segments = byTarget.get(target);
        if (!segments) byTarget.set(target, (segments = []));
        segments.push({ tile: id, list, indices, start, end });
      }
    }
    const targetIds = [...byTarget.keys()].sort((a, b) => a - b);
    const segments: Segment[] = [];
    const composition = new Map<number, ParticleTileViewPart>();
    let count = 0;
    for (const target of targetIds) {
      const parts = byTarget.get(target)!;
      const offsets = new Int32Array(parts.length + 1);
      for (let k = 0; k < parts.length; ++k) {
        offsets[k] = count;
        count += parts[k].end - parts[k].start;
        segments.push(parts[k]);
      }
      offsets[parts.length] = count;
      composition.set(target, { tiles: Int32Array.from(parts, (p) => p.tile), offsets });
    }

    const coordinates = concatPositions(segments, count);
    const rotations = concatRotations(segments, count);

    if (previous && previous.set === set && sameTiles(previous.tiles, ids) && sameComposition(previous, composition)) {
      const list = Particle.withPoses(previous.list, coordinates, rotations);
      const view: ParticleTileView = { set, frame, tiles: ids, list, composition: previous.composition };
      attach(list, view);
      return view;
    }

    const tileOf = set.tiles;
    const keys = new Int32Array(count);
    const segmentStarts = new Int32Array(segments.length + 1);
    let o = 0;
    for (let s = 0; s < segments.length; ++s) {
      const { tile, list, indices, start, end } = segments[s];
      segmentStarts[s] = o;
      const base = tileOf.get(tile)!.base;
      const keysOf = list.keys;
      for (let i = start; i < end; ++i, ++o) {
        // a tile's own keys where it has them (the whole's indices), else its base and the index in it
        keys[o] = keysOf ? keysOf[indices[i]] : base + indices[i];
      }
    }
    segmentStarts[segments.length] = count;

    const targetMapping = new Map<number, ParticleTarget>();
    const targetInfo = new Map<number, ParticleTargetInfo>();
    for (const id of targetIds) {
      const target = source.getTarget(id);
      if (!target) throw new Error(`particle tile view: target ${id} not loaded`);
      targetMapping.set(id, target);
      const info = set.catalog.targetInfo.get(id);
      if (info) targetInfo.set(id, info);
    }

    const first = segments[0]?.list;
    const all = (f: (l: ParticleList) => unknown) =>
      segments.length > 0 && segments.every((s) => f(s.list) !== undefined);
    const list: ParticleList = {
      entryId: first?.entryId,
      label: set.label,
      count,
      cell: set.catalog.cell,
      keys,
      // a target a run (the composition's ranges)
      targets: ParticleValues.ofRuns(
        Int32Array.from([...targetIds.map((t) => composition.get(t)!.offsets[0]), count]),
        Int32Array.from(targetIds),
      ),
      targetInfo,
      targetMapping,
      entities: all((l) => l.entities) ? concatValues(segments, count, (l) => l.entities!, Int32Array) : undefined,
      entityInfo: set.catalog.entityInfo,
      compartments: all((l) => l.compartments)
        ? concatValues(segments, count, (l) => l.compartments!, Int32Array)
        : undefined,
      compartmentInfo: set.catalog.compartmentInfo,
      coordinates,
      rotations,
      radii: all((l) => l.radii) ? concatValues(segments, count, (l) => l.radii!, Float32Array) : undefined,
      // a coarse node's particles drawn larger (`ParticleList.scale`, its runs' `scales`): by particle where a tile's are
      scales: concatScales(segments, count),
      attributes: concatAttributes(segments, count),
      getParticleLabel: (index) => {
        const s = findSegment(segmentStarts, index);
        const segment = segments[s];
        return segment.list.getParticleLabel(segment.indices[segment.start + index - segmentStarts[s]]);
      },
      // an empty view's the set's own
      sourceData: first?.sourceData ?? { kind: 'particle-tiles', name: set.label ?? '', data: set },
      customProperties: new CustomProperties(),
      _propertyData: Object.create(null),
    };
    const view: ParticleTileView = { set, frame, tiles: ids, list, composition };
    attach(list, view);
    return view;
  }

  /** Whether `target`'s part is the same in `a` and `b`: its tiles and their counts (its instances the same particles) */
  export function isSamePart(a: ParticleTileView, b: ParticleTileView, target: number) {
    // a slot view has no composition: its consumers update by its pages' changes
    if (a.slots || b.slots) return false;
    const p = a.composition.get(target);
    const q = b.composition.get(target);
    if (!p || !q) return p === q;
    if (p.tiles.length !== q.tiles.length) return false;
    for (let i = 0; i < p.tiles.length; ++i) {
      if (p.tiles[i] !== q.tiles[i]) return false;
      if (p.offsets[i + 1] - p.offsets[i] !== q.offsets[i + 1] - q.offsets[i]) return false;
    }
    return true;
  }

  /**
   * The slot view of `tiles` of `set` at `frame` (`ParticleSlots`): its slots the set's particles, a tile a page
   * at its base (its list as it is: no copy), a particle's slot its key. Given `previous`, a slot view of the same set:
   * its next version, the pages added, removed and posed between them (a page kept where its tile's list is the same);
   * an update by the tiles' count, not their particles'
   */
  /** per tile-ordered list its slot list of all its tiles (`ofWhole`) */
  const wholes = new WeakMap<ParticleList, ParticleList>();

  /**
   * Tile-ordered `list` (`ParticleTiles.getDirectory`: a pack written with tiles) as a slot list of all its tiles: a
   * tile a page, a view of its range of `list` (no copy), a slot a particle's index in `list`; drawn as a tile view
   * of all its tiles (by a page store, its targets' particles ranges a tile), not a target an index set. Undefined
   * without a directory; kept per list
   */
  export function ofWhole(list: ParticleList): ParticleList | undefined {
    const directory = ParticleTiles.getDirectory(list);
    if (!directory) return undefined;
    let whole = wholes.get(list);
    if (whole) return whole;
    const d = directory.dimensions;
    const set = ParticleTiles.fromList(list, [d, d, d]);
    const source: ParticleTileViewSource = {
      getTile: (id, frame) => set.provider.getTile(id, frame) as ParticleList,
      getTarget: (id) => list.targetMapping?.get(id),
      targets: list.targetMapping ?? new Map(),
    };
    whole = createSlots(set, source, set.tiles.keys(), 0).list;
    wholes.set(list, whole);
    return whole;
  }

  export function createSlots(
    set: ParticleTileSet,
    source: ParticleTileViewSource,
    tiles: Iterable<number>,
    frame: number,
    previous?: ParticleList,
  ): ParticleTileView {
    const ids = Int32Array.from(new Set(tiles)).sort();
    const from = previous && of(previous)?.set === set && previous.pages ? previous : undefined;
    const before = new Map(from?.pages!.pages.map((p) => [p.start, p]));
    const pages: ParticlePage[] = [];
    for (const id of ids) {
      const list = source.getTile(id, frame);
      if (!list) throw new Error(`particle tile view: tile ${id} at frame ${frame} not ready`);
      const start = set.tiles.get(id)!.base;
      const kept = before.get(start);
      pages.push(kept && kept.list === list ? kept : { start, list });
    }
    let targetMapping: Map<number, ParticleTarget>;
    // the last view's mapping where its source's targets are as they were (a step: thousands of targets not copied
    // again, the representation's check of them by identity)
    const last = from ? of(from) : undefined;
    const version = source.targetsVersion;
    if (
      source.targets &&
      last?.targetsVersion !== undefined &&
      last.targetsVersion === version &&
      from?.targetMapping
    ) {
      targetMapping = from.targetMapping as Map<number, ParticleTarget>;
    } else if (source.targets) {
      targetMapping = new Map(source.targets);
    } else {
      targetMapping = new Map();
      for (const id of ids) {
        for (const t of set.tiles.get(id)!.targets) {
          if (targetMapping.has(t)) continue;
          const target = source.getTarget(t);
          if (!target) throw new Error(`particle tile view: target ${t} not loaded`);
          targetMapping.set(t, target);
        }
      }
    }
    const { catalog } = set;
    const list = ParticleSlots.create(
      {
        label: set.label,
        cell: catalog.cell,
        targetInfo: catalog.targetInfo,
        targetMapping,
        entityInfo: catalog.entityInfo,
        compartmentInfo: catalog.compartmentInfo,
        sourceData: pages[0]?.list.sourceData ?? { kind: 'particle-tiles', name: set.label ?? '', data: set },
      },
      set.count,
      pages,
      from,
      set.boundingSphere,
    );
    const view: ParticleTileView = {
      set,
      frame,
      tiles: ids,
      list,
      composition: new Map(),
      slots: true,
      targetsVersion: source.targets ? version : undefined,
    };
    attach(list, view);
    return view;
  }

  function attach(list: ParticleList, view: ParticleTileView) {
    list.customProperties.add(ViewDescriptor);
    list._propertyData[ViewDescriptor.name] = view;
  }

  function sameTiles(a: Int32Array, b: Int32Array) {
    if (a.length !== b.length) return false;
    for (let i = 0; i < a.length; ++i) if (a[i] !== b[i]) return false;
    return true;
  }

  function sameComposition(view: ParticleTileView, composition: ReadonlyMap<number, ParticleTileViewPart>) {
    if (view.composition.size !== composition.size) return false;
    for (const [target, part] of composition) {
      const p = view.composition.get(target);
      if (!p || p.offsets.length !== part.offsets.length) return false;
      for (let i = 0; i < p.offsets.length; ++i) if (p.offsets[i] !== part.offsets[i]) return false;
    }
    return true;
  }
}
