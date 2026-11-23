/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CategoryColumn } from '@molstar/core/repr/tree';
import type { PairsList, Want } from '@molstar/seqx-core/data';
import type { AccessPlan, PairsData } from '@molstar/seqx-core/resolve';
import {
  dataRegions,
  fittedLevel,
  itemsGate,
  type LayerSource,
  pairsLevels,
  readTiles,
  type TileRead,
} from './layer-data.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import { withGenome } from '@molstar/seqx-core/resolve';
import { sidedName } from '@molstar/seqx-core/data';
import { PairDatum, type PairsPrep } from './pairs-lay.js';
import type { PairsReadOptions } from './pairs-reader.js';

// Pairs read by region (_spec/README.md "Render", The pairs layer: by region): the reader's (`pairs-reader.ts`) over a source's
// pairs, as held pairs are read. Two reads: bins at a pixel's level (`pairs.binPx`), always; items as items tiles at the
// finest level whose tiles fit, asked once the bin pairs with an end in view sum to at most `pairs.maxArcs` and kept
// while the items in view stay within `pairs.itemsHysteresis` × it. Bins draw until the items land (EC-169: pending is
// never empty), each bin pair one arc between its bins' middles, stroked by count. What is read becomes the frame's
// `PairsData` (ends on the layer's two systems, a pair found from its other end turned round; every row of the tiles
// read, each once by the anchor rule), kept while the entries read are the same, so a warm frame builds nothing.

/**
 * Pairs over an assembly's reader (`pairs-assembly.ts`), loaded with the genome chunk: a model with an assembly waits
 * for both, so a layer over one never reads before it lands. Why lazy: a viewer that never opens a genome pays nothing.
 */
const AssemblyPairs: { read?: typeof import('./pairs-assembly.js').assemblyPairs } = {};
withGenome(() =>
  import('./pairs-assembly.js').then((m) => {
    AssemblyPairs.read = m.assemblyPairs;
  }),
);

/** What a region read of pairs keeps on its owner (a lay, a reader) between frames. */
export interface RegionState {
  /** What the data was built from: the mode (0 bins, 1 items), the level, then each entry's tile and version. */
  readonly read: number[];
  prep: PairsPrep | undefined;
  /** Whether items are shown (the gate's state, for its hysteresis). */
  items: boolean;
  /** Bin pairs with an end in view (summed counts), as last read; items in view when shown. */
  binPairs: number;
  itemCount: number;
  /** The bins' tiles not loaded, `[start, end, …]` on the data's system (a bridge's pending hatch, EC-169). */
  pending?: readonly number[];
  /** Whether the last read waits on a bins' tile it wanted (a region's or an assembly's block's). */
  waiting: boolean;
  /**
   * The data's positions read this frame by its system (`dataRegions`: projected, normalized, merged near; over an
   * assembly each block's system): a bridge's drawn-once rule reads them.
   */
  readonly windows: Map<string, number[]>;
}

/** An owner's region state (a pairs lay's, a bridge side's, a reader's), made on first read. */
const States = new WeakMap<object, RegionState>();

export function regionState(owner: object): RegionState {
  let s = States.get(owner);
  if (!s)
    States.set(
      owner,
      (s = {
        read: [],
        prep: undefined,
        items: false,
        binPairs: 0,
        itemCount: 0,
        windows: new Map(),
        waiting: false,
      }),
    );
  return s;
}

/** Where a region pairs layer reads: its source, the names of its two systems as the source knows them. */
export interface PairsRegion {
  readonly src: LayerSource;
  readonly aName: string;
  readonly bName: string;
  /** Over an assembly: each block's plan by its system, and the unit (`PairsPrep.endUnits`) a listed name is. */
  readonly assembly?: AssemblyRegion;
}

/** Pairs over an assembly read by region (`pairs-assembly.ts`): a block a sequence. */
export interface AssemblyRegion {
  /** Each sequence's read (its plan's source and name, its side of a sided source), by its system. */
  readonly sources: ReadonlyMap<string, LayerSource>;
  /** One genome against itself over a sided source: a sequence read on either side (`sourceOn`). */
  readonly turned?: ReadonlyMap<string, LayerSource>;
  /** The sequences' systems: the ends' units (`PairsData.aUnit`, `bUnit`), in the assembly's order. */
  readonly units: readonly string[];
  /**
   * A listed name's unit (index into `units`), by its side and name over a sided source (`sidedName`: `chr1` listed on
   * both sides is two sequences), by its name alone else.
   */
  readonly unitOf: ReadonlyMap<string, number>;
  readonly sided: boolean;
}

/**
 * The region of pairs over an assembly: its sequences' plans by system as each block's reads; `plansB`, one genome
 * against itself over a sided source: each sequence's plan on side `b` (its name there).
 */
export function assemblyRegion(
  plans: ReadonlyMap<string, AccessPlan>,
  lengthOf: (system: string) => number | undefined,
  plansB?: ReadonlyMap<string, AccessPlan>,
): AssemblyRegion {
  const sources = new Map<string, LayerSource>();
  const turned = new Map<string, LayerSource>();
  const units: string[] = [];
  const unitOf = new Map<string, number>();
  let sided = false;
  const sourceOf = (system: string, plan: AccessPlan): LayerSource | undefined => {
    if (plan.kind !== 'source') return undefined;
    const length = lengthOf(system);
    unitOf.set(sidedName(plan.side, plan.system), units.length - 1);
    if (plan.side) sided = true;
    return {
      source: plan.source,
      system: plan.system,
      dataSystem: system,
      shape: 'pairs',
      ...(plan.side ? { side: plan.side } : {}),
      ...(length !== undefined ? { length, wraps: true } : {}),
      ratio: 1,
    };
  };
  for (const [system, plan] of plans) {
    if (plan.kind !== 'source') continue;
    units.push(system);
    sources.set(system, sourceOf(system, plan)!);
    const b = plansB?.get(system);
    const src = b && sourceOf(system, b);
    if (src) turned.set(system, src);
  }
  return { sources, ...(turned.size ? { turned } : {}), units, unitOf, sided };
}

/** A sequence's read on a side (`side`: the side asked; a sequence of one side read there whatever is asked). */
export function sourceOn(over: AssemblyRegion, system: string, side: 'a' | 'b' | undefined): LayerSource | undefined {
  const src = over.sources.get(system);
  return src?.side && side && src.side !== side ? (over.turned?.get(system) ?? src) : src;
}

/**
 * A bridge's read of a side's `R` (_spec/README.md "Render", Bridges: drawn once): the regions by data system, at the
 * near side's own level (`undefined`: items), never coarser than the far level, no gate; `wants`, the layer's others
 * this frame (`R`'s are asked only when they fit beside them).
 */
export interface RegionOver {
  regions: ReadonlyMap<string, readonly number[]>;
  level: number | undefined;
  wants: readonly Want[];
}

/** A dry read's wants, and those with the layer's others (`overFits`): kept, reused. */
const Asked: Want[] = [];
const Load: Want[] = [];
const Dry = {
  want: (w: Want) => {
    Asked.push(w);
  },
  config: undefined as unknown as LayoutContext['config'],
  key: undefined as string | undefined,
};

/** A layout context whose reads only collect their wants (`Asked`): an `R` read is asked once it is known to fit. */
export function dryOf(ctx: LayoutContext): LayoutContext {
  Asked.length = 0;
  Dry.config = ctx.config;
  Dry.key = ctx.key;
  return Dry as unknown as LayoutContext;
}

/**
 * Whether a dry read's wants fit the store's bounds beside the layer's others (`maxTilesPerWant` already said by the
 * read's `tooMany`; the layer's tiles and bytes here): then asked, else never (a layer refused for them would flicker).
 */
export function overFits(ctx: LayoutContext, data: StoreReader, others: readonly Want[]): boolean {
  Load.length = 0;
  for (let i = 0; i < others.length; i++) Load.push(others[i]);
  for (let i = 0; i < Asked.length; i++) Load.push(Asked[i]);
  if (data.layerLoad && !data.layerLoad(Load).fits) return false;
  for (let i = 0; i < Asked.length; i++) ctx.want(Asked[i]);
  return true;
}

/** Whether `[s, e)` meets any of `regions` (`[start, end, …]`). */
function inView(regions: readonly number[], s: number, e: number): boolean {
  for (let i = 0; i + 1 < regions.length; i += 2)
    if (s < regions[i + 1] && Math.max(e, s + 1) > regions[i]) return true;
  return false;
}

/**
 * The frame's pairs of a region layer: reads bins, gates items, and gives the prep the arcs lay out (`undefined` while
 * nothing is read). The read's state is its owner's (`regionState`). Over an assembly, a read a block
 * (`pairs-assembly.ts`; `side`, a bridge side's of one genome against itself). `over`: the regions read instead of the
 * view's at a level given (a bridge's near side reads its `R`); items pending, nothing. `options`: the reader's bin
 * width and gate (the pairs config's when absent).
 */
export function regionPairs(
  prep: PairsPrep,
  region: PairsRegion,
  ctx: LayoutContext,
  data: StoreReader,
  owner: object,
  over?: RegionOver,
  side?: 'a' | 'b',
  options?: PairsReadOptions,
): PairsPrep | undefined {
  if (region.assembly)
    return AssemblyPairs.read?.(prep, region, region.assembly, ctx, data, owner, over, side, options);
  const state = regionState(owner);
  const { src } = region;
  const cov = data.coverage(src.source);
  const levels = cov?.levels ?? [1];
  const regions = over
    ? ((over.regions.get(src.dataSystem) as number[] | undefined) ?? [])
    : dataRegions(src, ctx.axis);
  if (!over) {
    // one system: its entry replaced, no entry made a layout
    if (state.windows.size !== 1 || !state.windows.has(src.dataSystem)) state.windows.clear();
    state.windows.set(src.dataSystem, regions);
  }
  const span = (ctx.axis.u1 - ctx.axis.u0) * src.ratio;
  const cfg = ctx.config.pairs;
  const fit = pairsLevels(
    levels,
    src.ratio / ctx.axis.devicePx,
    span,
    data,
    options?.binPx ?? cfg.binPx,
    options?.gate ?? cfg.maxArcs,
  );
  const hysteresis = options?.hysteresis ?? cfg.itemsHysteresis;
  const items = fit.items;
  if (over) return overRead(prep, region, ctx, data, state, regions, over, items);
  // past the store's bounds on the layer's wants together: coarser bins, then zoom in
  const bins = fittedLevel(ctx, data, src, regions, fit.bins);
  const binRead = readTiles(ctx, data, src, regions, bins);
  state.pending = binRead.pending;
  state.waiting = binRead.pending.length > 0;
  // at the finest level the bins are the items themselves
  const binned = binRead.entries.some((e) => e.tables.pairs?.form === 'bins');
  state.binPairs = binned
    ? countBins(binRead.entries, bins, region, regions)
    : countItems(binRead.entries, region, regions);
  let mode = binned ? 0 : 1;
  let read = binRead.entries;
  if (binned) {
    const ask = itemsGate(state.items, state.binPairs, state.items ? state.itemCount : undefined, fit.gate, hysteresis);
    state.items = false;
    if (ask) {
      const itemRead = readTiles(ctx, data, src, regions, items, { items: true });
      if (!itemRead.pending.length && !itemRead.failed.length) {
        state.itemCount = countItems(itemRead.entries, region, regions);
        state.items = state.itemCount <= fit.gate * Math.max(1, hysteresis);
        if (state.items) {
          mode = 1;
          read = itemRead.entries;
        }
      } else state.items = true;
    }
  }
  const level = mode ? items : bins;
  if (sameRead(state.read, mode, level, read) && state.prep) return state.prep;
  const built = mode ? fromItems(read, region) : fromBins(read, bins, region);
  keepRead(state.read, mode, level, read);
  state.prep = withData(prep, built.data, built.weight, mode ? undefined : bins);
  return state.prep;
}

/**
 * A read of `R` at the near side's level (bins, or items: `undefined`), no gate; its wants asked only when they fit
 * beside the layer's; not whole (pending, failed, past a bound), nothing.
 */
function overRead(
  prep: PairsPrep,
  region: PairsRegion,
  ctx: LayoutContext,
  data: StoreReader,
  state: RegionState,
  regions: readonly number[],
  over: RegionOver,
  items: number,
): PairsPrep | undefined {
  const mode = over.level === undefined ? 1 : 0;
  const at = over.level ?? items;
  const r = readTiles(dryOf(ctx), data, region.src, regions, at, mode ? { items: true } : undefined);
  if (!overFits(ctx, data, over.wants) || r.pending.length || r.failed.length) return undefined;
  if (sameRead(state.read, mode, at, r.entries) && state.prep) return state.prep;
  const built = mode ? fromItems(r.entries, region) : fromBins(r.entries, at, region);
  keepRead(state.read, mode, at, r.entries);
  state.prep = withData(prep, built.data, built.weight, mode ? undefined : at);
  return state.prep;
}

export type Loaded = TileRead['entries'][number];

/** Whether the data was built from this read (its mode, level, tiles and versions): nothing made to tell. */
export function sameRead(r: readonly number[], mode: number, level: number, entries: readonly Loaded[]): boolean {
  if (r.length !== 2 + 2 * entries.length || r[0] !== mode || r[1] !== level) return false;
  for (let k = 0; k < entries.length; k++)
    if (r[2 + 2 * k] !== entries[k].region.tile || r[3 + 2 * k] !== entries[k].version) return false;
  return true;
}

export function keepRead(r: number[], mode: number, level: number, entries: readonly Loaded[]) {
  r.length = 0;
  r.push(mode, level);
  for (const e of entries) r.push(e.region.tile, e.version);
}

/**
 * Whether a row of entry `k` is that entry's to take, by the anchor rule (_spec/README.md "Render", The pairs layer: by
 * region): a tile holds every row with an end on the layer's system in it, so a row is taken by the first loaded tile
 * (by start) that holds it; `as`, `ae` its `a` end when on the system (else NaN), `bs`, `be` its `b` end.
 */
function owned(entries: readonly Loaded[], k: number, as: number, ae: number, bs: number, be: number): boolean {
  const start = entries[k].region.start;
  for (let o = 0; o < entries.length; o++) {
    const r = entries[o].region;
    if (o === k || r.start >= start) continue;
    if (as < r.end && Math.max(ae, as + 1) > r.start) return false;
    if (bs < r.end && Math.max(be, bs + 1) > r.start) return false;
  }
  return true;
}

/** The two systems' names of a row, `undefined` when it carries none (both on the request's). */
export function names(sys: CategoryColumn | undefined, i: number): string | undefined {
  return sys ? (sys.get(i) ?? undefined) : undefined;
}

/** How a row's ends fall on the layer's systems: 1 as they are, −1 turned round, 0 not on them. */
function orient(a: string | undefined, b: string | undefined, region: PairsRegion): number {
  const as = a ?? region.src.system;
  const bs = b ?? region.src.system;
  if (as === region.aName && bs === region.bName) return 1;
  if (as === region.bName && bs === region.aName) return -1;
  return 0;
}

/** The ends a read anchors by: a sided read's side's alone (1 `a`, 2 `b`), else either (3). */
function endsOf(region: PairsRegion): number {
  const side = region.src.side;
  return side === 'a' ? 1 : side === 'b' ? 2 : 3;
}

/** Bin pairs with an end in view, their counts summed, each once across tiles (the anchor rule: no key made). */
function countBins(entries: readonly Loaded[], level: number, region: PairsRegion, regions: readonly number[]): number {
  const sys = region.src.system;
  const ends = endsOf(region);
  let n = 0;
  for (let k = 0; k < entries.length; k++) {
    const t = entries[k].tables.pairs;
    if (t?.form !== 'bins' || 'aStart' in t) continue;
    for (let i = 0; i < t.count.length; i++) {
      const aOn = (ends & 1) !== 0 && (names(t.aSystem, i) ?? sys) === sys;
      const bOn = (ends & 2) !== 0 && (names(t.bSystem, i) ?? sys) === sys;
      const as = aOn ? t.aBin[i] * level : Number.NaN;
      const bs = bOn ? t.bBin[i] * level : Number.NaN;
      if (!(aOn && inView(regions, as, as + level)) && !(bOn && inView(regions, bs, bs + level))) continue;
      if (owned(entries, k, as, as + level, bs, bs + level)) n += t.count[i];
    }
  }
  return n;
}

/** A pair's identity across tiles: its row, else its id (not `.` or empty), with both ends; else its ends alone. */
export function itemKey(t: PairsList, i: number): string {
  const named = t.row ? `r${t.row[i]}` : t.id ? t.id.get(i) : null;
  const ends = `${names(t.aSystem, i)}:${t.aStart[i]}-${t.aEnd[i]}|${names(t.bSystem, i)}:${t.bStart[i]}-${t.bEnd[i]}`;
  return named && named !== '.' ? `${named}|${ends}` : ends;
}

/** Pairs with an end in view, each once across tiles (the anchor rule: no key made). */
function countItems(entries: readonly Loaded[], region: PairsRegion, regions: readonly number[]): number {
  const sys = region.src.system;
  const ends = endsOf(region);
  let n = 0;
  for (let k = 0; k < entries.length; k++) {
    const t = entries[k].tables.pairs;
    if (t?.form !== 'list') continue;
    for (let i = 0; i < t.aStart.length; i++) {
      const aOn = (ends & 1) !== 0 && (names(t.aSystem, i) ?? sys) === sys;
      const bOn = (ends & 2) !== 0 && (names(t.bSystem, i) ?? sys) === sys;
      const aIn = aOn && inView(regions, t.aStart[i], t.aEnd[i]);
      const bIn = bOn && inView(regions, t.bStart[i], t.bEnd[i]);
      if (!aIn && !bIn) continue;
      const as = aOn ? t.aStart[i] : Number.NaN;
      const bs = bOn ? t.bStart[i] : Number.NaN;
      if (owned(entries, k, as, t.aEnd[i], bs, t.bEnd[i])) n++;
    }
  }
  return n;
}

/** Rows gathered into a `PairsData` (and each one's weight). */
export class Gather {
  aStart: number[] = [];
  aEnd: number[] = [];
  bStart: number[] = [];
  bEnd: number[] = [];
  aStrand: number[] = [];
  bStrand: number[] = [];
  kind: number[] = [];
  kinds: string[] = [];
  kindAt = new Map<string, number>();
  keys: (string | number)[] = [];
  weight: number[] = [];
  push(
    as: number,
    ae: number,
    bs: number,
    be: number,
    sa: number,
    sb: number,
    kind: string | null,
    key: string,
    w: number,
  ) {
    this.aStart.push(as);
    this.aEnd.push(ae);
    this.bStart.push(bs);
    this.bEnd.push(be);
    this.aStrand.push(sa < 0 ? -1 : 1);
    this.bStrand.push(sb < 0 ? -1 : 1);
    let k = -1;
    if (kind !== null) {
      k = this.kindAt.get(kind) ?? -1;
      if (k < 0) this.kindAt.set(kind, (k = this.kinds.push(kind) - 1));
    }
    this.kind.push(k);
    this.keys.push(key);
    this.weight.push(w);
  }

  data(extra?: Record<string, readonly (string | number | null)[]>): PairsData {
    return {
      n: this.keys.length,
      aStart: Int32Array.from(this.aStart),
      aEnd: Int32Array.from(this.aEnd),
      bStart: Int32Array.from(this.bStart),
      bEnd: Int32Array.from(this.bEnd),
      kind: Int16Array.from(this.kind),
      kinds: this.kinds,
      keys: this.keys,
      aStrand: Int8Array.from(this.aStrand),
      bStrand: Int8Array.from(this.bStrand),
      ...(extra ? { extra } : {}),
    };
  }
}

/**
 * Items tiles as pairs: every row of the tiles read (a tile is read uncut, wider than the view: a pan inside the tiles
 * builds nothing, and an arc across the view is drawn), each once across tiles by the anchor rule, on the layer's two
 * systems (a pair found from its `b` end turned).
 */
function fromItems(entries: readonly Loaded[], region: PairsRegion) {
  const sys = region.src.system;
  const ends = endsOf(region);
  const g = new Gather();
  for (let k = 0; k < entries.length; k++) {
    const t = entries[k].tables.pairs;
    if (t?.form !== 'list') continue;
    for (let i = 0; i < t.aStart.length; i++) {
      const o = orient(names(t.aSystem, i), names(t.bSystem, i), region);
      if (!o) continue;
      const as = ends & 1 && (names(t.aSystem, i) ?? sys) === sys ? t.aStart[i] : Number.NaN;
      const bs = ends & 2 && (names(t.bSystem, i) ?? sys) === sys ? t.bStart[i] : Number.NaN;
      if (!owned(entries, k, as, t.aEnd[i], bs, t.bEnd[i])) continue;
      const key = itemKey(t, i);
      const sa = t.aStrand?.[i] ?? 1;
      const sb = t.bStrand?.[i] ?? 1;
      const type = t.type?.get(i) ?? null;
      if (o > 0) g.push(t.aStart[i], t.aEnd[i], t.bStart[i], t.bEnd[i], sa, sb, type, key, 1);
      else g.push(t.bStart[i], t.bEnd[i], t.aStart[i], t.aEnd[i], sb, sa, type, key, 1);
    }
  }
  return { data: g.data(), weight: undefined };
}

/**
 * Bins tiles as pairs: each bin pair once (the anchor rule), from its bins' starts to their ends, its relative strand on
 * its `b` end (so a reversed bin is twisted, GA's rule), its count its weight (and a `count` field a theme reads).
 */
function fromBins(entries: readonly Loaded[], level: number, region: PairsRegion) {
  const sys = region.src.system;
  const ends = endsOf(region);
  const g = new Gather();
  for (let k = 0; k < entries.length; k++) {
    const t = entries[k].tables.pairs;
    if (t?.form !== 'bins' || 'aStart' in t) continue;
    for (let i = 0; i < t.count.length; i++) {
      const a = names(t.aSystem, i);
      const b = names(t.bSystem, i);
      const o = orient(a, b, region);
      if (!o) continue;
      const as = ends & 1 && (a ?? sys) === sys ? t.aBin[i] * level : Number.NaN;
      const bs = ends & 2 && (b ?? sys) === sys ? t.bBin[i] * level : Number.NaN;
      if (!owned(entries, k, as, as + level, bs, bs + level)) continue;
      const st = t.strand?.[i] ?? 0;
      const cat = t.category?.get(i) ?? null;
      const key = `bin:${a ?? ''}:${t.aBin[i]}|${b ?? ''}:${t.bBin[i]}|${st}|${cat ?? ''}`;
      const [x, y] = o > 0 ? [t.aBin[i], t.bBin[i]] : [t.bBin[i], t.aBin[i]];
      g.push(x * level, (x + 1) * level, y * level, (y + 1) * level, 1, st < 0 ? -1 : 1, cat, key, t.count[i]);
    }
  }
  return { data: g.data({ count: g.weight }), weight: Float32Array.from(g.weight) };
}

/** The layer's prep over the frame's pairs: every pair shown, its index by key, its elements made on hit. */
export function withData(
  prep: PairsPrep,
  data: PairsData,
  weight: Float32Array | undefined,
  level: number | undefined,
): PairsPrep {
  const index = new Map<string | number, number>();
  for (let k = 0; k < data.n; k++) if (!index.has(data.keys[k])) index.set(data.keys[k], k);
  return {
    ...prep,
    data,
    datum: new PairDatum(data),
    shown: new Uint8Array(data.n).fill(1),
    filtered: 0,
    index,
    elements: [],
    extents: [],
    units: new WeakMap(),
    ...(weight ? { weight } : {}),
    ...(level !== undefined ? { level } : {}),
  };
}
