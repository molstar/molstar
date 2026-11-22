/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CategoryColumn } from '@molstar/core/repr/tree';
import type { PairsList } from '@molstar/seqx-core/data';
import type { PairsData } from '@molstar/seqx-core/resolve';
import { dataRegions, itemsGate, type LayerSource, pairsLevels, readTiles, type TileRead } from './layer-data.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import { PairDatum, type PairsLay, type PairsPrep } from './pairs-lay.js';

// Pairs read by region (_spec/README.md "Render", The pairs layer: by region): the layer's one reader over a source's
// pairs, as held pairs are read. Two reads: bins at a pixel's level (`pairs.binPx`), always; items as items tiles at the
// finest level whose tiles fit, asked once the bin pairs with an end in view sum to at most `pairs.maxArcs` and kept
// while the items in view stay within `pairs.itemsHysteresis` × it. Bins draw until the items land (EC-169: pending is
// never empty), each bin pair one arc between its bins' middles, stroked by count. What is read becomes the frame's
// `PairsData` (ends on the layer's two systems, a pair found from its other end turned round; every row of the tiles
// read, each once by the anchor rule), kept while the entries read are the same, so a warm frame builds nothing.

/** What a region read of pairs keeps on its lay between frames. */
export interface RegionState {
  /** What the data was built from: the mode (0 bins, 1 items), the level, then each entry's tile and version. */
  readonly read: number[];
  prep: PairsPrep | undefined;
  /** Whether items are shown (the gate's state, for its hysteresis). */
  items: boolean;
  /** Bin pairs with an end in view (summed counts), as last read; items in view when shown. */
  binPairs: number;
  itemCount: number;
}

/** A lay's region state, made on first read. */
const States = new WeakMap<PairsLay, RegionState>();

export function regionState(lay: PairsLay): RegionState {
  let s = States.get(lay);
  if (!s) States.set(lay, (s = { read: [], prep: undefined, items: false, binPairs: 0, itemCount: 0 }));
  return s;
}

/** Where a region pairs layer reads: its source, the names of its two systems as the source knows them. */
export interface PairsRegion {
  readonly src: LayerSource;
  readonly aName: string;
  readonly bName: string;
}

/** Whether `[s, e)` meets any of `regions` (`[start, end, …]`). */
function inView(regions: readonly number[], s: number, e: number): boolean {
  for (let i = 0; i + 1 < regions.length; i += 2)
    if (s < regions[i + 1] && Math.max(e, s + 1) > regions[i]) return true;
  return false;
}

/**
 * The frame's pairs of a region layer: reads bins, gates items, and gives the prep the arcs lay out (`undefined` while
 * nothing is read). The read's state is the lay's.
 */
export function regionPairs(
  prep: PairsPrep,
  region: PairsRegion,
  ctx: LayoutContext,
  data: StoreReader,
  lay: PairsLay,
): PairsPrep | undefined {
  const state = regionState(lay);
  const { src } = region;
  const cov = data.coverage(src.source);
  const levels = cov?.levels ?? [1];
  const regions = dataRegions(src, ctx.axis);
  const span = (ctx.axis.u1 - ctx.axis.u0) * src.ratio;
  const cfg = ctx.config.pairs;
  const { bins, items } = pairsLevels(levels, src.ratio / ctx.axis.devicePx, span, data, cfg.binPx);
  const binRead = readTiles(ctx, data, src, regions, bins);
  // at the finest level the bins are the items themselves
  const binned = binRead.entries.some((e) => e.tables.pairs?.form === 'bins');
  state.binPairs = binned
    ? countBins(binRead.entries, bins, region, regions)
    : countItems(binRead.entries, region, regions);
  let mode = binned ? 0 : 1;
  let read = binRead.entries;
  if (binned) {
    const ask = itemsGate(
      state.items,
      state.binPairs,
      state.items ? state.itemCount : undefined,
      cfg.maxArcs,
      cfg.itemsHysteresis,
    );
    state.items = false;
    if (ask) {
      const itemRead = readTiles(ctx, data, src, regions, items, { items: true });
      if (!itemRead.pending.length && !itemRead.failed.length) {
        state.itemCount = countItems(itemRead.entries, region, regions);
        state.items = state.itemCount <= cfg.maxArcs * Math.max(1, cfg.itemsHysteresis);
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

type Loaded = TileRead['entries'][number];

/** Whether the data was built from this read (its mode, level, tiles and versions): nothing made to tell. */
function sameRead(r: readonly number[], mode: number, level: number, entries: readonly Loaded[]): boolean {
  if (r.length !== 2 + 2 * entries.length || r[0] !== mode || r[1] !== level) return false;
  for (let k = 0; k < entries.length; k++)
    if (r[2 + 2 * k] !== entries[k].region.tile || r[3 + 2 * k] !== entries[k].version) return false;
  return true;
}

function keepRead(r: number[], mode: number, level: number, entries: readonly Loaded[]) {
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
function names(sys: CategoryColumn | undefined, i: number): string | undefined {
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

/** Bin pairs with an end in view, their counts summed, each once across tiles (the anchor rule: no key made). */
function countBins(entries: readonly Loaded[], level: number, region: PairsRegion, regions: readonly number[]): number {
  const sys = region.src.system;
  let n = 0;
  for (let k = 0; k < entries.length; k++) {
    const t = entries[k].tables.pairs;
    if (t?.form !== 'bins' || 'aStart' in t) continue;
    for (let i = 0; i < t.count.length; i++) {
      const aOn = (names(t.aSystem, i) ?? sys) === sys;
      const bOn = (names(t.bSystem, i) ?? sys) === sys;
      const as = aOn ? t.aBin[i] * level : Number.NaN;
      const bs = bOn ? t.bBin[i] * level : Number.NaN;
      if (!(aOn && inView(regions, as, as + level)) && !(bOn && inView(regions, bs, bs + level))) continue;
      if (owned(entries, k, as, as + level, bs, bs + level)) n += t.count[i];
    }
  }
  return n;
}

/** A pair's identity across tiles: its row, else its id (not `.` or empty), with both ends; else its ends alone. */
function itemKey(t: PairsList, i: number): string {
  const named = t.row ? `r${t.row[i]}` : t.id ? t.id.get(i) : null;
  const ends = `${names(t.aSystem, i)}:${t.aStart[i]}-${t.aEnd[i]}|${names(t.bSystem, i)}:${t.bStart[i]}-${t.bEnd[i]}`;
  return named && named !== '.' ? `${named}|${ends}` : ends;
}

/** Pairs with an end in view, each once across tiles (the anchor rule: no key made). */
function countItems(entries: readonly Loaded[], region: PairsRegion, regions: readonly number[]): number {
  const sys = region.src.system;
  let n = 0;
  for (let k = 0; k < entries.length; k++) {
    const t = entries[k].tables.pairs;
    if (t?.form !== 'list') continue;
    for (let i = 0; i < t.aStart.length; i++) {
      const aOn = (names(t.aSystem, i) ?? sys) === sys;
      const bOn = (names(t.bSystem, i) ?? sys) === sys;
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
class Gather {
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
  const g = new Gather();
  for (let k = 0; k < entries.length; k++) {
    const t = entries[k].tables.pairs;
    if (t?.form !== 'list') continue;
    for (let i = 0; i < t.aStart.length; i++) {
      const o = orient(names(t.aSystem, i), names(t.bSystem, i), region);
      if (!o) continue;
      const as = (names(t.aSystem, i) ?? sys) === sys ? t.aStart[i] : Number.NaN;
      const bs = (names(t.bSystem, i) ?? sys) === sys ? t.bStart[i] : Number.NaN;
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
  const g = new Gather();
  for (let k = 0; k < entries.length; k++) {
    const t = entries[k].tables.pairs;
    if (t?.form !== 'bins' || 'aStart' in t) continue;
    for (let i = 0; i < t.count.length; i++) {
      const a = names(t.aSystem, i);
      const b = names(t.bSystem, i);
      const o = orient(a, b, region);
      if (!o) continue;
      const as = (a ?? sys) === sys ? t.aBin[i] * level : Number.NaN;
      const bs = (b ?? sys) === sys ? t.bBin[i] * level : Number.NaN;
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
function withData(
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
