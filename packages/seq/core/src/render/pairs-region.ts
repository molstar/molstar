/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CategoryColumn } from '@molstar/core/repr/tree';
import type { PairsBins, PairsList } from '@molstar/seqx-core/data';
import type { PairsData } from '@molstar/seqx-core/resolve';
import { dataRegions, itemsGate, type LayerSource, pairsLevels, readTiles } from './layer-data.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import { PairDatum, type PairsLay, type PairsPrep } from './pairs-lay.js';

// Pairs read by region (_spec/README.md "Render", The pairs layer: by region): the layer's one reader over a source's
// pairs, as held pairs are read. Two reads: bins at a pixel's level (`pairs.binPx`), always; items as items tiles at the
// finest level whose tiles fit, asked once the bin pairs with an end in view sum to at most `pairs.maxArcs` and kept
// while the items in view stay within `pairs.itemsHysteresis` × it. Bins draw until the items land (EC-169: pending is
// never empty), each bin pair one arc between its bins' middles, stroked by count. What is read becomes the frame's
// `PairsData` (ends on the layer's two systems, a pair found from its other end turned round), kept while the entries
// read are the same, so a warm frame builds nothing.

/** What a region read of pairs keeps on its lay between frames. */
export interface RegionState {
  /** The entries' key the data was built for (their versions and the mode). */
  key: string;
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
  if (!s) States.set(lay, (s = { key: '', prep: undefined, items: false, binPairs: 0, itemCount: 0 }));
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
  const binTables = binRead.entries.map((e) => e.tables.pairs).filter((t) => t !== undefined);
  // at the finest level the bins are the items themselves
  const lists = binTables.filter((t): t is PairsList => t.form === 'list');
  const binned = binTables.filter((t): t is PairsBins => t.form === 'bins');
  state.binPairs = countBins(binned, bins, region, regions) + countItems(lists, region, regions);
  let mode: 'bins' | 'items' = binned.length ? 'bins' : 'items';
  let itemTables: PairsList[] = lists;
  let itemKeys = '';
  if (binned.length) {
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
        itemTables = itemRead.entries.map((e) => e.tables.pairs).filter((t): t is PairsList => t?.form === 'list');
        state.itemCount = countItems(itemTables, region, regions);
        state.items = state.itemCount <= cfg.maxArcs * Math.max(1, cfg.itemsHysteresis);
        if (state.items) {
          mode = 'items';
          itemKeys = itemRead.entries.map((e) => `${e.region.tile}:${e.version}`).join(',');
        }
      } else state.items = true;
    }
  }
  const binKeys = binRead.entries.map((e) => `${e.region.tile}:${e.version}`).join(',');
  const key = `${mode}|${bins}|${mode === 'items' && itemKeys ? itemKeys : binKeys}`;
  if (state.key === key && state.prep) return state.prep;
  const built = mode === 'items' ? fromItems(itemTables, region, regions) : fromBins(binned, bins, region);
  state.key = key;
  state.prep = withData(prep, built.data, built.weight, mode === 'bins' ? bins : undefined);
  return state.prep;
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

function countBins(
  tables: readonly PairsBins[],
  level: number,
  region: PairsRegion,
  regions: readonly number[],
): number {
  const seen = new Set<string>();
  let n = 0;
  for (const t of tables)
    for (let i = 0; i < t.count.length; i++) {
      const a = names(t.aSystem, i);
      const b = names(t.bSystem, i);
      const sys = region.src.system;
      const aIn = (a ?? sys) === sys && inView(regions, t.aBin[i] * level, (t.aBin[i] + 1) * level);
      const bIn = (b ?? sys) === sys && inView(regions, t.bBin[i] * level, (t.bBin[i] + 1) * level);
      if (!aIn && !bIn) continue;
      const k = `${a}|${t.aBin[i]}|${b}|${t.bBin[i]}|${t.strand?.[i] ?? 0}|${t.category?.codes[i] ?? -1}`;
      if (seen.has(k)) continue;
      seen.add(k);
      n += t.count[i];
    }
  return n;
}

/** A pair's identity across tiles: its row, else its id (not `.` or empty), with both ends; else its ends alone. */
function itemKey(t: PairsList, i: number): string {
  const named = t.row ? `r${t.row[i]}` : t.id ? t.id.get(i) : null;
  const ends = `${names(t.aSystem, i)}:${t.aStart[i]}-${t.aEnd[i]}|${names(t.bSystem, i)}:${t.bStart[i]}-${t.bEnd[i]}`;
  return named && named !== '.' ? `${named}|${ends}` : ends;
}

function countItems(tables: readonly PairsList[], region: PairsRegion, regions: readonly number[]): number {
  const seen = new Set<string>();
  const sys = region.src.system;
  for (const t of tables)
    for (let i = 0; i < t.aStart.length; i++) {
      const aIn = (names(t.aSystem, i) ?? sys) === sys && inView(regions, t.aStart[i], t.aEnd[i]);
      const bIn = (names(t.bSystem, i) ?? sys) === sys && inView(regions, t.bStart[i], t.bEnd[i]);
      if (aIn || bIn) seen.add(itemKey(t, i));
    }
  return seen.size;
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
 * Items tiles as pairs: those with an end in view (a tile is read uncut, wider than the view), each once across tiles,
 * on the layer's two systems (a pair found from its `b` end turned).
 */
function fromItems(tables: readonly PairsList[], region: PairsRegion, regions: readonly number[]) {
  const sys = region.src.system;
  const g = new Gather();
  const seen = new Set<string>();
  for (const t of tables)
    for (let i = 0; i < t.aStart.length; i++) {
      const o = orient(names(t.aSystem, i), names(t.bSystem, i), region);
      if (!o) continue;
      const aIn = (names(t.aSystem, i) ?? sys) === sys && inView(regions, t.aStart[i], t.aEnd[i]);
      if (!aIn && !((names(t.bSystem, i) ?? sys) === sys && inView(regions, t.bStart[i], t.bEnd[i]))) continue;
      const key = itemKey(t, i);
      if (seen.has(key)) continue;
      seen.add(key);
      const sa = t.aStrand?.[i] ?? 1;
      const sb = t.bStrand?.[i] ?? 1;
      const type = t.type?.get(i) ?? null;
      if (o > 0) g.push(t.aStart[i], t.aEnd[i], t.bStart[i], t.bEnd[i], sa, sb, type, key, 1);
      else g.push(t.bStart[i], t.bEnd[i], t.aStart[i], t.aEnd[i], sb, sa, type, key, 1);
    }
  return { data: g.data(), weight: undefined };
}

/**
 * Bins tiles as pairs: each bin pair once, from its bins' starts to their ends, its relative strand on its `b` end (so
 * a reversed bin is twisted, GA's rule), its count its weight (and a `count` field a theme reads).
 */
function fromBins(tables: readonly PairsBins[], level: number, region: PairsRegion) {
  const g = new Gather();
  const seen = new Set<string>();
  for (const t of tables)
    for (let i = 0; i < t.count.length; i++) {
      const a = names(t.aSystem, i);
      const b = names(t.bSystem, i);
      const o = orient(a, b, region);
      if (!o) continue;
      const st = t.strand?.[i] ?? 0;
      const cat = t.category?.get(i) ?? null;
      const key = `bin:${a ?? ''}:${t.aBin[i]}|${b ?? ''}:${t.bBin[i]}|${st}|${cat ?? ''}`;
      if (seen.has(key)) continue;
      seen.add(key);
      const [x, y] = o > 0 ? [t.aBin[i], t.bBin[i]] : [t.bBin[i], t.aBin[i]];
      g.push(x * level, (x + 1) * level, y * level, (y + 1) * level, 1, st < 0 ? -1 : 1, cat, key, t.count[i]);
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
