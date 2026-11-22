/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Want } from '@molstar/seqx-core/data';
import type { AxisTransform } from './axis.js';
import { dataRegions, itemsGate, type LayerSource, pairsLevels, readTiles } from './layer-data.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import type { PairsLay, PairsPrep } from './pairs-lay.js';
import {
  type AssemblyRegion,
  Gather,
  type Loaded,
  names,
  type PairsRegion,
  regionState,
  withData,
} from './pairs-region.js';

// Pairs over an assembly read by region (_spec/README.md "Render", The pairs layer: over an assembly): each block in
// view (a strip: blocks under a device pixel have none, the joint small-blocks rule) reads its own sequence's tiles,
// bins at one level for every block (a bin pair across two blocks is one bin pair from either end), items once the bin
// pairs with an end in view sum to at most `pairs.maxArcs`, as on one system. A row is found from either end's block:
// it is taken once (by its id and ends, a bin pair by its bins), each end on the unit its listed name is (`PairsData.aUnit`, `bUnit` over the
// assembly's systems), so the pairs layer places it on its blocks (a chord between two chromosomes on a ring). The
// layer's wants over every block carry its key, so the store bounds them together (`maxTilesPerLayer`); the level is
// coarser while the blocks' wants together do not fit (`StoreReader.layerLoad`), else the band's "zoom in".

/** A block read this frame: its source, its positions in view, its unit. */
interface BlockRead {
  src: LayerSource;
  regions: number[];
  unit: number;
}

/** What an assembly read keeps between frames, beside the region state. */
interface AssemblyState {
  /** The blocks in view by unit, this frame (`-1`: none); the reads, reused. */
  blocks: BlockRead[];
  inView: number[][];
  /** The level last fitted to the store's bounds and the view's tiles it was fitted for. */
  fitSig: number[];
  fitted: number;
}

const States = new WeakMap<PairsLay, AssemblyState>();

function stateOf(lay: PairsLay): AssemblyState {
  let s = States.get(lay);
  if (!s) States.set(lay, (s = { blocks: [], inView: [], fitSig: [], fitted: 0 }));
  return s;
}

/** Whether `[s, e)` meets any of `regions` (`[start, end, …]`). */
function meets(regions: readonly number[] | undefined, s: number, e: number): boolean {
  if (!regions) return false;
  for (let i = 0; i + 1 < regions.length; i += 2)
    if (s < regions[i + 1] && Math.max(e, s + 1) > regions[i]) return true;
  return false;
}

/** The blocks in view (by the strips laid out), each once with its positions in view merged. */
function blocksInView(over: AssemblyRegion, ctx: LayoutContext, st: AssemblyState): number {
  const strips = ctx.place!.strips;
  st.inView.length = over.units.length;
  for (let u = 0; u < over.units.length; u++) st.inView[u] = st.inView[u] ?? [];
  for (let u = 0; u < over.units.length; u++) st.inView[u].length = 0;
  let n = 0;
  for (let i = 0; i < strips.n; i++) {
    const axis = strips.axes[i] as AxisTransform;
    const src = over.sources.get(axis.system);
    if (!src) continue;
    const unit = over.unitOf.get(src.system)!;
    const regions = dataRegions(src, axis);
    if (!regions.length) continue;
    const seen = st.inView[unit];
    if (!seen.length) {
      const b = st.blocks[n] ?? (st.blocks[n] = { src, regions: [], unit });
      b.src = src;
      b.unit = unit;
      b.regions = seen;
      n++;
    }
    seen.push(...regions);
  }
  return n;
}

/** The finest level whose blocks' wants together fit the store's bounds (coarser first, then "zoom in"). */
function fitBlocks(
  ctx: LayoutContext,
  data: StoreReader,
  st: AssemblyState,
  n: number,
  levels: readonly number[],
  level: number,
) {
  if (!data.layerLoad) return level;
  const bins = data.settings?.tileBins ?? 1024;
  const sig = st.fitSig;
  let k = 0;
  const put = (v: number) => {
    if (sig[k] !== v) sig.length = k;
    sig[k++] = v;
  };
  const was = sig.length;
  put(level);
  for (let b = 0; b < n; b++) {
    put(st.blocks[b].unit);
    for (const r of st.blocks[b].regions) put(Math.floor(r / (bins * level)));
  }
  if (k === was && sig.length === k && st.fitted) return st.fitted;
  sig.length = k;
  let fitted = level;
  for (let i = Math.max(0, levels.indexOf(level)); i < levels.length; i++) {
    const l = levels[i];
    if (l < level) continue;
    const wants: Want[] = [];
    for (let b = 0; b < n; b++) {
      const { src, regions } = st.blocks[b];
      for (let r = 0; r + 1 < regions.length; r += 2)
        wants.push({
          source: src.source,
          system: src.system,
          start: regions[r],
          end: regions[r + 1],
          level: l,
          layer: ctx.key,
        });
    }
    if (data.layerLoad(wants).fits) {
      fitted = l;
      break;
    }
  }
  st.fitted = fitted;
  return fitted;
}

/** Tiles of every block in view at `level`, each entry with its block's unit (`units`, parallel). */
function readBlocks(
  ctx: LayoutContext,
  data: StoreReader,
  st: AssemblyState,
  n: number,
  level: number,
  items: boolean,
  units: number[],
): { entries: Loaded[]; pending: boolean } {
  const entries: Loaded[] = [];
  let pending = false;
  units.length = 0;
  for (let b = 0; b < n; b++) {
    const { src, regions, unit } = st.blocks[b];
    const r = readTiles(ctx, data, src, regions, level, items ? { items: true } : undefined);
    if (r.pending.length || r.failed.length) pending = true;
    for (const e of r.entries) {
      entries.push(e);
      units.push(unit);
    }
  }
  return { entries, pending };
}

/** The unit of a row's end: its listed name's, else (no name) the block's it was read on. */
function unitOfEnd(over: AssemblyRegion, name: string | undefined, block: number): number {
  return name === undefined ? block : (over.unitOf.get(name) ?? -1);
}

/** What a read gives: rows once each (by key), their ends' units; counts of those with an end in view. */
function gather(
  over: AssemblyRegion,
  st: AssemblyState,
  entries: readonly Loaded[],
  units: readonly number[],
  level: number,
  bins: boolean,
  build: boolean,
): { count: number; g?: Gather; ua: number[]; ub: number[] } {
  const seen = new Set<string>();
  const g = build ? new Gather() : undefined;
  const ua: number[] = [];
  const ub: number[] = [];
  let count = 0;
  for (let k = 0; k < entries.length; k++) {
    const t = entries[k].tables.pairs;
    if (!t) continue;
    if (bins ? t.form !== 'bins' || 'aStart' in t : t.form !== 'list') continue;
    const rows = t.form === 'bins' ? t.count.length : t.aStart.length;
    for (let i = 0; i < rows; i++) {
      const a = names(t.aSystem, i);
      const b = names(t.bSystem, i);
      const au = unitOfEnd(over, a, units[k]);
      const bu = unitOfEnd(over, b, units[k]);
      if (au < 0 || bu < 0) continue;
      let as: number;
      let ae: number;
      let bs: number;
      let be: number;
      let key: string;
      if (t.form === 'bins' && !('aStart' in t)) {
        as = t.aBin[i] * level;
        ae = as + level;
        bs = t.bBin[i] * level;
        be = bs + level;
        key = `bin:${a ?? over.units[au]}:${t.aBin[i]}|${b ?? over.units[bu]}:${t.bBin[i]}|${t.strand?.[i] ?? 0}|${t.category?.get(i) ?? ''}`;
      } else if (t.form === 'list') {
        as = t.aStart[i];
        ae = t.aEnd[i];
        bs = t.bStart[i];
        be = t.bEnd[i];
        // its id and ends, never its row: a source numbers the rows of each sequence's tile apart, and a link is
        // listed on both its sequences
        const id = t.id?.get(i);
        key = `${id && id !== '.' ? id : ''}|${au}:${as}-${ae}|${bu}:${bs}-${be}`;
      } else continue;
      if (!meets(st.inView[au], as, ae) && !meets(st.inView[bu], bs, be)) continue;
      if (seen.has(key)) continue;
      seen.add(key);
      if (t.form === 'bins' && !('aStart' in t)) {
        count += t.count[i];
        const strand = t.strand?.[i] ?? 0;
        g?.push(as, ae, bs, be, 1, strand < 0 ? -1 : 1, t.category?.get(i) ?? null, key, t.count[i]);
      } else if (t.form === 'list') {
        count++;
        g?.push(as, ae, bs, be, t.aStrand?.[i] ?? 1, t.bStrand?.[i] ?? 1, t.type?.get(i) ?? null, key, 1);
      }
      if (g) {
        ua.push(au);
        ub.push(bu);
      }
    }
  }
  return { count, ...(g ? { g } : {}), ua, ub };
}

/** Whether the data was built from this read (its mode, level, and each entry's block, tile and version). */
function sameRead(
  r: readonly number[],
  mode: number,
  level: number,
  entries: readonly Loaded[],
  units: readonly number[],
) {
  if (r.length !== 2 + 3 * entries.length || r[0] !== mode || r[1] !== level) return false;
  for (let k = 0; k < entries.length; k++)
    if (r[2 + 3 * k] !== units[k] || r[3 + 3 * k] !== entries[k].region.tile || r[4 + 3 * k] !== entries[k].version)
      return false;
  return true;
}

const Units: number[] = [];
const ItemUnits: number[] = [];

/**
 * The frame's pairs of a layer over an assembly: each block in view read on its sequence, bins then items (the gate),
 * the rows found once with their ends' units. `undefined` while nothing is read.
 */
export function assemblyPairs(
  prep: PairsPrep,
  region: PairsRegion,
  over: AssemblyRegion,
  ctx: LayoutContext,
  data: StoreReader,
  lay: PairsLay,
): PairsPrep | undefined {
  if (!ctx.place) return undefined;
  const state = regionState(lay);
  const st = stateOf(lay);
  const n = blocksInView(over, ctx, st);
  const levels = data.coverage(region.src.source)?.levels ?? [1];
  // one level for every block: the finest any block in view needs, its span the widest
  let ppp = Number.POSITIVE_INFINITY;
  let span = 0;
  const strips = ctx.place.strips;
  for (let i = 0; i < strips.n; i++) {
    const axis = strips.axes[i];
    if (!over.sources.has(axis.system)) continue;
    ppp = Math.min(ppp, 1 / axis.devicePx);
    span = Math.max(span, axis.u1 - axis.u0);
  }
  if (!n || !Number.isFinite(ppp)) return state.prep;
  const cfg = ctx.config.pairs;
  const fit = pairsLevels(levels, ppp, span, data, cfg.binPx);
  const bins = fitBlocks(ctx, data, st, n, levels, fit.bins);
  const binRead = readBlocks(ctx, data, st, n, bins, false, Units);
  const binned = binRead.entries.some((e) => e.tables.pairs?.form === 'bins');
  state.binPairs = gather(over, st, binRead.entries, Units, bins, binned, false).count;
  let mode = binned ? 0 : 1;
  let read = binRead.entries;
  let units = Units;
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
      const itemRead = readBlocks(ctx, data, st, n, fit.items, true, ItemUnits);
      if (!itemRead.pending) {
        state.itemCount = gather(over, st, itemRead.entries, ItemUnits, fit.items, false, false).count;
        state.items = state.itemCount <= cfg.maxArcs * Math.max(1, cfg.itemsHysteresis);
        if (state.items) {
          mode = 1;
          read = itemRead.entries;
          units = ItemUnits;
        }
      } else state.items = true;
    }
  }
  const level = mode ? fit.items : bins;
  if (sameRead(state.read, mode, level, read, units) && state.prep) return state.prep;
  const built = gather(over, st, read, units, bins, mode === 0, true);
  const g = built.g!;
  const pairs = {
    ...g.data(mode ? undefined : { count: g.weight }),
    aUnit: Int32Array.from(built.ua),
    bUnit: Int32Array.from(built.ub),
  };
  state.read.length = 0;
  state.read.push(mode, level);
  for (let k = 0; k < read.length; k++) state.read.push(units[k], read[k].region.tile, read[k].version);
  state.prep = {
    ...withData(prep, pairs, mode ? undefined : Float32Array.from(g.weight), mode ? undefined : bins),
    endUnits: over.units,
  };
  return state.prep;
}
