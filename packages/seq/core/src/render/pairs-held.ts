/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import {
  binPairs,
  coarsenPairBins,
  DefaultSummarySettings,
  type PairsBins,
  type PairsList,
} from '@molstar/seqx-core/data';
import type { AxisTransform } from './axis.js';
import { boundedLevel } from './layer-data.js';
import type { LayoutContext, StoreReader } from './layer-type.js';
import { endsIn, type PairsPrep } from './pairs-lay.js';
import type { PairsReadOptions } from './pairs-reader.js';
import { Gather, withData } from './pairs-region.js';

// Held pairs by the pairs level (_spec/README.md "Render", The pairs layer: held pairs past the gate): a held table
// (inline pairs, a small generated genome's links, inline links over an assembly on its blocks) whose pairs the arcs
// reach (`inView`) pass `pairs.maxArcs` is binned as the store's summaries bin a level-1 source's pairs: levels
// `summaries.factor` apart up to the longest system, each from the one below, only those that halve the rows below
// kept; the level read the coarsest kept no wider than `pairs.binPx` device pixels, each bin pair one arc stroked by
// its count. So a held table and the same pairs from a file are one picture. Over an assembly each end is on its row's
// unit (a sequence), binned on it. The pyramid is a prep's (its filter's shown pairs, rebuilt with the prep), each
// level built when first read and kept, each level's prep once.

interface HeldPyramid {
  readonly list: PairsList;
  /** The longest system: no level at or past it. */
  readonly length: number;
  /** Rows a unit of the ends: the bins' systems are units' indices. */
  readonly units: boolean;
  /** The numeric field summed into the bins (`PairsReadOptions.sum`), when one is. */
  readonly sum?: string;
  /** The kept levels built so far, ascending, and their bins. */
  readonly bins: Map<number, PairsBins>;
  readonly preps: Map<number, PairsPrep>;
  /** The last level built (kept or not) and its bins; the rows below the next; the next level to build. */
  below?: PairsBins;
  belowLevel: number;
  rows: number;
  next: number;
}

/** A prep's pyramids, by the field summed (`''`: counts alone). */
const Pyramids = new WeakMap<PairsPrep, Map<string, HeldPyramid>>();

/**
 * Pairs of `prep` the arcs layout reaches at the axis's window (a pair once), its filter's shown ones only: an end in
 * it, its ends on either side of it (an arc across the window), or both past one edge within `reach` of the window (in
 * the arcs' scale before they show). So the gate counts what the layout draws, never fewer: a pair across the view
 * never turns the items to device-column bins the gate did not ask for.
 */
function inView(prep: PairsPrep, ctx: LayoutContext, axis: AxisTransform): number {
  const ends = endsIn(prep, ctx.section, axis);
  const u = ends.u;
  const r = prep.reach * (axis.u1 - axis.u0);
  const lo = axis.u0 - r;
  const hi = axis.u1 + r;
  let n = 0;
  let last = -1;
  for (let e = 0; e < ends.n; e++) {
    const k = ends.pair[e];
    if (k === last || !prep.shown[k]) continue;
    const a = u[2 * e];
    const b = u[2 * e + 1];
    const inA = a >= lo && a <= hi;
    const inB = b >= lo && b <= hi;
    if (inA || inB || (a < axis.u0 && b > axis.u1) || (b < axis.u0 && a > axis.u1)) {
      n++;
      last = k;
    }
  }
  return n;
}

/**
 * Held pairs over an assembly's own parts (`pairs-assembly.ts`), loaded with the genome chunk as the region reader is:
 * each end's unit as the list's systems and the longest unit, the blocks' scale, the bins' units back. Why lazy: a
 * viewer that never opens a genome pays nothing.
 */
export const OverAssembly: {
  list?(prep: PairsPrep, rows: readonly number[]): Pick<PairsList, 'aSystem' | 'bSystem'> & { length: number };
  scale?(ctx: LayoutContext): { ppp: number; span: number } | undefined;
  units?(bins: PairsBins): { aUnit: Int32Array; bUnit: Int32Array };
} = {};

/** The shown pairs as a list (over an assembly, each end's unit its system), the pyramid's levels yet unbuilt. */
function pyramidOf(prep: PairsPrep, sum: string | undefined): HeldPyramid {
  let byField = Pyramids.get(prep);
  if (!byField) Pyramids.set(prep, (byField = new Map()));
  const d = prep.data!;
  // a field summed only when the data holds it as numbers (else the counts alone, said by the reader's owner)
  const summed = sum !== undefined && d.extra?.[sum]?.some((v) => typeof v === 'number') ? sum : undefined;
  let p = byField.get(summed ?? '');
  if (p) return p;
  const rows: number[] = [];
  for (let k = 0; k < d.n; k++) if (prep.shown[k]) rows.push(k);
  const take = <A extends Int32Array | Int8Array>(c: A): A =>
    (c.constructor as unknown as { from(a: number[]): A }).from(rows.map((k) => c[k]));
  const units = prep.endUnits && OverAssembly.list!(prep, rows);
  const lengthOf = (s: string) => prep.model.systems.get(s)?.length ?? 0;
  const list: PairsList = {
    form: 'list',
    ...units,
    aStart: take(d.aStart),
    aEnd: take(d.aEnd),
    bStart: take(d.bStart),
    bEnd: take(d.bEnd),
    aStrand: take(d.aStrand),
    bStrand: take(d.bStrand),
    type: categoryColumn(rows.map((k) => (d.kind[k] < 0 ? null : d.kinds[d.kind[k]]))),
    ...(summed
      ? {
          fields: [{ name: summed, column: Float64Array.from(rows, (k) => Number(d.extra![summed][k]) || 0) as never }],
        }
      : {}),
  };
  p = {
    list,
    ...(summed ? { sum: summed } : {}),
    length: units?.length ?? Math.max(lengthOf(prep.aSystem), lengthOf(prep.bSystem)),
    units: !!units,
    bins: new Map(),
    preps: new Map(),
    belowLevel: 1,
    rows: rows.length,
    next: DefaultSummarySettings.factor,
  };
  byField.set(summed ?? '', p);
  return p;
}

/**
 * Builds the pyramid's levels until one kept is at least `target` positions (the store's summaries' rule: each from the
 * one below, kept when it halves the rows below), each once: a view that never zooms out builds none.
 */
function buildTo(p: HeldPyramid, target: number): void {
  let top = 1;
  for (const l of p.bins.keys()) top = l;
  while (p.next < p.length && (p.next <= target || top < target)) {
    const level = p.next;
    const b = p.below
      ? coarsenPairBins(p.below, level / p.belowLevel)
      : binPairs(p.list, level, { category: true, systems: p.units, ...(p.sum ? { field: p.sum } : {}) });
    if (b.count.length * 2 <= p.rows) p.bins.set((top = level), b);
    p.below = b;
    p.belowLevel = level;
    p.rows = b.count.length;
    p.next = level * DefaultSummarySettings.factor;
  }
}

/**
 * The binned prep of held pairs past the gate, or undefined when they draw as items (under it, at the items' level,
 * a complex's links, a grid's dots: not binned). `options`: the reader's bin width and gate (the pairs config's when
 * absent).
 */
export function heldBins(
  prep: PairsPrep,
  ctx: LayoutContext,
  data: StoreReader,
  options?: PairsReadOptions,
): PairsPrep | undefined {
  const d = prep.data;
  if (!d || prep.region || prep.weight || ctx.axisY) return undefined;
  // a complex's links stay items (their units are a structure's, binned by device column past the bound); over an
  // assembly, once its part has loaded
  if (prep.endUnits && (!OverAssembly.list || prep.model.systems.get(prep.endUnits[0])?.kind === 'unit'))
    return undefined;
  const cfg = ctx.config.pairs;
  const gate = options?.gate ?? cfg.maxArcs;
  if (d.n <= gate || inView(prep, ctx, ctx.axis) <= gate) return undefined;
  const p = pyramidOf(prep, options?.sum);
  const { ppp, span } = OverAssembly.scale?.(ctx) ?? { ppp: 1 / ctx.axis.devicePx, span: ctx.axis.u1 - ctx.axis.u0 };
  const need = (options?.binPx ?? cfg.binPx) * ppp;
  // as a source's bins are read (`boundedLevel`): no wider than `binPx` pixels, coarser while the view's tiles at it
  // would pass the store's bound, so a held table and the same pairs from a file read one level
  const s = data.settings;
  buildTo(p, Math.max(need, s ? span / (s.tileBins * Math.max(1, s.maxTilesPerWant - 1)) : 0));
  const level = boundedLevel([1, ...p.bins.keys()], need, span, data);
  if (level <= 1) return undefined;
  let out = p.preps.get(level);
  if (out) return out;
  const bins = p.bins.get(level)!;
  const g = new Gather();
  for (let i = 0; i < bins.count.length; i++) {
    const x = bins.aBin[i];
    const y = bins.bBin[i];
    const st = bins.strand?.[i] ?? 0;
    const cat = bins.category?.get(i) ?? null;
    const a = bins.aSystem?.get(i) ?? '';
    const b = bins.bSystem?.get(i) ?? '';
    const key = `bin:${a}:${x}|${b}:${y}|${st}|${cat ?? ''}`;
    g.push(x * level, (x + 1) * level, y * level, (y + 1) * level, 1, st < 0 ? -1 : 1, cat, key, bins.count[i]);
  }
  // a summed field rides along as a column of its name (the bins' sums, in the bins' order)
  const pairs = g.data(p.sum && bins.sum ? { count: g.weight, [p.sum]: Array.from(bins.sum) } : { count: g.weight });
  const binned = withData(
    prep,
    p.units ? { ...pairs, ...OverAssembly.units!(bins) } : pairs,
    Float32Array.from(g.weight),
    level,
  );
  out = { ...binned, filtered: prep.filtered };
  p.preps.set(level, out);
  return out;
}
