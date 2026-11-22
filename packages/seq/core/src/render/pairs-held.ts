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
import { Gather, withData } from './pairs-region.js';

// Held pairs by the pairs level (_spec/README.md "Render", The pairs layer: held pairs past the gate): a held table
// (inline pairs, a small generated genome's links) whose pairs with an end in view pass `pairs.maxArcs` is binned as
// the store's summaries bin a level-1 source's pairs: levels `summaries.factor` apart up to the system's length, each
// from the one below, only those that halve the rows below kept; the level read the coarsest kept no wider than
// `pairs.binPx` device pixels (coarser while a view's tiles would pass the store's bound), each bin pair one arc
// stroked by its count. So a held table and the same pairs from a file are one picture. The pyramid is built once a
// prep (from the pairs its filter shows), each level's prep once.

interface HeldPyramid {
  /** 1 and the kept levels, ascending. */
  readonly levels: number[];
  readonly bins: Map<number, PairsBins>;
  readonly preps: Map<number, PairsPrep>;
}

const Pyramids = new WeakMap<PairsPrep, HeldPyramid>();

/** Pairs of `prep` with an end in the axis's window (a pair once), its filter's shown ones only. */
function inView(prep: PairsPrep, ctx: LayoutContext, axis: AxisTransform): number {
  const ends = endsIn(prep, ctx.section, axis);
  const u = ends.u;
  let n = 0;
  let last = -1;
  for (let e = 0; e < ends.n; e++) {
    const k = ends.pair[e];
    if (k === last || !prep.shown[k]) continue;
    const a = u[2 * e];
    const b = u[2 * e + 1];
    if ((a >= axis.u0 && a <= axis.u1) || (b >= axis.u0 && b <= axis.u1)) {
      n++;
      last = k;
    }
  }
  return n;
}

/** The pyramid of the shown pairs, as the store's summaries build a source's (`pairsPyramid`'s levels rule). */
function pyramidOf(prep: PairsPrep): HeldPyramid {
  let p = Pyramids.get(prep);
  if (p) return p;
  const d = prep.data!;
  const rows: number[] = [];
  for (let k = 0; k < d.n; k++) if (prep.shown[k]) rows.push(k);
  const take = <A extends Int32Array | Int8Array>(c: A): A =>
    (c.constructor as unknown as { from(a: number[]): A }).from(rows.map((k) => c[k]));
  const list: PairsList = {
    form: 'list',
    aStart: take(d.aStart),
    aEnd: take(d.aEnd),
    bStart: take(d.bStart),
    bEnd: take(d.bEnd),
    aStrand: take(d.aStrand),
    bStrand: take(d.bStrand),
    type: categoryColumn(rows.map((k) => (d.kind[k] < 0 ? null : d.kinds[d.kind[k]]))),
  };
  const length = Math.max(
    prep.model.systems.get(prep.aSystem)?.length ?? 0,
    prep.model.systems.get(prep.bSystem)?.length ?? 0,
  );
  const f = DefaultSummarySettings.factor;
  const levels = [1];
  const bins = new Map<number, PairsBins>();
  let below: PairsBins | undefined;
  let belowLevel = 1;
  let belowRows = rows.length;
  for (let level = f; level < length; level *= f) {
    const b: PairsBins = below ? coarsenPairBins(below, level / belowLevel) : binPairs(list, level, { category: true });
    if (b.count.length * 2 <= belowRows) {
      levels.push(level);
      bins.set(level, b);
    }
    below = b;
    belowLevel = level;
    belowRows = b.count.length;
  }
  p = { levels, bins, preps: new Map() };
  Pyramids.set(prep, p);
  return p;
}

/**
 * The binned prep of held pairs past the gate, or undefined when they draw as items (under it, at the items' level,
 * a complex's links, a grid's dots, blocks: not binned).
 */
export function heldBins(prep: PairsPrep, ctx: LayoutContext, data: StoreReader): PairsPrep | undefined {
  const d = prep.data;
  if (!d || prep.region || prep.endUnits || prep.weight || ctx.axisY || ctx.section.blocks) return undefined;
  const cfg = ctx.config.pairs;
  if (d.n <= cfg.maxArcs || inView(prep, ctx, ctx.axis) <= cfg.maxArcs) return undefined;
  const p = pyramidOf(prep);
  const level = boundedLevel(p.levels, cfg.binPx / ctx.axis.devicePx, ctx.axis.u1 - ctx.axis.u0, data);
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
    const key = `bin::${x}|:${y}|${st}|${cat ?? ''}`;
    g.push(x * level, (x + 1) * level, y * level, (y + 1) * level, 1, st < 0 ? -1 : 1, cat, key, bins.count[i]);
  }
  const binned = withData(prep, g.data({ count: g.weight }), Float32Array.from(g.weight), level);
  out = { ...binned, filtered: prep.filtered };
  p.preps.set(level, out);
  return out;
}
