/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Per-position cells (seqspec:viewer#extension-points `letters`, `color`; seqspec:sources#levels-and-the-zoom-gate),
// part of the renderer API: a cell per shown position where it is a device px or wider, else per device column the
// color theme's reduction. A `classes` theme reduces by the plurality of its block's (color, alpha) classes
// (`cellClasses`, built on the first aggregate draw per prepared state, theme and block); a `mean` theme by
// `palette(mean)` over a sum-and-count pyramid of its values. Raw positions below the pyramid's base, in a filtered
// block, in a column across a break, and in a block past 255 classes. A column takes the bins whose start lies in it,
// so a copied strip and a full redraw count the same bins. Columns are found in closed form. Allocation-free once warm:
// module scratch, grown only.

import { type Axis, AxisGeometry } from '../axis.js';
import { binClass, binValues, PYRAMID_BASE, type ValuePyramid, valuePyramid } from '../pyramids.js';
import { type CellClasses, type CellColoring, cellClasses } from '../reductions.js';
import { OPAQUE, type Theme, ThemeLocation } from '../themes.js';
import { NO_COLOR, type RectBatch } from './backend.js';
import type { DrawView } from './view.js';

/**
 * What cells read: a coloring by block index (`alpha` reads the location `color` filled for the same position), and
 * for a `mean` theme a block's values by `i − base` (NaN: none) with their pyramid when prepared; null values: classes.
 */
export interface CellSource extends CellColoring {
  meanValues(b: number): Float64Array | null;
  meanPyramid(b: number): ValuePyramid | null;
  /** Read positions raw, no index over the block (a streamed block holds only its regions' letters, ST3a). */
  readonly raw?: boolean;
  /**
   * Block indices [out[0], out[1]) holding every position block b colors (a hull; none: the block's own): the classes
   * are built over it and the shown piece, never the whole domain (a protein's letters on its chromosome, GM-02;
   * seqspec:viewer#performance).
   */
  extent?(b: number, out: Int32Array): void;
  /**
   * What the coloring reads beyond the theme (a render prop: the letters' `unobservedGray`), as a number: the classes
   * are built again when it changes. None: 0.
   */
  variant?(v: DrawView): number;
}

/** The first shown index of block b in the device column of shown k (linear). */
function columnStart(v: DrawView, b: number, k: number): number {
  const c = v.xd(b, k);
  let j = Math.max(0, Math.min(k, Math.ceil(((c + v.scrollDev - 0.5) / v.dpr - v.blockX0[b]) / (v.span[b] * v.cell))));
  while (j > 0 && v.xd(b, j - 1) === c) j--;
  while (j < k && v.xd(b, j) < c) j++;
  return j;
}
/** The first shown index past the device column of shown k (linear), at most the block's shown count. */
export function columnEnd(v: DrawView, b: number, k: number): number {
  const n = v.axis ? v.axis.shown[b] : 0;
  const c = v.xd(b, k);
  let j = Math.max(
    k + 1,
    Math.min(n, Math.ceil(((c + v.scrollDev + 0.5) / v.dpr - v.blockX0[b]) / (v.span[b] * v.cell))),
  );
  while (j > k + 1 && v.xd(b, j - 1) > c) j--;
  while (j < n && v.xd(b, j) <= c) j++;
  return j;
}
/**
 * Block b's shown range drawn here into `out`: wrapped, the view's line; linear, the view's, widened to whole device
 * columns where cells are thinner than one, so a reduction counts every position under a pixel; [0, 0): not here.
 */
export function visibleRange(v: DrawView, b: number, out: Int32Array): void {
  out[0] = out[1] = 0;
  if (v.onlyBlock >= 0 && b !== v.onlyBlock) return;
  if (v.geometry === AxisGeometry.wrapped) {
    out[0] = v.lineStart;
    out[1] = v.lineStop;
    return;
  }
  const lo = v.visLo[b];
  const hi = v.visHi[b];
  if (lo >= hi) return;
  const thin = v.cell * v.span[b] * v.dpr < 1;
  out[0] = thin ? columnStart(v, b, lo) : lo;
  out[1] = thin ? columnEnd(v, b, hi - 1) : hi;
}

const runs = new Int32Array(4);
/**
 * Bins [out[0], out[1]) at `level` of the column of shown [k0, k1): those whose start lies in its one block coordinate
 * system run; false when the column crosses a break or a cyclic origin (read it raw).
 */
export function columnBins(v: DrawView, b: number, k0: number, k1: number, level: number, out: Int32Array): boolean {
  if (v.shownToRuns(b, k0, k1, runs) !== 1) return false;
  const base = v.axis ? v.axis.base[b] : 0;
  const size = (1 << level) - 1;
  out[0] = (runs[0] - base + size) >> level;
  out[1] = (runs[1] - base + size) >> level;
  return true;
}
/** The block's hull of block indices over shown [k0, k1) into `out`; a cyclic origin inside: its whole domain. */
export function shownHull(v: { axis: Axis | null }, b: number, k0: number, k1: number, out: Int32Array): void {
  const ax = v.axis;
  if (!ax || k1 <= k0) {
    out[0] = out[1] = 0;
    return;
  }
  const w = ax.wrapAt(b);
  if (w > k0 && w < k1) {
    out[0] = ax.base[b];
    out[1] = ax.base[b] + ax.length[b];
    return;
  }
  const a = ax.indexOf(b, k0 as never);
  const c = ax.indexOf(b, (k1 - 1) as never);
  out[0] = Math.min(a, c);
  out[1] = Math.max(a, c) + 1;
}
let hullRuns = new Int32Array(32);
/**
 * The shown hull of block b's run [start, stop) of indices (`shownHull`'s inverse): its first shown position and one
 * past its last, into `out`; false: none of it is shown. Through the axis's `runToShown` (a filter's runs, a reversed
 * piece, a cyclic origin), module scratch grown only. Ints in and out.
 */
export function runHull(
  v: { runToShown(b: number, start: number, stop: number, out: Int32Array): number },
  b: number,
  start: number,
  stop: number,
  out: Int32Array,
): boolean {
  let n = v.runToShown(b, start, stop, hullRuns);
  if (2 * n > hullRuns.length) {
    hullRuns = new Int32Array(4 * n);
    n = v.runToShown(b, start, stop, hullRuns);
  }
  if (n === 0) return false;
  let lo = hullRuns[0];
  let hi = hullRuns[1];
  for (let r = 1; r < n; r++) {
    if (hullRuns[2 * r] < lo) lo = hullRuns[2 * r];
    if (hullRuns[2 * r + 1] > hi) hi = hullRuns[2 * r + 1];
  }
  out[0] = lo;
  out[1] = hi;
  return hi > lo;
}
/** The pyramid level for block b's device columns (positions per column), of `levels`; -1: raw. */
export function columnLevel(v: DrawView, b: number, levels: number, base = PYRAMID_BASE): number {
  const per = 1 / (v.cell * v.span[b] * v.dpr);
  if (!(per >= 2 ** base) || levels === 0) return -1;
  return Math.min(Math.floor(Math.log2(per)), base + levels - 1);
}

/**
 * Per (prepared state, theme), per block: its classes (null past 255) over block indices [base + offset, base + offset
 * + its length), its mean values and their pyramid.
 */
interface CellIndex {
  theme: Theme;
  /** The source's `variant` it was built for. */
  variant: number;
  classes: (CellClasses | null | undefined)[];
  /** Per block the classes' first index from the block's base (a multiple of their top bin) and their length. */
  span: Int32Array;
  /** Per block the classes' levels a column may read: those whose bins start at a multiple of their size. */
  levels: Int32Array;
  values: (Float64Array | null | undefined)[];
  pyramids: (ValuePyramid | null | undefined)[];
  /** Per block its values' [min, max]: a mean is clamped to them (a sum of equal values may pass them by an ulp). */
  range: Float64Array;
}
const indexes = new WeakMap<object, CellIndex>();
function cellIndex(owner: object, theme: Theme, blocks: number, variant: number): CellIndex {
  const x = indexes.get(owner);
  return x?.theme === theme && x.variant === variant && x.classes.length === blocks
    ? x
    : newIndex(owner, theme, blocks, variant);
}
/** Edge-rate, apart: a closure in `cellIndex` would make its context on every call. */
function newIndex(owner: object, theme: Theme, blocks: number, variant: number): CellIndex {
  const none = () => new Array<undefined>(blocks).fill(undefined);
  const x = {
    theme,
    variant,
    classes: none(),
    span: new Int32Array(2 * blocks),
    levels: new Int32Array(blocks),
    values: none(),
    pyramids: none(),
    range: new Float64Array(2 * blocks),
  };
  indexes.set(owner, x);
  return x;
}

/** Block b's mean values and their range (edge-rate). */
function valuesOf(x: CellIndex, src: CellSource, b: number) {
  const values = src.meanValues(b);
  x.values[b] = values;
  let [lo, hi] = [Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY];
  for (let i = 0; values && i < values.length; i++) if (values[i] < lo) lo = values[i];
  for (let i = 0; values && i < values.length; i++) if (values[i] > hi) hi = values[i];
  x.range[2 * b] = lo;
  x.range[2 * b + 1] = hi;
}

const ext = new Int32Array(2);
/**
 * Into `ext` the block indices block b's classes must cover, from its base: its shown piece's hull, within the
 * source's `extent`; [0, 0): none. Integers only (each aggregate frame checks the index still covers it).
 */
function classRange(v: DrawView, ax: Axis, src: CellSource, b: number): void {
  const base = ax.base[b];
  let lo = 0;
  let hi = ax.length[b];
  if (ax.shown[b] > 0) {
    shownHull(v, b, 0, ax.shown[b], ext);
    lo = Math.max(lo, ext[0] - base);
    hi = Math.min(hi, ext[1] - base);
  }
  if (src.extent) {
    src.extent(b, ext);
    lo = Math.max(lo, ext[0] - base);
    hi = Math.min(hi, ext[1] - base);
  }
  ext[0] = lo;
  ext[1] = hi > lo ? hi : lo;
}
/**
 * Block b's classes over `ext` (edge rate): from an offset aligned to 2^m, the least power of two at least their extent,
 * so they cover under three times it wherever it lies (an offset aligned to a power of two past the whole range fell
 * to 0 across a high one: LC-25); the levels a column reads are capped at m, whose bins start at multiples of their
 * size, so a column's bins shift by `offset >> level` (a coarser column sums several level-m bins).
 */
function classesOf(x: CellIndex, ax: Axis, src: CellSource, b: number): void {
  const lo = ext[0];
  const hi = ext[1];
  let m = PYRAMID_BASE;
  while (2 ** m < hi - lo) m++;
  const off = lo - (lo % 2 ** m);
  x.span[2 * b] = off;
  x.span[2 * b + 1] = hi - off;
  x.levels[b] = m - PYRAMID_BASE + 1;
  x.classes[b] = hi > lo ? cellClasses(src, b, ax.base[b] + off, hi - off) : null;
}

/** The cells' band (device px): top, height, and the gap right of each cell. Set by the caller. */
export const cellBand = new Int32Array(3);
const bins = new Int32Array(2);
const acc = new Float64Array(4);
/** A raw column's distinct (color, alpha) pairs counted: a protein's themes have up to 20 colors (GM-09: 8 dropped some). */
const SLOTS = 32;
const pairKey = new Float64Array(SLOTS);
const pairN = new Int32Array(SLOTS);

/** The sum and count of shown [k0, k1)'s values into `acc[2..3]` (a column read raw). */
function rawMean(v: DrawView, b: number, k0: number, k1: number, mean: Float64Array) {
  const ax = v.axis;
  const base = ax ? ax.base[b] : 0;
  for (let q = k0; q < k1; q++) {
    const val = mean[(ax ? ax.indexOf(b, q as never) : q) - base];
    if (!Number.isNaN(val)) {
      acc[2] += val;
      acc[3]++;
    }
  }
}

/** Shown [lo, hi) of block b as cells by the view's color theme, into `rects` (the caller hands it to the backend). */
export function drawCells(
  v: DrawView,
  b: number,
  lo: number,
  hi: number,
  src: CellSource,
  owner: object,
  rects: RectBatch,
) {
  const t = v.color;
  const ax = v.axis;
  if (!t || !ax || lo >= hi) return;
  const top = cellBand[0];
  const h = cellBand[1];
  if (v.geometry === AxisGeometry.wrapped || v.cell * v.span[b] * v.dpr >= 1) {
    for (let k = lo; k < hi; k++) {
      const i = ax.indexOf(b, k as never);
      const color = src.color(b, i);
      // A skin token (below -1) is a color too: the backend resolves it.
      if (color === NO_COLOR) continue;
      rects.begin(color, src.alpha(b, i));
      const x = v.xd(b, k);
      rects.push(x, top, Math.max(1, v.xe(b, k) - cellBand[2] - x), h);
    }
    return;
  }
  const x = cellIndex(owner, t, ax.count, src.variant ? src.variant(v) : 0);
  let mean: Float64Array | null = null;
  let pyr: ValuePyramid | null = null;
  let cls: CellClasses | null = null;
  if (t.reduction === 'mean') {
    if (x.values[b] === undefined) valuesOf(x, src, b);
    mean = x.values[b] as Float64Array | null;
    if (mean && x.pyramids[b] === undefined) x.pyramids[b] = src.meanPyramid(b) ?? valuePyramid(mean, false);
    pyr = mean ? (x.pyramids[b] as ValuePyramid) : null;
  }
  if (!mean && src.raw !== true) {
    classRange(v, ax, src, b);
    // Built, or built again where the piece now shows indices the classes do not cover.
    const o = x.span[2 * b];
    if (x.classes[b] === undefined || ext[0] < o || ext[1] > o + x.span[2 * b + 1]) classesOf(x, ax, src, b);
    cls = x.classes[b] as CellClasses | null;
  }
  const shift = cls ? x.span[2 * b] : 0;
  const levels =
    ax.positions[b] !== null ? 0 : pyr ? pyr.sum.length : cls ? Math.min(cls.pyramid.cls.length, x.levels[b]) : 0;
  const level = columnLevel(v, b, levels);
  for (let k = lo; k < hi; ) {
    const k1 = columnEnd(v, b, k);
    const x0 = v.xd(b, k);
    let color = -1;
    let alpha = OPAQUE;
    const binned = level >= 0 && columnBins(v, b, k, k1, level, bins);
    if (mean) {
      acc[2] = acc[3] = 0;
      if (binned) binValues(pyr as ValuePyramid, level, bins[0], bins[1], acc);
      else rawMean(v, b, k, k1, mean);
      if (acc[3] > 0) {
        const loc = ThemeLocation.at(v.loc, ax, v.layer, b, ax.indexOf(b, k as never));
        loc.annotation = loc.row = -1;
        loc.value = Math.min(x.range[2 * b + 1], Math.max(x.range[2 * b], acc[2] / acc[3]));
        color = t.value(loc);
        alpha = t.alpha ? t.alpha(loc) : OPAQUE;
        loc.value = Number.NaN;
      }
    } else if (binned && cls) {
      const c = binClass(cls.pyramid, level, bins[0] - (shift >> level), bins[1] - (shift >> level), acc);
      if (c > 0) {
        color = cls.color[c];
        alpha = cls.alpha[c];
      }
    } else {
      // The plurality (color, alpha) of the column's positions.
      let n = 0;
      let best = -1;
      for (let q = k; q < k1; q++) {
        const i = ax.indexOf(b, q as never);
        const c = src.color(b, i);
        if (c === NO_COLOR) continue;
        // A skin token's key is negative and decodes the same (floor division).
        const key = c * 256 + src.alpha(b, i);
        let s = 0;
        while (s < n && pairKey[s] !== key) s++;
        if (s === n && n < SLOTS) {
          pairKey[n] = key;
          pairN[n++] = 0;
        }
        if (s < n && ++pairN[s] > (best < 0 ? 0 : pairN[best])) best = s;
      }
      if (best >= 0) {
        color = Math.floor(pairKey[best] / 256);
        alpha = pairKey[best] - color * 256;
      }
    }
    if (color !== NO_COLOR) {
      rects.begin(color, alpha);
      rects.push(x0, top, Math.max(1, v.xd(b, k1) - x0), h);
    }
    k = k1;
  }
}
