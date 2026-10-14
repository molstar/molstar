/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `values` renderer (seqspec:viewer#extension-points): bars, a line or an area over the resolved domain, a cell
// each; below a device px per position, per device column the min–max band (light) and the mean as the bar, the line or
// the area (seqspec:sources#levels-and-the-zoom-gate; BigWig's summaries), from the values' prepared pyramid at the
// level the column holds (the bins whose start lies in it), raw below its base, in a filtered block and across a break.
// Every encoding in the color theme at the value drawn, a column's mean where columns aggregate
// (seqspec:format#layer-types): a bar, an area's column and the band their own; a line each position's over its own
// cell, so it changes color halfway between two positions, where their cells meet. With no theme, the ink.
// The overlay outlines the hovered position's bar, or rings its point on a line or an area (seqspec:interaction#hover).
import {
  AxisGeometry,
  binValues,
  columnBins,
  columnEnd,
  columnLevel,
  type DrawView,
  hoverEmphasis,
  hoverWidth,
  kernelEntry,
  NO_COLOR,
  NO_SIZE,
  outlineItem,
  PathBatch,
  PathMode,
  PathVerb,
  type PerPosition,
  RectBatch,
  ringItem,
  SIZE_ONE,
  SkinToken,
  sizeTheme,
  streamColumn,
  streamedBlock,
  ThemeLocation,
  visibleRange,
} from "../draw/api.ts";

const rects = new RectBatch();
/** The band and the area, translucent, under the bars (a backend groups a batch by color and alpha). */
const under = new RectBatch();
const paths = new PathBatch();
const range = new Int32Array(2);
const bins = new Int32Array(2);
/** A column's or cell's minimum, maximum, sum, count and mean (a double crosses no call: V8 boxes it). */
const stats = new Float64Array(5);
/** The band's bottom and height, device px, and the domain [lo, hi). */
const band = new Float64Array(4);
/** The line's last point (device px): held here, not in loop locals, which an unoptimized tier boxes. */
const prev = new Float64Array(2);
const BAND = 89;
const AREA = 77;
// Enum-like members read once (a per-frame read allocates under the test runner's module wrapper).
const LINE = PathVerb.line;
const MOVE = PathVerb.move;
const STROKE = PathMode.stroke;
const INK = SkinToken.ink;

/** `stats[j]`'s height above the band's bottom, device px (one at least). */
const heightOf = (j: number) =>
  Math.max(1, Math.round(Math.min(1, Math.max(0, (stats[j] - band[2]) / (band[3] - band[2]))) * band[1]));
/** The theme's color at `stats[4]` (none: the ink), its location filled for block b's index i. */
function colorAt(v: DrawView, b: number, i: number): number {
  const t = v.color;
  if (!v.axis) return INK;
  const loc = ThemeLocation.at(v.loc, v.axis, v.layer, b, i);
  loc.value = stats[4];
  if (!t) return INK;
  const c = t.value(loc);
  return c === NO_COLOR ? INK : c;
}

/** A bar `w` device px wide narrowed to the thickness channel's share of its cell (`size`, a theme's int), centered. */
function narrowed(w: number, size: number): number {
  return size < 0 || size >= SIZE_ONE ? w : Math.max(1, (w * size + (SIZE_ONE >> 1)) >> 8);
}

function draw(v: DrawView, s: PerPosition, p: unknown) {
  const ax = v.axis;
  if (!ax || v.geometry !== AxisGeometry.linear) return;
  const encoding = (p as { encoding?: string }).encoding;
  const bars = encoding !== "line" && encoding !== "area";
  const dpr = v.dpr;
  const top = Math.round((v.top + (v.height > 10 ? 3 : 0)) * dpr);
  band[0] = Math.round((v.top + v.height - (v.height > 10 ? 2 : 0)) * dpr);
  band[1] = band[0] - top;
  band[2] = s.lo;
  band[3] = s.hi;
  const bottom = band[0];
  const thickness = sizeTheme(v, 0, 1);
  rects.reset();
  under.reset();
  paths.reset();
  ThemeLocation.reset(v.loc);
  const width = Math.max(1, Math.round(3 * dpr));
  for (let b = 0; b < s.value.length; b++) {
    const arr = s.value[b];
    visibleRange(v, b, range);
    // ST1: a streamed block reads the stream manager's regions, a column or a cell at a time.
    if ((!arr && !streamedBlock(v, b)) || range[0] >= range[1]) continue;
    const base = ax.base[b];
    const thin = v.cell * v.span[b] * dpr < 1;
    let pen = false;
    let last = -1;
    // The line's last color (its last point is `prev`).
    let lc = 0;
    for (let k = range[0]; k < range[1]; ) {
      const k1 = thin ? columnEnd(v, b, k) : k + 1;
      const x = v.xd(b, k);
      const w = thin ? 1 : v.xe(b, k) - x;
      const i = ax.indexOf(b, k as never);
      if (!arr) streamColumn(v, b, k, k1, stats);
      else if (thin) columnStats(v, s, b, k, k1);
      else {
        stats[0] = stats[1] = stats[2] = arr[i - base];
        stats[3] = Number.isNaN(stats[0]) ? 0 : 1;
      }
      k = k1;
      if (stats[3] === 0) {
        pen = false;
        continue;
      }
      stats[4] = stats[2] / stats[3];
      const h = heightOf(4);
      const fill = colorAt(v, b, i);
      // The thickness channel (a bar's share of its cell): the theme reads the location `colorAt` left at the value.
      const size = bars && thickness && v.axis ? thickness.value(v.loc) : NO_SIZE;
      if (thin) {
        under.begin(fill, BAND);
        under.push(x, bottom - heightOf(1), 1, Math.max(1, heightOf(1) - heightOf(0)));
      }
      if (bars) {
        rects.begin(fill);
        const bw = Math.max(1, w - (w > 5 * dpr ? 1 : 0));
        const tw = thin ? bw : narrowed(bw, size);
        rects.push(x + ((bw - tw) >> 1), bottom - h, tw, h);
        continue;
      }
      if (encoding === "area") {
        under.begin(fill, AREA);
        under.push(x, bottom - h, Math.max(1, w), h);
      }
      // The line: a new subpath after a gap (a missing value, or columns more than 2 px apart), and where its color
      // changes, halfway to the point (where the two cells meet). Coordinates are written in place: no double crosses
      // a call.
      const px = x + w / 2;
      const py = bottom - h;
      const join = pen && (!thin || x - last <= 2 * dpr);
      if (join && fill !== lc) {
        const o = paths.verb(LINE);
        paths.coords[o] = (prev[0] + px) / 2;
        paths.coords[o + 1] = (prev[1] + py) / 2;
        paths.begin(fill, 255, STROKE, width);
        const m = paths.verb(MOVE);
        paths.coords[m] = paths.coords[o];
        paths.coords[m + 1] = paths.coords[o + 1];
      } else if (!join) paths.begin(fill, 255, STROKE, width);
      const o = paths.verb(join ? LINE : MOVE);
      paths.coords[o] = px;
      paths.coords[o + 1] = py;
      pen = true;
      last = x;
      prev[0] = px;
      prev[1] = py;
      lc = fill;
    }
  }
  if (under.n > 0) v.backend.rects(under);
  if (rects.n > 0) v.backend.rects(rects);
  if (paths.n > 0) v.backend.paths(paths);
}

/**
 * `stats` for the device column of shown [k0, k1) of block b: minimum, maximum, sum and count, from the pyramid's
 * bins whose start lies in it, or raw positions (below its base, filtered, or across a break).
 */
function columnStats(v: DrawView, s: PerPosition, b: number, k0: number, k1: number) {
  const ax = v.axis;
  const pyr = s.pyramid[b];
  const level = ax && pyr && ax.positions[b] === null ? columnLevel(v, b, pyr.min.length) : -1;
  if (pyr && level >= 0 && columnBins(v, b, k0, k1, level, bins)) {
    binValues(pyr, level, bins[0], bins[1], stats);
    return;
  }
  const arr = s.value[b] as Float64Array;
  const base = ax ? ax.base[b] : 0;
  stats[0] = stats[1] = Number.NaN;
  stats[2] = stats[3] = 0;
  for (let k = k0; k < k1; k++) {
    const val = arr[(ax ? ax.indexOf(b, k as never) : k) - base];
    if (Number.isNaN(val)) continue;
    if (!(val >= stats[0])) stats[0] = val;
    if (!(val <= stats[1])) stats[1] = val;
    stats[2] += val;
    stats[3]++;
  }
}

const el = new Int32Array(4);
/**
 * The hovered position (the layer's own hit, `marks.element`) where a position is a device px or more: its bar outlined
 * inside in its color's hover ink, or on a line or an area a ring at its point where the ring fits its column
 *.
 */
function overlay(v: DrawView, s: PerPosition, p: unknown) {
  const ax = v.axis;
  if (!hoverEmphasis(v) || !ax || v.geometry !== AxisGeometry.linear || !v.marks.element(v.layer, el)) return;
  const b = el[0];
  const arr = b < s.value.length ? s.value[b] : null;
  const k = v.shownOf(b, el[1]);
  if (!arr || k < 0 || v.cell * v.span[b] * v.dpr < 1) return;
  const value = arr[el[1] - ax.base[b]];
  if (Number.isNaN(value)) return;
  const dpr = v.dpr;
  band[0] = Math.round((v.top + v.height - (v.height > 10 ? 2 : 0)) * dpr);
  band[1] = band[0] - Math.round((v.top + (v.height > 10 ? 3 : 0)) * dpr);
  band[2] = s.lo;
  band[3] = s.hi;
  ThemeLocation.reset(v.loc);
  stats[4] = value;
  const h = heightOf(4) | 0;
  const x = v.xd(b, k);
  const w = v.xe(b, k) - x;
  const bottom = band[0] | 0;
  const encoding = (p as { encoding?: string }).encoding;
  if (encoding !== "line" && encoding !== "area") {
    const key = colorAt(v, b, el[1]);
    const t = sizeTheme(v, 0, 1);
    const bw = Math.max(1, w - (w > 5 * dpr ? 1 : 0));
    const tw = narrowed(bw, t ? t.value(v.loc) : NO_SIZE);
    outlineItem(v, key, x + ((bw - tw) >> 1), bottom - h, tw, h);
    return;
  }
  // A ring at the point in its color, only where it stays inside the position's column.
  const r = Math.max(2, Math.round(2 * dpr)) | 0;
  if (r + hoverWidth(dpr) + 1 <= w >> 1) ringItem(v, colorAt(v, b, el[1]), 0, x + (w >> 1), bottom - h, r, w >> 1);
}

export const ValuesRenderer = kernelEntry<PerPosition>("values", 40, {
  draw,
  tiny: draw,
  overlay,
  sized: ["thickness"],
});
