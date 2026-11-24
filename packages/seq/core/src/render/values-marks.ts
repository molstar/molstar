/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { perPosition } from './config.js';
import { NoColor, type PassWriter } from './draw-lists.js';
import { drawFailed, drawSpans } from './marks.js';
import type { ValuesLay } from './values.js';

// The values layer's dots and its `points` encoding (_spec/README.md "Render", Layer types, `values`): a dot is a small
// polygon (a circle to the eye: the squircle's corner is no circle) of the value's color, in a ring of the well's
// ground that sets it off the line and the other dots. A line or an area shows its dots from `DotMinPx` CSS px a value;
// `points` draws them at every value of a per-position step, and below it (a value under a pixel) the device column's
// extent as a thin vertical mark, min to max. A missing value has no dot and no mark. Nothing allocated once warm: the
// numbers go through `out.args`, the polygon through `ring`.

/** The share of the height the domain's top takes: the rest is headroom (10%). */
export const Fill = 0.9;

/**
 * A value on the y scale: its log10 on a log scale (a value at or below 0 under the domain, at its low end), else
 * itself. Heights are linear in it, so a log scale's decades are equal steps (EC-215).
 */
export function ys(lay: Pick<ValuesLay, 'yLog'>, v: number): number {
  return lay.yLog ? (v > 0 ? Math.log10(v) : Number.NEGATIVE_INFINITY) : v;
}

/** The y domain's low and high ends (`[0, 1]` before one is computed), in values. */
export const domainLo = (lay: Pick<ValuesLay, 'domain'>) => (lay.domain ? lay.domain[0] : 0);
export const domainHi = (lay: Pick<ValuesLay, 'domain'>) => (lay.domain ? lay.domain[1] : 1);

/** The baseline, in values: zero inside the domain, else its low end. Heights are `ys` of values over the scaled span. */
export function baseOf(lay: Pick<ValuesLay, 'domain'>): number {
  const d0 = domainLo(lay);
  return d0 <= 0 && domainHi(lay) >= 0 ? 0 : d0;
}

/** A line and an area show a dot a value from this many CSS px a value. */
export const DotMinPx = 6;
/** A dot's radius (3 px across) and the ring's width round it, CSS px. */
export const DotRadius = 1.5;
export const DotRing = 1;
const Sides = 16;
const Cos = Float64Array.from({ length: Sides }, (_, i) => Math.cos((2 * Math.PI * i) / Sides));
const Sin = Float64Array.from({ length: Sides }, (_, i) => Math.sin((2 * Math.PI * i) / Sides));
const ring = new Float64Array(2 * Sides);

/** The polygon of `Sides` corners around (x, y) at radius r, written to `ring`. */
function polygon(x: number, y: number, r: number) {
  for (let i = 0; i < Sides; i++) {
    ring[2 * i] = x + r * Cos[i];
    ring[2 * i + 1] = y + r * Sin[i];
  }
}

const scale = new Float64Array(3);

/** The y scale's low end and high end (`ys` of the domain's) and its span with the headroom, to `into` (no return). */
function yScale(lay: ValuesLay, into: Float64Array) {
  const lo = ys(lay, domainLo(lay));
  const hi = ys(lay, domainHi(lay));
  into[0] = lo;
  into[1] = (hi - lo || 1) / Fill;
  into[2] = hi;
}

/** The hairline at the baseline (zero, else the domain's low), a device pixel thick. No number is passed to a call. */
export function drawBaseline(lay: ValuesLay, out: PassWriter) {
  const { box } = lay;
  const lo = ys(lay, domainLo(lay));
  const hi = ys(lay, domainHi(lay));
  const base = baseOf(lay);
  const span = (hi - lo || 1) / Fill;
  const yBase = box.y + box.height - ((Math.max(lo, Math.min(hi, ys(lay, base))) - lo) / span) * box.height;
  const px = 1 / lay.axis.dpr;
  const a = out.args;
  a[0] = box.x;
  a[1] = yBase - px;
  a[2] = box.width;
  a[3] = px;
  a[4] = lay.tokens.colors['seq-sequence'];
  a[5] = NoColor;
  a[6] = 0;
  a[7] = 0;
  a[8] = 0;
  out.rectArgs();
}

/** A dot a value, in x order: its ring of the well's ground, then the dot in the value's color. */
export function drawDots(lay: ValuesLay, out: PassWriter) {
  const { box } = lay;
  const lo = ys(lay, domainLo(lay));
  const hi = ys(lay, domainHi(lay));
  const span = (hi - lo || 1) / Fill;
  const y1 = box.y + box.height;
  const bh = box.height;
  const ground = lay.tokens.colors['seq-track'];
  const a = out.args;
  const ord = lay.ord;
  for (let k = 0; k < lay.n; k++) {
    const j = ord[k];
    const v = lay.v[j];
    if (Number.isNaN(v)) continue;
    const y = y1 - ((Math.max(lo, Math.min(hi, ys(lay, v))) - lo) / span) * bh;
    const x = (lay.x0[j] + lay.x1[j]) / 2;
    polygon(x, y, DotRadius + DotRing);
    a[4] = ground;
    a[5] = NoColor;
    a[6] = 0;
    a[7] = 0;
    out.pathArgs(ring, Sides, true);
    polygon(x, y, DotRadius);
    a[4] = lay.col[j];
    out.pathArgs(ring, Sides, true);
  }
}

/** The data pass of `points`: placeholders, the baseline, then dots (a per-position step) or the columns' extents. */
export function drawPoints(lay: ValuesLay, out: PassWriter) {
  drawSpans(out, lay.pending, lay.box.y, lay.box.height, lay.tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, lay.box.y, lay.box.height, lay.tokens);
  drawBaseline(lay, out);
  if (perPosition(lay.step)) {
    drawDots(lay, out);
    return;
  }
  const { box } = lay;
  // the scale's low end and span through an array, never a call's return (a settled domain's double, returned by a call
  // that is not inlined here, is a heap number a frame)
  yScale(lay, scale);
  const lo = scale[0];
  const span = scale[1];
  const hi = scale[2];
  const y1 = box.y + box.height;
  const bh = box.height;
  const px = 1 / lay.axis.dpr;
  const a = out.args;
  for (let j = 0; j < lay.n; j++) {
    const v = lay.v[j];
    if (Number.isNaN(v)) continue;
    // a column's extent, min to max, at least a device pixel tall (a single value is a pixel, not nothing)
    const top = y1 - ((Math.max(lo, Math.min(hi, ys(lay, v))) - lo) / span) * bh;
    const bottom = y1 - ((Math.max(lo, Math.min(hi, ys(lay, lay.lo[j]))) - lo) / span) * bh;
    a[0] = lay.x0[j];
    a[1] = top;
    a[2] = lay.x1[j] - lay.x0[j];
    a[3] = Math.max(bottom - top, px);
    a[4] = lay.col[j];
    a[5] = NoColor;
    a[6] = 0;
    a[7] = 0;
    a[8] = 0;
    out.rectArgs();
  }
}

/** The hover's ring on the dot of bar `j`, outside it by its ring and the offset. */
export function hoverDot(lay: ValuesLay, out: PassWriter, j: number) {
  const { box } = lay;
  const lo = ys(lay, domainLo(lay));
  const hi = ys(lay, domainHi(lay));
  const span = (hi - lo || 1) / Fill;
  const y = box.y + box.height - ((Math.max(lo, Math.min(hi, ys(lay, lay.v[j]))) - lo) / span) * box.height;
  polygon((lay.x0[j] + lay.x1[j]) / 2, y, DotRadius + DotRing + lay.hoverOffset + lay.hoverWidth / 2);
  const a = out.args;
  a[4] = NoColor;
  a[5] = lay.tokens.colors['seq-hover'];
  a[6] = lay.hoverWidth;
  a[7] = 0;
  out.pathArgs(ring, Sides, true);
}
