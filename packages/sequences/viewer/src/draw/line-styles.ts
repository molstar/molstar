/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Line styles in the draws (seqspec:format#themes; seqspec:viewer#the-render-backend), part of the renderer API: an
// item's style from its layer's line style theme, and a styled segment over a band of device px. `solid` is the
// caller's own rectangle; the others go to a `StyledLines`' batches, one per style kind (as colors are batched): dashed
// and dotted a stroke with a dash pattern, wavy a sine polyline, both anchored to the content (the pattern and the
// wave move with the positions, never with a piece cut by the view), double two thin lines, box an outline, highlight
// a translucent fill drawn under the rest. Device px are integers in and doubles only inside; allocation-free once warm.

import { AxisGeometry } from '../axis.js';
import type { Features, StreamFields } from '../prepare.js';
import { LineStyle, type Theme } from '../themes.js';
import { PathBatch, PathMode, PathVerb, RectBatch, type RenderBackend } from './backend.js';
import { featureLoc } from './density.js';
import type { DrawView } from './view.js';

// Read once: a per-frame read of a module object's member allocated under the test runner's module wrapper.
const SOLID = LineStyle.solid;
const DASHED = LineStyle.dashed;
const DOTTED = LineStyle.dotted;
const WAVY = LineStyle.wavy;
const DOUBLE = LineStyle.double;
const BOX = LineStyle.box;
const HIGHLIGHT = LineStyle.highlight;
const STROKE = PathMode.stroke;
const WRAPPED = AxisGeometry.wrapped;
const MOVE = PathVerb.move;
const LINE = PathVerb.line;
/** A highlight's alpha over its color's (of 256). */
const HIGHLIGHT_ALPHA = 72;
/** A wave's length and a dashed or dotted line's widest, CSS px; samples per wave. */
const WAVE = 8;
const LINE_MAX = 3;
const WAVE_SAMPLES = 8;
const TAU = 2 * Math.PI;

/**
 * Feature k's line style (a `LineStyle`): its layer's line style theme at its location, a streamed item's from its
 * region's table (`Theme.raw`); solid without one. Leaves `v.loc` at k.
 */
export function featureLineStyle(v: DrawView, f: Features, k: number): number {
  const t = v.lineStyle;
  if (!t || !v.axis) return SOLID;
  const s = f.fields;
  const x =
    s != null && t.raw !== undefined
      ? s.styleTheme === t
        ? (s.styles as Int32Array)[k]
        : styleValues(f, s, t, k)
      : t.value(featureLoc(v, f, k));
  return x < 0 ? SOLID : x;
}
/** Edge rate (a region's first draw under a theme): its items' styles by `t`; feature k's. */
function styleValues(f: Features, s: StreamFields, t: Theme, k: number): number {
  const raw = t.raw as NonNullable<Theme['raw']>;
  s.styles = Int32Array.from({ length: f.count }, (_, q) => raw(s.table, f.row[q]));
  s.styleTheme = t;
  return s.styles[k];
}

/** A renderer module's styled marks: one of each batch, reused by every draw. */
export class StyledLines {
  /** Highlights, drawn first: under the layer's other marks. */
  readonly under = new RectBatch(16);
  /** Doubles' lines and boxes' outlines. */
  readonly rects = new RectBatch(16);
  /** Dashed, dotted and wavy lines. */
  readonly strokes = new PathBatch(16);
  reset(): void {
    this.under.reset();
    this.rects.reset();
    this.strokes.reset();
  }
  /** The highlights, before the layer's own marks. */
  drawUnder(backend: RenderBackend): void {
    if (this.under.n > 0) backend.rects(this.under);
  }
  /** The other styles, after the layer's own rectangles. */
  draw(backend: RenderBackend): void {
    if (this.rects.n > 0) backend.rects(this.rects);
    if (this.strokes.n > 0) backend.paths(this.strokes);
  }
}

/**
 * Scratch, so no double crosses a call: [0] the content's anchor x (device px) of the block drawn, [1] a dash's phase,
 * [2] a wave's length (device px).
 */
const num = new Float64Array(3);
/** Block b's content anchor into `num[0]` (device px): its position 0 linear, the line's left wrapped. */
function anchor(v: DrawView, b: number): void {
  num[0] = v.geometry === WRAPPED ? (v.left + v.lineLeft) * v.dpr : v.blockX0[b] * v.dpr - v.scrollDev;
}

/**
 * Feature ink `ink` (alpha `alpha`) in `style` over [x0, x0 + w) × [y0, y1) of block b (device px, integers): a styled
 * line along the band (its middle for dashed, dotted and wavy, its edges for double), its outline (box) or its fill
 * (highlight). False: solid, or a band too thin for the style: the caller draws its rectangle.
 */
export function styledSegment(
  v: DrawView,
  out: StyledLines,
  style: number,
  b: number,
  x0: number,
  w: number,
  y0: number,
  y1: number,
  ink: number,
  alpha: number,
): boolean {
  if (style === SOLID || style < 0 || w <= 0) return false;
  const h = y1 - y0;
  const px = v.dpr >= 1 ? (v.dpr + 0.5) | 0 : 1;
  if (style === HIGHLIGHT) {
    out.under.begin(ink, (alpha * HIGHLIGHT_ALPHA) >> 8);
    out.under.push(x0, y0, w, h);
    return true;
  }
  if (style === BOX) {
    if (h < 2 * px + 1 || w < 2 * px + 1) return false;
    out.rects.begin(ink, alpha, px);
    out.rects.push(x0, y0, w, h);
    return true;
  }
  if (style === DOUBLE) {
    if (h < 3 * px) return false;
    out.rects.begin(ink, alpha);
    out.rects.push(x0, y0, w, px);
    out.rects.push(x0, y1 - px, w, px);
    return true;
  }
  if (style !== DASHED && style !== DOTTED && style !== WAVY) return false;
  anchor(v, b);
  const s = out.strokes;
  if (style === WAVY) {
    // A sine of a fixed length in px, its amplitude what the band holds past the line's width (at most 3 CSS px).
    const lw = h >= 3 * px ? px : 1;
    let amp = (h - lw) >> 1;
    if (amp > 3 * px) amp = 3 * px;
    if (amp < 1) return false;
    s.begin(ink, alpha, STROKE, 2 * lw);
    num[2] = WAVE * v.dpr;
    wave(s, x0, x0 + w, y0 + y1, amp);
    return true;
  }
  // A line through the band's middle, at most 3 CSS px wide; the pattern in units of its width.
  let lw = h < LINE_MAX * px ? h : LINE_MAX * px;
  if (lw > 127) lw = 127;
  if (lw < 1) lw = 1;
  const on = style === DASHED ? 3 * lw : lw;
  const gap = style === DASHED ? 2 * lw : lw;
  num[1] = (x0 - num[0]) % (on + gap);
  if (num[1] < 0) num[1] += on + gap;
  s.begin(ink, alpha, STROKE, 2 * lw);
  // `dashed`'s fields written here: its phase a double.
  s.dash[s.n - 1] = on;
  s.gap[s.n - 1] = gap;
  s.phase[s.n - 1] = num[1];
  let o = s.verb(MOVE);
  s.coords[o] = x0;
  s.coords[o + 1] = (y0 + y1) / 2;
  o = s.verb(LINE);
  s.coords[o] = x0 + w;
  s.coords[o + 1] = (y0 + y1) / 2;
  return true;
}

/**
 * A sine polyline from x0 to x1 about `twiceY / 2` (device px), `num[2]` px a wave, its phase from the content's anchor
 * (`num[0]`), sampled at a fixed grid of the wave (`WAVE_SAMPLES` a length) and at both ends.
 */
function wave(s: PathBatch, x0: number, x1: number, twiceY: number, amp: number): void {
  const a = num[0];
  const length = num[2];
  const y = twiceY / 2;
  const step = length / WAVE_SAMPLES;
  let o = s.verb(MOVE);
  s.coords[o] = x0;
  s.coords[o + 1] = y + amp * Math.sin((TAU * (x0 - a)) / length);
  let x = a + Math.ceil((x0 - a) / step) * step;
  if (x <= x0) x += step;
  for (; x < x1; x += step) {
    o = s.verb(LINE);
    s.coords[o] = x;
    s.coords[o + 1] = y + amp * Math.sin((TAU * (x - a)) / length);
  }
  o = s.verb(LINE);
  s.coords[o] = x1;
  s.coords[o + 1] = y + amp * Math.sin((TAU * (x1 - a)) / length);
}
