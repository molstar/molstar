/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `frame` renderer (seqspec:viewer#frames, seqspec:format#layer-types): named regions outlined across the tracks
// its layer spans, an outline per item and run (`FrameFeatures`), drawn by `across` after the entries; its own band
// holds each item's title on the line (or in the view) its item starts. A side is drawn only at its item's own end:
// where a line break or the view's edge cuts the item, that side stays open. Edges on the cells' edges; a dashed or
// dotted style's pattern fixed to the content. Batches are this module's, reset per draw: nothing is allocated a frame.

import {
  type AxisView,
  AxisGeometry,
  type DrawView,
  type FrameFeatures,
  featureInk,
  featureKey,
  featureLineStyle,
  fitIn,
  fitted,
  keyAlpha,
  kernelEntry,
  type LineTable,
  LineStyle,
  PathBatch,
  PathMode,
  PathVerb,
  runHull,
  TextAlign,
  TextBaseline,
  TextBatch,
  ThemeLocation,
} from '../draw/api.js';

// Read once: a module member read per frame allocates under the test runner's module wrapper.
const WRAPPED = AxisGeometry.wrapped;
const STROKE = PathMode.stroke;
const MOVE = PathVerb.move;
const LINE = PathVerb.line;
const DASHED = LineStyle.dashed;
const DOTTED = LineStyle.dotted;

const outline = new PathBatch(16);
const titles = new TextBatch(16);
/** An item's shown hull [first, end) on its block. */
const hull = new Int32Array(2);
const fit = new Float64Array(2);

interface Params {
  readonly labels?: boolean;
}

/** A title's font px for a band of `h` CSS px: 0.7 of the text's (`fontPx`), within the band. */
const titlePx = (fontPx: number, h: number) => Math.max(7, Math.min(Math.round(fontPx * 0.7), Math.floor(h * 0.8)));
/** The title's row in the text: the track gap, then its line height. */
const titleRow = (fontPx: number, gap: number) => gap + Math.ceil(Math.max(7, Math.round(fontPx * 0.7)) * 1.4);

/** Item k's shown hull on its block into `hull`; false: none of it shown in the block. */
const hullOf = (v: AxisView, f: FrameFeatures, k: number) =>
  runHull(v, f.block[k], f.start[k], Math.max(f.stop[k], f.start[k] + 1), hull);

/**
 * Lane offset k's display row (`v.spanDisplay`, the rows' `RowModel.display`): inside a collapsed clade, its wedge's row
 * (the clade's first lane holds it; the others -1, contiguous after it); -1 where none is shown.
 */
function displayOf(disp: Int32Array, k: number): number {
  let q = k < disp.length ? k : disp.length - 1;
  while (q > 0 && disp[q] < 0) q--;
  return q >= 0 ? disp[q] : -1;
}

/** The shown range drawn here: wrapped the line's, linear the whole view's (a pan's strip decides as a full draw). */
const loOf = (v: DrawView, b: number) => (v.geometry === WRAPPED ? v.lineStart : v.viewLo[b]);
const hiOf = (v: DrawView, b: number) => (v.geometry === WRAPPED ? v.lineStop : v.viewHi[b]);
/** The right edge of shown k − 1's cell (device px): wrapped its advance's, linear the next edge. */
const rightOf = (v: DrawView, b: number, k: number) => (v.geometry === WRAPPED ? v.xe(b, k - 1) : v.xd(b, k));

/** The titles: each item's label where it starts here, at its left edge, up to the next title or the line's end. */
function draw(v: DrawView, f: FrameFeatures, p: Params) {
  if (p.labels === false || !v.axis) return;
  titles.reset();
  ThemeLocation.reset(v.loc);
  const dpr = v.dpr;
  const font = v.text.font('ui', titlePx(v.geometry === WRAPPED ? v.fontPx : v.props.labelFontPx, v.height), 600);
  titles.style(font, TextAlign.left, TextBaseline.bottom);
  fit[1] = v.text.label('…', font).width;
  const y = Math.round((v.top + v.height - 1) * dpr) | 0;
  for (let k = 0; k < f.count; k++) {
    const label = f.label[k];
    const b = f.block[k];
    if (label === null || (v.onlyBlock >= 0 && b !== v.onlyBlock) || !hullOf(v, f, k)) continue;
    const lo = loOf(v, b);
    const hi = hiOf(v, b);
    if (hull[0] < lo || hull[0] >= hi) continue;
    const x = v.xd(b, hull[0]) + Math.round(2 * dpr);
    // Up to the next title on this block, else the line's end (the view's right edge). The items are sorted by start
    // (prepare.ts): on a forward block their shown starts ascend, so past the line none later is in it; a reversed
    // block's shown start is its item's stop side, in no order.
    let end = rightOf(v, b, hi);
    const forward = v.dir[b] > 0;
    for (let j = k + 1; j < f.count && f.block[j] === b; j++) {
      if (f.label[j] === null || !hullOf(v, f, j)) continue;
      if (hull[0] >= hi) {
        if (forward) break;
        continue;
      }
      if (hull[0] > lo) {
        end = Math.min(end, v.xd(b, hull[0]) - Math.round(4 * dpr));
        break;
      }
    }
    const pl = v.text.label(label, font);
    fit[0] = (end - x) / dpr;
    const m = fitIn(pl, fit);
    if (m === 0 || m < Math.min(3, pl.ends.length)) continue;
    const key = featureKey(v, f, k);
    titles.begin(featureInk(key), keyAlpha(key));
    const t = titles.add(fitted(pl, m));
    titles.x[t] = x;
    titles.y[t] = y;
  }
  if (titles.n > 0) v.backend.text(titles);
}

/** A segment of the outline (device px, integers: no double crosses the call). */
function segment(
  ink: number,
  alpha: number,
  w: number,
  style: number,
  phase: number,
  ax: number,
  ay: number,
  bx: number,
  by: number,
) {
  outline.begin(ink, alpha, STROKE, w);
  if (style === DASHED || style === DOTTED) {
    const lw = Math.max(1, w >> 1);
    const on = style === DASHED ? 3 * lw : lw;
    const gap = style === DASHED ? 2 * lw : lw;
    outline.dashed(on, gap, phase % (on + gap));
  }
  const a = outline.verb(MOVE);
  outline.coords[a] = ax;
  outline.coords[a + 1] = ay;
  const z = outline.verb(LINE);
  outline.coords[z] = bx;
  outline.coords[z + 1] = by;
}

/**
 * Each item's outline over the span (`v.top`, `v.height`), its top or bottom at the rows its item targets where a rows
 * band is that edge (`v.spanRows`); a side only at its item's own end (seqspec:viewer#frames).
 */
function across(v: DrawView, f: FrameFeatures) {
  if (!v.axis) return;
  outline.reset();
  ThemeLocation.reset(v.loc);
  const dpr = v.dpr;
  // 1.5 CSS px rounded to whole device px (a width in half device px: even).
  const w = (2 * Math.max(1, Math.round(1.5 * dpr))) | 0;
  const sr = v.spanRows;
  const disp = v.spanDisplay;
  for (let k = 0; k < f.count; k++) {
    const b = f.block[k];
    if ((v.onlyBlock >= 0 && b !== v.onlyBlock) || !hullOf(v, f, k)) continue;
    const lo = loOf(v, b);
    const hi = hiOf(v, b);
    if (hull[1] <= lo || hull[0] >= hi) continue;
    const s0 = Math.max(lo, hull[0]);
    const s1 = Math.min(hi, hull[1]);
    const x0 = v.xd(b, s0) | 0;
    const x1 = rightOf(v, b, s1) | 0;
    let top = v.top;
    let bottom = v.top + v.height;
    const dLo = f.laneLo[k] >= 0 && disp !== null ? displayOf(disp, f.laneLo[k]) : f.laneLo[k];
    const dHi = f.laneHi[k] >= 0 && disp !== null ? displayOf(disp, f.laneHi[k]) : f.laneHi[k];
    if (sr[2] === 1 && dLo >= 0) top = sr[0] + dLo * sr[1];
    if (sr[3] === 1 && dHi >= 0) bottom = sr[0] + (dHi + 1) * sr[1];
    const y0 = Math.round(top * dpr) | 0;
    const y1 = Math.round(bottom * dpr) | 0;
    if (y1 <= y0) continue;
    const key = featureKey(v, f, k);
    const ink = featureInk(key);
    const alpha = keyAlpha(key);
    const style = featureLineStyle(v, f, k);
    // The pattern from the item's start: linear, the cut's offset along the item (a pan keeps it on the content).
    const phase = v.geometry === WRAPPED ? 0 : Math.round((s0 - hull[0]) * v.cell * v.span[b] * dpr) | 0;
    segment(ink, alpha, w, style, phase, x0, y0, x1, y0);
    segment(ink, alpha, w, style, phase, x0, y1, x1, y1);
    if (hull[0] >= lo) segment(ink, alpha, w, style, 0, x0, y0, x0, y1);
    if (hull[1] <= hi) segment(ink, alpha, w, style, 0, x1, y0, x1, y1);
  }
  if (outline.n > 0) v.backend.paths(outline);
}

/**
 * Wrapped: a title's row on each line where a titled item starts (`labels`), else none (seqspec:viewer#frames).
 * Edge rate: at a layout.
 */
function lineHeights(s: unknown, p: unknown, axis: AxisView, lines: LineTable, out: Float32Array) {
  const f = s as FrameFeatures | undefined;
  for (let i = 0; i < lines.count; i++) out[i] = 0;
  if (!f || (p as Params).labels === false) return;
  const h = titleRow(axis.fontPx, axis.trackGap);
  for (let k = 0; k < f.count; k++) {
    if (f.label[k] === null || !hullOf(axis, f, k)) continue;
    for (let i = 0; i < lines.count; i++)
      if (lines.block[i] === f.block[k] && hull[0] >= lines.start[i] && hull[0] < lines.stop[i]) out[i] = h;
  }
}

/** In tracks, its titles' row; without titles (`labels: false`) none (the layout gives the track no band). */
const height = (_s: unknown, p: unknown, props: { rowHeight: number }) =>
  (p as Params).labels === false ? 0 : props.rowHeight;

export const FrameRenderer = kernelEntry<FrameFeatures>('frame', height, {
  geometries: ['linear', 'wrapped'],
  content: false,
  miniature: false,
  draw: (v, f, p) => draw(v, f, p as Params),
  tiny: (v, f, p) => draw(v, f, p as Params),
  across: (v, f) => across(v, f),
  lineHeights,
});
