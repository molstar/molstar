/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the feature draws share (seqspec:viewer#extension-points `ranges`, `points`, `variants`, `secondary`;
// seqspec:sources#levels-and-the-zoom-gate), part of the renderer API: a feature's ink, a lane's rows in device px (the
// hit's `laneHeight`, `lanePad`,
// `laneInset`), the density band over an entry, and a wrapped line's underlines on the lanes `lineLanes` compacts to
// (the same compaction as the band's height); the hovered item's outlines on those underlines
// (seqspec:interaction#hover). Allocation-free once warm: module scratch, grown only.

import { type Axis, AxisGeometry, type AxisView } from "../axis.js";
import type { Hit } from "../hit.js";
import { type Features, featureAt } from "../prepare.js";
import { segWindow } from "../pyramids.js";
import { hitItem, laneHeight, laneInset, lanePad } from "../renderer-helpers.js";
import { groupOf, type LineTable, lineLanes } from "../runs.js";
import { NO_COLOR, RectBatch, SkinToken } from "./backend.js";
import { shownHull } from "./cells.js";
import { densityKindOf, drawDensity, featureKey, keyAlpha, keyColor } from "./density.js";
import { HoverOutline, hoveredItem, type ItemOrder } from "./hover.js";
import type { DrawView } from "./view.js";

/** A feature's ink: its color key's, the faint ink where its theme gives none. */
export const featureInk = (key: number): number => (key === NO_COLOR ? SkinToken.inkFaint : keyColor(key));
/** What feature k is drawn in, as `hoverInk` takes it: its color key (alpha folded in), the faint ink without one. */
export const featureHoverKey = (v: DrawView, f: Features, k: number): number => {
  const key = featureKey(v, f, k);
  return key === NO_COLOR ? SkinToken.inkFaint : key;
};

/**
 * The j-th feature of `x` in `hoveredItem`'s span when it is the hovered item's own on this entry (its annotation's, on
 * the entry's row): its index; -1: not.
 */
export function hoveredAt(v: DrawView, f: Features, x: ItemOrder, j: number): number {
  const ann = v.marks.itemAnnotation(v.layer);
  const k = x.order[j];
  return (ann >= 0 && x.annotation[j] !== ann) || (f.on[k] >= 0 && f.on[k] !== v.row) ? -1 : k;
}
const span = new Int32Array(2);

let rows = new Int32Array(64);
/**
 * Per lane its rows [y0, y1) in device px at `2 × lane` (the lane's inset off both sides); `flat`: every lane the whole
 * entry, the tiny path; `lanes`: a streamed layer's (more than one region's table has). Returned scratch, valid until the next call.
 */
export function laneRows(v: DrawView, f: Features, flat: boolean, lanes = f.lanes): Int32Array {
  if (rows.length < 2 * lanes) rows = new Int32Array(4 * lanes);
  const lh = flat ? v.height : laneHeight(v.height, lanes);
  const inset = flat ? 0 : laneInset(lh);
  for (let lane = 0; lane < lanes; lane++) {
    const top = v.top + (flat ? 0 : lanePad(v.height) + lane * lh);
    const y0 = Math.round((top + inset) * v.dpr) | 0;
    rows[2 * lane] = y0;
    rows[2 * lane + 1] = Math.max(y0 + 1, Math.round((top + lh - inset) * v.dpr) | 0);
  }
  return rows;
}

/** A feature layer's density over its entry (2 px off both sides of a row over 8 px): the shared `drawDensity`. */
export function densityBand(v: DrawView, f: Features, p: unknown, point: boolean, rects: RectBatch): void {
  const pad = v.height > 8 ? 2 : 0;
  const y0 = Math.round((v.top + pad) * v.dpr) | 0;
  const y1 = Math.round((v.top + v.height - pad) * v.dpr) | 0;
  drawDensity(v, f, f, point, densityKindOf(p), SkinToken.inkFaint, rects, y0, y1);
}

const line: LineTable = {
  count: 1,
  block: new Int32Array(1),
  start: new Int32Array(1),
  stop: new Int32Array(1),
  top: new Float64Array(1),
};
let used = new Int32Array(32);
const hull = new Int32Array(2);
const win = new Int32Array(2);
let shown = new Int32Array(32);

/**
 * A wrapped line's features as underlines, a lane's height (`wrappedLaneHeight`) apart from the entry's bottom and 2 px
 * less high (3 px of 5), one lane each of those with a feature on
 * the line (`lineLanes`, so the entry's height is `laneLineHeights`'), clamped to the line's positions and cut at the
 * gaps between groups (what is drawn is over positions: it hits back to its item).
 */
export function drawLineLanes(v: DrawView, f: Features, rects: RectBatch): void {
  const b = v.onlyBlock;
  if (!lineFrame(v, f)) return;
  shownHull(v, b, v.lineStart, v.lineStop, hull);
  for (let lane = 0; lane < f.lanes; lane++) {
    if (used[lane] < 0) continue;
    const g = b * f.lanes + lane;
    segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], hull[0], hull[1], win);
    for (let k = win[0]; k < win[1]; k++) {
      const stop = Math.max(f.stop[k], f.start[k] + 1);
      if (stop <= hull[0] || (f.on[k] >= 0 && f.on[k] !== v.row)) continue;
      const key = featureKey(v, f, k);
      rects.begin(featureInk(key), keyAlpha(key));
      lineUnderline(v, f, k, rects);
    }
  }
}
/** The line's lanes in `used` (`lineLanes`); false: no line or no features on its block. */
function lineFrame(v: AxisView, f: Features): boolean {
  const b = v.onlyBlock;
  if (b < 0 || v.lineStop <= v.lineStart || (b + 1) * f.lanes >= f.seg.length) return false;
  line.block[0] = b;
  line.start[0] = v.lineStart;
  line.stop[0] = v.lineStop;
  if (used.length < f.lanes) used = new Int32Array(2 * f.lanes);
  lineLanes(f, v, line, 0, used);
  return true;
}
// Read once: a member read per pointer event allocated under the test runner's module wrapper.
const WRAPPED = AxisGeometry.wrapped;
/**
 * Wrapped: the hit on a line's underlines, as `drawLineLanes` draws them (the hit's `onlyBlock` and line set): lane L's
 * underline a lane's height apart from the entry's bottom, its slot the lane's height, then the feature over the
 * hit's position on that lane, as its item. True: the layer's hit, an item or the position (between lanes, on the
 * band's top padding, or under no feature); false: not on a line here.
 */
export function hitLineLanes(v: AxisView, f: Features, out: Hit): boolean {
  if (v.geometry !== WRAPPED || v.stream !== null || out.block !== v.onlyBlock || out.index < 0) return false;
  if (!lineFrame(v, f)) return false;
  const lane = Math.floor((v.top + v.height - out.point[1] - 1) / v.laneHeight);
  for (let p = 0; lane >= 0 && p < f.lanes; p++) {
    if (used[p] !== lane) continue;
    const k = featureAt(f, out.block, p, out.index);
    if (k >= 0 && (f.on[k] < 0 || f.on[k] === v.row)) hitItem(out, (v.axis as Axis).set, f.annotation[k], f.row[k]);
    break;
  }
  return true;
}
/** Feature k's underlines on the line (its lane's, from `lineFrame`), cut at the gaps between groups, into `rects`. */
function lineUnderline(v: DrawView, f: Features, k: number, rects: RectBatch): void {
  const b = v.onlyBlock;
  const lane = used[f.lane[k]];
  if (lane < 0) return;
  const d = Math.max(1, Math.round(v.dpr)) | 0;
  const lh = v.laneHeight;
  const h = Math.max(1, Math.round((lh > 3 ? lh - 2 : 1) * v.dpr)) | 0;
  const y = Math.round((v.top + v.height - lh - lane * lh) * v.dpr) | 0;
  const stop = Math.max(f.stop[k], f.start[k] + 1);
  let n = v.runToShown(b, f.start[k], stop, shown);
  if (2 * n > shown.length) {
    shown = new Int32Array(4 * n);
    n = v.runToShown(b, f.start[k], stop, shown);
  }
  const run = v.runs[b];
  for (let r = 0; r < n; r++) {
    const end = Math.min(v.lineStop, shown[2 * r + 1]);
    for (let a = Math.max(v.lineStart, shown[2 * r]); a < end; ) {
      const g = groupOf(run, a);
      const e = g + 1 < run.segments ? Math.min(end, run.startShown[g + 1]) : end;
      const x0 = v.xd(b, a);
      rects.push(x0 + d, y, Math.max(2 * d, v.xe(b, e - 1) - x0 - 2 * d), h);
      a = e;
    }
  }
}
const under = new RectBatch(16);
/** Wrapped: the hovered item's underlines on the line, each outlined (`HoverOutline`). */
export function outlineLineLanes(v: DrawView, f: Features): void {
  const x = hoveredItem(v, f, span);
  if (x === null || v.onlyBlock < 0 || !lineFrame(v, f)) return;
  HoverOutline.reset();
  for (let j = span[0]; j < span[1]; j++) {
    const k = hoveredAt(v, f, x, j);
    if (k < 0 || f.block[k] !== v.onlyBlock) continue;
    under.reset();
    lineUnderline(v, f, k, under);
    const key = featureHoverKey(v, f, k);
    for (let i = 0; i < under.n; i++) HoverOutline.box(v, key, under.x[i], under.y[i], under.w[i], under.h[i]);
  }
  HoverOutline.hand(v);
}
