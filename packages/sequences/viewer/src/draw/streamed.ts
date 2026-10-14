/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Streamed items and placeholders (seqspec:sources#levels-and-the-zoom-gate; seqspec:sources#the-viewers-side,
// `StreamView.items`), part of the renderer API: a streamed `ranges` or `variants` layer draws region by region, each
// region's feature table through the renderer's own item draw with the regions' ownership rule (`owned`: an item
// crossing regions is drawn by the first in view that holds it), or past its threshold of items in view, and for a
// level's `count` bins, density per device column (the shared `densityDraw`). A region in view not yet loaded draws a
// quiet placeholder band in the skin's pending shade over its extent (ST3b adds the shimmer); past the letters gate, a
// faint hatched band with a centered hint. Every
// path is allocation-free once warm: module scratch, grown only. ST3b: a streamed item hit (`hitStreamed`: the hover's
// `part` is its number in the layer) and outlined (`outlineStreamed`); skeletons and errors are the overlay's
// (progress.ts). ST7: through a path, a region's bins and an item's range are its system's positions, taken to the
// block's indices (`StreamView.toBlock`).

import { AxisGeometry, type AxisView } from "../axis.ts";
import type { Hit } from "../hit.ts";
import { type Features, featureAt } from "../prepare.ts";
import { segWindow } from "../pyramids.ts";
import { laneHeight, laneInset, lanePad } from "../renderer-helpers.ts";
import { LETTERS_PER_PX, RegionState, type StreamItems } from "../stream-contract.ts";
import { streamStats } from "../stream-sources.ts";
import {
  NO_COLOR,
  type PathBatch,
  PathMode,
  PathVerb,
  type RectBatch,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from "./backend.ts";
import { shownHull, visibleRange } from "./cells.ts";
import {
  densityAdd,
  densityCover,
  densityDraw,
  densityKindOf,
  densityMode,
  densityShown,
  densitySlot,
  densityStart,
  featureKey,
} from "./density.ts";
import { featureHoverKey, laneRows } from "./features.ts";
import { HoverOutline, hoverEmphasis } from "./hover.ts";
import { placeRun } from "./place.ts";
import type { DrawView } from "./view.ts";

/** The region drawn: its items' first positions, and in `from` the index below which they are an earlier region's. */
const owner = { first: null as Int32Array | null };
const from = new Int32Array(1);
/** Feature k of the table drawn is its region's to draw (always, outside a streamed draw). */
export const owned = (k: number): boolean => owner.first === null || owner.first[k] >= from[0];

/**
 * The lanes a streamed layer draws in a band of `height` (CSS px): its own, at most those of `MIN_LANE` px (ST4b: a
 * gene with dozens of transcripts crushed every lane to a pixel); the items of the lanes past it are not drawn.
 */
export const MIN_LANE = 12;
/** The lanes of `MIN_LANE` px a band of `height` holds (at least one): what `StreamView.fold` is told. */
export const fitLanes = (height: number): number =>
  Math.max(1, Math.floor((height - 2 * lanePad(height)) / MIN_LANE)) | 0;
export const laneCap = (height: number, lanes: number): number => Math.min(lanes, fitLanes(height));
/** A region's items drawn by the renderer (`lanes`: the layer's, capped by `laneCap`, for the lanes' rows). */
export type RegionDraw = (v: DrawView, f: Features, p: unknown, lanes: number) => void;

const range = new Int32Array(2);
const hull = new Int32Array(2);
const win = new Int32Array(2);
const out = new Int32Array(7);
/** `toBlock`'s answer. */
const blk = new Int32Array(2);
const piece = new Float64Array(64);

/** Block b's [start, stop) as pieces in device px, a band [y0, y1) in `color` (the placeholder's, the hatch's base). */
export function bandOver(
  v: DrawView,
  b: number,
  start: number,
  stop: number,
  y0: number,
  y1: number,
  color: number,
  rects: RectBatch,
) {
  const n = placeRun(v, b, start, stop, piece);
  if (n === 0) return;
  rects.begin(color);
  for (let q = 0; q < n; q++) {
    const x0 = piece[2 * q] | 0;
    rects.push(x0, y0, Math.max(1, (piece[2 * q + 1] | 0) - x0), y1 - y0);
  }
}

/** Positions a device px of block b, as the stream reads take it (a small integer: a double crosses no call). */
const perPx = (v: DrawView, b: number) => Math.floor(1 / (v.cell * v.span[b] * v.dpr));

/**
 * A streamed items layer (`v.stream.items`): per block and region in view, its features drawn by `each` (`point`:
 * variants, density past 64 in view below two device px a cell; else ranges, past the data area's device columns), or
 * density; pending regions as placeholder bands over the entry. False: the layer is not streamed.
 */
export function drawStreamedItems(
  v: DrawView,
  p: unknown,
  point: boolean,
  tiny: boolean,
  each: RegionDraw,
  rects: RectBatch,
): boolean {
  const st = v.stream;
  const ax = v.axis;
  if (!st?.items || !ax) return false;
  const pad = v.height > 8 ? 2 : 0;
  const y0 = Math.round((v.top + pad) * v.dpr) | 0;
  const y1 = Math.round((v.top + v.height - pad) * v.dpr) | 0;
  // A layer that folds its items (`collapse`) is told whether its lanes fit, before any read of them.
  st.fold?.(fitLanes(v.height));
  // The items in view (each region's window), and whether any region holds bins: items or density, once a frame.
  let count = 0;
  let binned = false;
  for (let b = 0; b < ax.count; b++) {
    if (!view(v, b)) continue;
    out[0] = perPx(v, b);
    for (let k = 0; ; k++) {
      const r = st.items(b, hull[0], hull[1], k, out);
      if (r === null) break;
      if (r.bins !== null) binned = true;
      else if (r.features !== null) count += windowed(r.features, b);
    }
  }
  const columns = Math.max(1, Math.round(v.right * v.dpr) - Math.round(v.left * v.dpr));
  let hint = 0;
  const dense =
    binned ||
    tiny ||
    (v.reductions && (point ? v.cell * v.dpr < 2 && densityMode(v, count, 64) : densityMode(v, count, columns)));
  if (dense) densityStart(v);
  for (let b = 0; b < ax.count; b++) {
    if (!view(v, b)) continue;
    out[0] = perPx(v, b);
    for (let k = 0; ; k++) {
      const r = st.items(b, hull[0], hull[1], k, out);
      if (r === null) break;
      if (out[4] === RegionState.pending) bandOver(v, b, out[1], out[2], y0, y1, SkinToken.pending, rects);
      if (r.features === null && r.bins === null) continue;
      if (r.hint > hint) hint = r.hint;
      owner.first = r.own ?? r.first;
      from[0] = out[3];
      if (dense) densityOf(v, b, r, point);
      else each(v, r.features as Features, p, laneCap(v.height, out[5]));
    }
  }
  owner.first = null;
  if (dense) densityDraw(rects, y0, y1, densityKindOf(p), SkinToken.inkFaint);
  if (hint > 0) hinted(v, y0, y1, hint === 2 ? (point ? GATED_POINTS : GATED) : point ? ZOOM_POINTS : ZOOM);
  else if (!dense && !point && laneCap(v.height, out[5]) < out[5]) hinted(v, y1 - 12 * v.dpr, y1, CAPPED);
  return true;
}

// The zoom gate's hints (ST4b): over estimated counts, and over a gated data region (a click there loads it).
const ZOOM = "zoom in to see features";
const CAPPED = "more lanes than fit: raise the track's height";
const ZOOM_POINTS = "zoom in to see variants";
const GATED = "zoom in to see features, or click to load them";
const GATED_POINTS = "zoom in to see variants, or click to load them";
const hints = new TextBatch(2);
/** The hint centered over the data area between device rows y0 and y1, where it fits. */
function hinted(v: DrawView, y0: number, y1: number, words: string) {
  const dpr = v.dpr;
  const left = Math.round(v.left * dpr);
  const right = Math.round(v.right * dpr);
  const font = v.text.font("ui", Math.max(7, Math.min(11, Math.floor(((y1 - y0) * 0.6) / dpr))), 600);
  if (right - left < v.text.label(words, font).width * dpr + 16 * dpr || y1 - y0 < 9 * dpr) return;
  hints.reset();
  hints.style(font, TextAlign.center, TextBaseline.middle);
  hints.begin(SkinToken.inkQuiet);
  const t = hints.add(words);
  hints.x[t] = (left + right) / 2;
  hints.y[t] = Math.round((y0 + y1) / 2) + 0.5;
  v.backend.text(hints);
}

/** Block b's range drawn and its hull into `range`, `hull`; false: none here. */
function view(v: DrawView, b: number): boolean {
  visibleRange(v, b, range);
  if (range[0] >= range[1]) return false;
  shownHull(v, b, range[0], range[1], hull);
  return true;
}
/** Features of block b reaching the hull (every lane's window). */
function windowed(f: Features, b: number): number {
  let n = 0;
  for (let lane = 0; lane < f.lanes && (b + 1) * f.lanes < f.seg.length; lane++) {
    const g = b * f.lanes + lane;
    segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], hull[0], hull[1], win);
    n += win[1] - win[0];
  }
  return n;
}

/** A region into the density columns: its bins at their start's column, or its owned items (a point at its center). */
function densityOf(v: DrawView, b: number, r: StreamItems, point: boolean) {
  const bins = r.bins;
  const st = v.stream;
  if (bins !== null) {
    const s = densitySlot(NO_COLOR);
    for (let j = 0; j < bins.length; j++) {
      const a = out[6] + j * r.width;
      // Its start on the block (through a path, the first index it reaches; a view without `toBlock`: itself).
      if (bins[j] <= 0 || (st?.toBlock && !st.toBlock(b, a, a + r.width, blk))) continue;
      const at = st?.toBlock ? blk[0] : a;
      if (at < hull[0] || at >= hull[1]) continue;
      const q = v.shownOf(b, at);
      if (q >= range[0] && q < range[1]) densityAdd(v.xd(b, q), s, bins[j] | 0);
    }
    return;
  }
  const f = r.features as Features;
  for (let lane = 0; lane < f.lanes && (b + 1) * f.lanes < f.seg.length; lane++) {
    const g = b * f.lanes + lane;
    segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], hull[0], hull[1], win);
    for (let k = win[0]; k < win[1]; k++) {
      if (!owned(k)) continue;
      const s = densitySlot(featureKey(v, f, k));
      if (point) {
        const q = v.shownOf(b, (f.start[k] + Math.max(f.stop[k], f.start[k] + 1) - 1) >> 1);
        if (q >= range[0] && q < range[1]) densityAdd(v.xd(b, q), s, 1);
        continue;
      }
      const n = placeRun(v, b, f.start[k], Math.max(f.stop[k], f.start[k] + 1), piece);
      for (let q = 0; q < n; q++)
        densityCover(piece[2 * q] | 0, Math.max(piece[2 * q] | 0, (piece[2 * q + 1] | 0) - 1), s);
    }
  }
}

const hio = new Int32Array(7);
/**
 * A streamed items layer's item under the hit (`v.stream`): in the region holding its index, by the lane under y
 * (`lanes`: ranges') or the first at the index (variants); its number in `out.part`, what the hover keeps (its label,
 * its outline). False: the layer is not streamed here.
 */
export function hitStreamed(v: AxisView, out: Hit, lanes: boolean): boolean {
  const st = v.stream;
  if (!st?.items || out.block < 0 || v.geometry !== AxisGeometry.linear || v.row !== -1) return false;
  hio[0] = Math.floor(1 / (v.cell * v.span[out.block] * v.dpr));
  // The regions over the index in order: through a path a split codon's parts may lie in two (ST7).
  let r: StreamItems | null = null;
  let k = -1;
  let q = 0;
  for (; k < 0; q++) {
    const x = st.items(out.block, out.index, out.index + 1, q, hio);
    if (x === null) break;
    const f = x.features;
    if (f === null || x.item === null) continue;
    r = x;
    if (lanes) {
      const cap = laneCap(v.height, hio[5]);
      const h = laneHeight(v.height, cap);
      const at = (out.point[1] - v.top - lanePad(v.height)) / h;
      const lane = Math.floor(at);
      if (lane >= cap) return true;
      if (at - lane >= laneInset(h) / h && at - lane < 1 - laneInset(h) / h)
        k = featureAt(f, out.block, lane, out.index);
    } else for (let lane = 0; lane < f.lanes && k < 0; lane++) k = featureAt(f, out.block, lane, out.index);
  }
  // None over it: not streamed here; else the hit is the layer's, an item or not.
  if (r === null) return q > 0;
  if (k >= 0 && r.ends !== null && r.item !== null) {
    out.part = r.item[k];
    // Its range on the block: through a path, the hull of what it reaches.
    const ok = st.toBlock?.(out.block, r.first[k], Math.max(r.ends[k], r.first[k] + 1), blk) ?? false;
    out.span[0] = ok ? blk[0] : r.first[k];
    out.span[1] = ok ? blk[1] : Math.max(r.ends[k], r.first[k] + 1);
  }
  return true;
}

const el = new Int32Array(4);
/**
 * The hovered streamed item (the layer's own hit, its `part`) outlined over its features in the regions in view, in
 * its lane (`lanes`: ranges') or the entry's band; nothing in density. False: the layer is not streamed here.
 */
export function outlineStreamed(v: DrawView, lanes: boolean): boolean {
  const st = v.stream;
  const ax = v.axis;
  if (!st?.items || !ax) return false;
  if (!hoverEmphasis(v) || !v.marks.element(v.layer, el) || el[3] < 0 || (v.reductions && densityShown(v))) return true;
  HoverOutline.reset();
  for (let b = 0; b < ax.count; b++) {
    if (!view(v, b)) continue;
    out[0] = perPx(v, b);
    for (let k = 0; ; k++) {
      const r = st.items(b, hull[0], hull[1], k, out);
      if (r === null) break;
      const f = r.features;
      if (f === null || r.item === null) continue;
      const cap = laneCap(v.height, out[5]);
      const ys = laneRows(v, f, !lanes, cap);
      for (let q = 0; q < f.count; q++) {
        if (r.item[q] !== el[3] || (lanes && f.lane[q] >= cap)) continue;
        const y0 = ys[2 * (lanes ? f.lane[q] : 0)];
        const y1 = ys[2 * (lanes ? f.lane[q] : 0) + 1];
        const n = placeRun(v, b, f.start[q], Math.max(f.stop[q], f.start[q] + 1), piece);
        for (let j = 0; j < n; j++) {
          const x0 = piece[2 * j] | 0;
          HoverOutline.box(v, featureHoverKey(v, f, q), x0, y0, Math.max(1, (piece[2 * j + 1] | 0) - x0), y1 - y0);
        }
      }
    }
  }
  HoverOutline.hand(v);
  return true;
}

/**
 * Block b's letters gate (`StreamView.gate`): false where a device px holds at most that many positions (letters are
 * asked); else, over shown [lo, hi), a faint hatched band between device rows y0 and y1 and the hint "zoom in to see
 * letters" centered where it fits.
 */
export function lettersGate(
  v: DrawView,
  b: number,
  lo: number,
  hi: number,
  y0: number,
  y1: number,
  paths: PathBatch,
  text: TextBatch,
): boolean {
  const gate = v.stream?.gate ?? LETTERS_PER_PX;
  if (v.cell * v.span[b] * v.dpr * gate >= 1) return false;
  const dpr = v.dpr;
  const left = Math.round(v.left * dpr);
  const right = Math.round(v.right * dpr);
  const x0 = Math.max(left, v.xd(b, lo));
  const x1 = Math.min(right, v.xd(b, hi));
  const h = y1 - y0;
  if (x1 <= x0 || h <= 0) return true;
  // 45° strokes every 6 px, clipped to the band.
  const step = 6 * dpr;
  paths.begin(SkinToken.inkFaint, 110, PathMode.stroke, Math.max(1, dpr));
  for (let x = x0 - h - ((x0 - h) % step); x < x1; x += step) {
    const a = Math.max(x, x0);
    const e = Math.min(x + h, x1);
    if (e <= a) continue;
    const o = paths.verb(PathVerb.move);
    paths.coords[o] = a;
    paths.coords[o + 1] = y1 - (a - x);
    const q = paths.verb(PathVerb.line);
    paths.coords[q] = e;
    paths.coords[q + 1] = y1 - (e - x);
  }
  const font = v.text.font("ui", Math.max(7, Math.min(11, Math.floor((h * 0.6) / dpr))), 400);
  const w = v.text.label(HINT, font).width * dpr;
  if (x1 - x0 < w + 16 * dpr || h < 9 * dpr) return true;
  text.style(font, TextAlign.center, TextBaseline.middle);
  text.begin(SkinToken.inkQuiet);
  const t = text.add(HINT);
  text.x[t] = (x0 + x1) / 2;
  text.y[t] = Math.round((y0 + y1) / 2) + 0.5;
  return true;
}
const HINT = "zoom in to see letters";

let runs = new Int32Array(64);
const part = new Float64Array(5);
/**
 * A streamed column's statistics into `out` (as `streamStats`), shown [k0, k1) read run by run: across a cyclic origin
 * or a filter's gap the stream is asked each run's positions, never everything between its two ends (VT-06's draw
 * side). Its minimum, maximum, sum and count merged; `out[0]` and `out[1]` NaN where no run has a value.
 */
export function streamColumn(v: DrawView, b: number, k0: number, k1: number, out: Float64Array): void {
  let n = v.shownToRuns(b, k0, k1, runs);
  if (n <= 1) {
    streamStats(v, b, k0, k1, out);
    return;
  }
  if (2 * n > runs.length) {
    runs = new Int32Array(4 * n);
    n = v.shownToRuns(b, k0, k1, runs);
  }
  out[0] = out[1] = Number.NaN;
  out[2] = out[3] = 0;
  for (let r = 0, k = k0; r < n; r++) {
    const e = k + runs[2 * r + 1] - runs[2 * r];
    streamStats(v, b, k, e, part);
    k = e;
    if (part[3] === 0) continue;
    if (!(part[0] >= out[0])) out[0] = part[0];
    if (!(part[1] <= out[1])) out[1] = part[1];
    out[2] += part[2];
    out[3] += part[3];
  }
}
