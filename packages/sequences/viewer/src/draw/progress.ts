/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Progressive loading on the overlay (seqspec:sources#levels-and-the-zoom-gate; seqspec:sources#the-viewers-side;
// ST3b), the composer's: a streamed layer's regions in view at the view's level, drawn over the main pass's cells,
// which are never redrawn for it. A pending region with nothing drawn in its place is a skeleton: the skin's pending
// shade and, while motion is on, a slow shimmer (a soft band of the surface sweeping the data area once every
// `SWEEP_MS`). A failed one is an error band: the shade tinted in the conflict color, hatched, "failed to load"
// centered where it fits. One first drawn less than `FADE_MS` ago is the pending shade fading out over its data (its
// data fading in; regions first drawn in one frame fade together). Lazy rows take the same per tile, a coarser tile
// drawn in a pending one's place counting as drawn. `progressEnd` says whether another overlay frame is needed (a
// skeleton shimmering, a region fading). Allocation-free: module scratch.

import { type RowModel, SectionLayouts } from "../layout.js";
import { blockLevels, RegionState, type StreamView } from "../streams.js";
import { PathBatch, PathMode, PathVerb, RectBatch, SkinToken, TextAlign, TextBaseline, TextBatch } from "./backend.js";
import { shownHull, visibleRange } from "./cells.js";
import { placeRun } from "./place.js";
import type { DrawView } from "./view.js";

/** A shimmer's sweep across the data area, and a region's fade-in (ms). */
export const SWEEP_MS = 1000;
export const FADE_MS = 150;
/** The shimmer: half its band (CSS px), its strips a side, the surface's alpha at its middle. */
const HALF = 90;
const STEPS = 6;
const PEAK = 120;
const TINT = 40;
const HINT = "failed to load";

const fills = new RectBatch(64);
const sheen = new RectBatch(64);
const hatch = new PathBatch(64);
const words = new TextBatch(4);
/** Reads: the frame's time in, then a region's [start, stop), state and age. */
const io = new Float64Array(6);
/** The frame's time, the sweep's center and half width (device px); motion on (1); another frame needed (1). */
const at = new Float64Array(3);
const flags = new Int32Array(2);
const range = new Int32Array(2);
const hull = new Int32Array(2);
const piece = new Float64Array(64);
const lv = new Int32Array(4);
let runs = new Int32Array(16);

/**
 * A frame's start: its time (ms, `now[0]`: a double passed to a call is boxed), whether motion is on (a live surface,
 * `Viewport.ReducedMotion` not set) and whether skeletons shimmer (`Stream.ShimmerFps` above 0; ST4b).
 */
export function progressBegin(v: DrawView, now: Float64Array, motion: boolean, shimmers: boolean): void {
  at[0] = motion ? now[0] : -1;
  const left = Math.round(v.left * v.dpr);
  const right = Math.round(v.right * v.dpr);
  at[2] = HALF * v.dpr;
  at[1] = left - at[2] + ((now[0] % SWEEP_MS) / SWEEP_MS) * (right - left + 2 * at[2]);
  flags[0] = (motion ? 1 : 0) | (motion && shimmers ? 2 : 0);
  flags[1] = 0;
}
/** What moves, so another overlay frame is needed: `Moving.shimmer`, `Moving.fade` (bits; 0: nothing). */
export const Moving = { shimmer: 1, fade: 2 } as const;
/** The frame's end: what moved (`Moving` bits). */
export function progressEnd(): number {
  return flags[1];
}

/** A 1-D streamed layer's regions in view (`StreamView.region`) over the entry's band (2 px off a band over 8 px). */
export function progress1D(v: DrawView, st: StreamView): void {
  const ax = v.axis;
  if (!ax || !st.region) return;
  const pad = v.height > 8 ? 2 : 0;
  const y0 = Math.round((v.top + pad) * v.dpr) | 0;
  const y1 = Math.round((v.top + v.height - pad) * v.dpr) | 0;
  for (let b = 0; b < ax.count; b++) {
    visibleRange(v, b, range);
    if (range[0] >= range[1]) continue;
    shownHull(v, b, range[0], range[1], hull);
    io[0] = at[0];
    for (let k = 0; st.region(b, hull[0], hull[1], k, io); k++) {
      const s = io[3] | 0;
      if ((s === RegionState.empty || s === RegionState.loaded) && !(io[4] < FADE_MS)) continue;
      const n = placeRun(v, b, io[1] | 0, io[2] | 0, piece);
      for (let q = 0; q < n; q++) mark(v, piece[2 * q] | 0, piece[2 * q + 1] | 0, y0, y1, s);
    }
  }
}

/**
 * A lazy rows band's tiles in view (`StreamView.tileAt`), at the view's levels (`blockLevels`), by runs of display rows
 * on consecutive lanes of one tile row: a pending tile whose coarser tile is loaded counts as drawn (its fade, the
 * coarser one's).
 */
export function progressRows(v: DrawView, st: StreamView, r: RowModel): void {
  const ax = v.axis;
  const ln = v.lanes;
  if (!ax || !st.tileAt || ln[1] <= ln[0]) return;
  for (let b = 0; b < ax.count; b++) {
    const lo = v.visLo[b];
    const hi = v.visHi[b];
    if (hi <= lo || (v.onlyBlock >= 0 && b !== v.onlyBlock)) continue;
    blockLevels(v, b, st.budget, ln, lv, st.served);
    let n = v.shownToRuns(b, lo, hi, runs);
    if (2 * n > runs.length) {
      runs = new Int32Array(4 * n);
      n = v.shownToRuns(b, lo, hi, runs);
    }
    const shift = lv[0] + 8;
    for (let d = ln[2] | 0; d <= ln[3]; ) {
      const k = SectionLayouts.laneOf(r, d);
      let e = d + 1;
      while (k >= 0 && e <= ln[3] && SectionLayouts.laneOf(r, e) === k + (e - d) && (k + e - d) >> shift === k >> shift)
        e++;
      if (k >= 0) {
        const ya = Math.round((ln[4] + d * ln[5]) * v.dpr) | 0;
        const yb = Math.round((ln[4] + e * ln[5]) * v.dpr) | 0;
        for (let q = 0, ks = lo; q < n; q++) {
          const c0 = runs[2 * q];
          const c1 = runs[2 * q + 1];
          // A reversed run (its first column at its end) is left out: its tiles are not contiguous in x.
          if (c1 - c0 === 1 || v.indexOf(b, ks as never) === c0) tiles(v, st, b, k, ks, c0, c1, ya, yb);
          ks += c1 - c0;
        }
      }
      d = e;
    }
  }
}

/** Columns [c0, c1) of lane k (shown from `ks`) by tile of the view's level, between device rows ya and yb. */
function tiles(
  v: DrawView,
  st: StreamView,
  b: number,
  k: number,
  ks: number,
  c0: number,
  c1: number,
  ya: number,
  yb: number,
) {
  const view = st as Required<StreamView>;
  const ty = k >> (lv[0] + 8);
  for (let tx = c0 >> (lv[1] + 8); tx <= (c1 - 1) >> (lv[1] + 8); tx++) {
    const a = Math.max(c0, tx << (lv[1] + 8));
    const z = Math.min(c1, (tx + 1) << (lv[1] + 8));
    io[0] = at[0];
    let s = view.tileAt(lv[0], lv[1], ty, tx, io);
    const behind = s === RegionState.pending || s === RegionState.failed;
    const coarse = behind ? view.tileAt(lv[2], lv[3], k >> (lv[2] + 8), a >> (lv[3] + 8), io) : -1;
    // The coarser tile drawn in its place: a pending one is drawn, a failed one only marked (`io[5]`).
    io[5] = coarse === RegionState.loaded ? 1 : 0;
    if (coarse === RegionState.loaded && s === RegionState.pending) s = coarse;
    io[4] = io[1];
    if ((s === RegionState.empty || s === RegionState.loaded) && !(io[4] < FADE_MS)) continue;
    mark(v, v.xd(b, ks + (a - c0)), v.xd(b, ks + (z - c0)), ya, yb, s);
  }
}

/** One region's pieces [x0, x1) × [y0, y1) (device px) by its state; its age in `io[4]`. */
function mark(v: DrawView, x0: number, x1: number, y0: number, y1: number, s: number) {
  if (x1 <= x0 || y1 <= y0) return;
  if (s === RegionState.pending) {
    fills.begin(SkinToken.pending);
    fills.push(x0, y0, x1 - x0, y1 - y0);
    if ((flags[0] & 2) !== 0) shimmer(x0, x1, y0, y1);
  } else if (s === RegionState.failed && io[5] === 1) {
    // A coarser level drawn in its place (ST4b): no band over the data, a thin mark along the band's foot.
    const h = Math.min(y1 - y0, Math.max(2, Math.round(2 * v.dpr))) | 0;
    fills.begin(SkinToken.conflict);
    fills.push(x0, y1 - h, x1 - x0, h);
  } else if (s === RegionState.failed) {
    fills.begin(SkinToken.conflict, TINT);
    fills.push(x0, y0, x1 - x0, y1 - y0);
    hatched(v, x0, x1, y0, y1);
  } else if ((flags[0] & 1) !== 0) {
    // Loaded or empty less than `FADE_MS` ago: the shade fading out over it.
    fills.begin(SkinToken.pending, Math.round(255 * (1 - io[4] / FADE_MS)) | 0);
    fills.push(x0, y0, x1 - x0, y1 - y0);
    flags[1] |= Moving.fade;
  }
}

/** The sweep's strips over [x0, x1): the surface, strongest at the band's middle. */
function shimmer(x0: number, x1: number, y0: number, y1: number) {
  flags[1] |= Moving.shimmer;
  const w = at[2] / STEPS;
  for (let q = 0; q < 2 * STEPS; q++) {
    const a = Math.max(x0, Math.round(at[1] - at[2] + q * w)) | 0;
    const z = Math.min(x1, Math.round(at[1] - at[2] + (q + 1) * w)) | 0;
    if (z <= a) continue;
    sheen.begin(SkinToken.surface, Math.round(PEAK * (1 - Math.abs(q + 0.5 - STEPS) / STEPS)) | 0);
    sheen.push(a, y0, z - a, y1 - y0);
  }
}

/** 45° strokes every 6 px in the conflict color over the band, and the hint centered where it fits. */
function hatched(v: DrawView, x0: number, x1: number, y0: number, y1: number) {
  const h = y1 - y0;
  const step = 6 * v.dpr;
  hatch.begin(SkinToken.conflict, 150, PathMode.stroke, Math.max(1, v.dpr));
  for (let x = x0 - h - ((x0 - h) % step); x < x1; x += step) {
    const a = Math.max(x, x0);
    const e = Math.min(x + h, x1);
    if (e <= a) continue;
    const o = hatch.verb(PathVerb.move);
    hatch.coords[o] = a;
    hatch.coords[o + 1] = y1 - (a - x);
    const p = hatch.verb(PathVerb.line);
    hatch.coords[p] = e;
    hatch.coords[p + 1] = y1 - (e - x);
  }
  const font = v.text.font("ui", Math.max(7, Math.min(11, Math.floor((h * 0.6) / v.dpr))), 600);
  if (x1 - x0 < v.text.label(HINT, font).width * v.dpr + 16 * v.dpr || h < 9 * v.dpr) return;
  words.style(font, TextAlign.center, TextBaseline.middle);
  words.begin(SkinToken.conflict);
  const t = words.add(HINT);
  words.x[t] = (x0 + x1) / 2;
  words.y[t] = Math.round((y0 + y1) / 2) + 0.5;
}

/** What was gathered for an entry, drawn under its clip: shades, the sweep, hatches, hints. */
export function progressFlush(v: DrawView): void {
  if (fills.n > 0) v.backend.rects(fills);
  if (sheen.n > 0) v.backend.rects(sheen);
  if (hatch.n > 0) v.backend.paths(hatch);
  if (words.n > 0) v.backend.text(words);
  fills.reset();
  sheen.reset();
  hatch.reset();
  words.reset();
}
