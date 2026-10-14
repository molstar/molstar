/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `points` and `variants` draws (seqspec:viewer#extension-points; seqspec:sources#levels-and-the-zoom-gate), part
// of the renderer API: a glyph at each item's position in view (variants stacked per position by their place in it, an
// alternative's letter on a mark that can hold it), wrapped underlines; thinner than two device px a cell with more
// than 64 items in view, or on the tiny path, density per column and color, each item once at its center (the shared
// `drawDensity`). Every item counts. The overlay outlines the hovered item's glyph or its wrapped underline inside its
// own cell; in density, nothing (seqspec:interaction#hover).

import { AxisGeometry, type AxisView, type ShownIndex } from "../axis.js";
import type { Hit } from "../hit.js";
import type { Features } from "../prepare.js";
import { segWindow } from "../pyramids.js";
import { hitItem } from "../renderer-helpers.js";
import { NO_SIZE, SIZE_ONE, ThemeLocation } from "../themes.js";
import { CHARS, GlyphBatch, RectBatch, SkinToken, TextAlign, TextBaseline, TextBatch } from "./backend.js";
import { shownHull, visibleRange } from "./cells.js";
import { densityMode, densityShown, featureGlyph, featureKey, itemsInView, keyAlpha } from "./density.js";
import { densityBand, drawLineLanes, featureHoverKey, featureInk, hoveredAt, outlineLineLanes } from "./features.js";
import { HoverOutline, hoverEmphasis, hoveredItem } from "./hover.js";
import { boundaryEdge, placeRun } from "./place.js";
import {
  DOT_R,
  DOT_R_MAX,
  DOT_STACK_R,
  DOT_STACK_R_MAX,
  featureSize,
  hitLayer,
  hitReach,
  hitWindowAt,
  sizeTheme,
} from "./sizes.js";
import { drawStreamedItems, owned } from "./streamed.js";
import type { DrawView } from "./view.js";

const bars = new RectBatch(256);
const dots = new GlyphBatch(256);
const alts = new TextBatch(32);
const range = new Int32Array(2);
const hull = new Int32Array(2);
const win = new Int32Array(2);
const at = new Float64Array(2);
const ctr = new Float64Array(1);
const runs = new Int32Array(8);
const piece = new Float64Array(8);

/** "1" to "99": an insertion's length on its variant's marker. */
const COUNTS = Array.from({ length: 100 }, (_, n) => String(n));

// The dot's radii read once: a module member read per frame allocates under the test runner's module wrapper.
const LEAST = DOT_R;
const MOST = DOT_R_MAX;
const STACK_LEAST = DOT_STACK_R;
const STACK_MOST = DOT_STACK_R_MAX;

/** The glyphs' radius (CSS px), the stack's pitch (device px), its rows and the first row's center (device px). */
const frame = new Float64Array(4);
function glyphFrame(v: AxisView, stacked: boolean) {
  // The viewer's dot, grown with the cell up to a most, held by its band's height; never under 2.5 px.
  const least = stacked ? STACK_LEAST : LEAST;
  const most = stacked ? STACK_MOST : MOST;
  let r = v.cell * 0.45;
  if (r > most) r = most;
  if (r < least) r = least;
  const fit = (v.height - 2) / (stacked ? 2.2 : 3.5);
  if (r > fit) r = fit;
  if (r < 2.5) r = 2.5;
  frame[0] = r;
  frame[1] = (2 * r + 1) * v.dpr;
  frame[2] = stacked ? Math.max(1, Math.floor((v.height - 4) / (2 * r + 1))) : 1;
  frame[3] = Math.round((stacked ? v.top + 2 + r : v.top + v.height / 2) * v.dpr);
}

/**
 * Block b's item [start, stop)'s glyph center into `out[0]` (CSS px): the middle of the longest piece of it in view (a
 * shown run clipped to the view: one run, reversed or filtered alike; two across a cyclic origin, the first on a tie),
 * so a residue over its codon is centered on it, and a deletion partly scrolled out keeps its marker on what is shown.
 * False: none of it in view. The draw and the hit share it (VL-05).
 */
function itemCenter(v: AxisView, b: number, start: number, stop: number, out: Float64Array): boolean {
  const n = v.runToShown(b, start, stop, runs);
  let k0 = 0;
  let k1 = 0;
  for (let j = 0; j < n && j < 4; j++) {
    const a = Math.max(runs[2 * j], v.visLo[b]);
    const e = Math.min(runs[2 * j + 1], v.visHi[b]);
    if (e - a > k1 - k0) {
      k0 = a;
      k1 = e;
    }
  }
  if (k1 <= k0) return false;
  out[0] = (v.x(b, k0 as ShownIndex) + v.x(b, k1 as ShownIndex)) / 2;
  return true;
}

/**
 * Item k's glyph center into `at` (device px): the middle of its positions' cells in view (`itemCenter`); a boundary (an insertion) on the edge before its position in reading order, or after the piece's
 * last position (seqspec:format#placement). False: not shown here.
 */
function glyphAt(v: DrawView, b: number, f: Features, k: number): boolean {
  const s = f.start[k];
  if (v.onlyBlock >= 0 && b !== v.onlyBlock) return false;
  if (f.stop[k] > s) {
    if (!itemCenter(v, b, s, f.stop[k], ctr)) return false;
    at[0] = ctr[0] * v.dpr;
    at[1] = 0;
    return true;
  }
  const e = boundaryEdge(v, b, s);
  if (e < 0 || e < v.visLo[b] || e > v.visHi[b]) return false;
  at[0] = v.xd(b, e as never);
  at[1] = 0;
  return true;
}

/** Glyphs at the items in view: `stacked` (variants) a row each by their place in their position, in a mark of 6 px. */
function glyphsAt(v: DrawView, f: Features, stacked: boolean) {
  const ax = v.axis;
  if (!ax) return;
  const dpr = v.dpr;
  glyphFrame(v, stacked);
  const r = frame[0];
  const pitch = frame[1];
  const rows = frame[2];
  const cy = frame[3];
  const letters = stacked && r >= 4.5;
  // The radius channel (points): a glyph's size in px, up to what its row holds; the thickness channel (variants): a
  // share of its place in the stack.
  const sizes = sizeTheme(v, 0, 1);
  const most = Math.max(r, (v.height - 2) / 2);
  // An alternative's letter: a capital's ink centered on its mark (seqspec:viewer#the-render-backend, text placement).
  let lift = 0;
  if (letters) {
    const font = v.text.font("ui", Math.max(7, Math.floor(r * 1.4)), 600);
    alts.style(font, TextAlign.center, TextBaseline.alphabetic);
    lift = v.text.ink(font, Math.round(dpr * 64))[0] / 2;
  }
  for (let b = 0; b < ax.count && (b + 1) * f.lanes < f.seg.length; b++) {
    visibleRange(v, b, range);
    if (range[0] >= range[1]) continue;
    shownHull(v, b, range[0], range[1], hull);
    for (let lane = 0; lane < f.lanes; lane++) {
      const g = b * f.lanes + lane;
      // One past the hull: a boundary after the last position in view (at the piece's end) is drawn on its edge.
      segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], hull[0], hull[1] + 1, win);
      for (let k = win[0]; k < win[1]; k++) {
        if (
          (f.stop[k] <= f.start[k] ? f.start[k] < hull[0] : f.stop[k] <= hull[0] || f.start[k] >= hull[1]) ||
          (f.on[k] >= 0 && f.on[k] !== v.row) ||
          !owned(k)
        )
          continue;
        const shown = glyphAt(v, b, f, k);
        const wide = stacked && f.stop[k] - f.start[k] > 1;
        if (!shown && !wide) continue;
        const key = featureKey(v, f, k);
        const y = cy + Math.min(f.kind[k], rows - 1) * pitch;
        if (wide) {
          // A variant over several positions (a deletion): its one marker on a thin line across their cells, a piece
          // per shown run in view (reversed, across a cyclic origin, partly scrolled out alike).
          const n = placeRun(v, b, f.start[k], f.stop[k], piece);
          if (n > 0) bars.begin(featureInk(key), keyAlpha(key));
          for (let j = 0; j < n; j++)
            bars.push(
              piece[2 * j],
              Math.round(y - dpr / 2),
              Math.max(1, piece[2 * j + 1] - piece[2 * j]),
              Math.max(1, dpr),
            );
          if (!shown) continue;
        }
        const id = featureGlyph(v, f, k, 0);
        dots.begin(featureInk(key), keyAlpha(key));
        const i = dots.add(id < 0 ? 0 : id);
        dots.x[i] = at[0];
        dots.y[i] = y;
        const size = featureSize(v, f, k, sizes);
        const share = Math.max(1, (r * Math.min(size, SIZE_ONE)) / SIZE_ONE);
        dots.radius[i] = (size < 0 ? r : stacked ? share : Math.min(most, Math.max(1, size / SIZE_ONE))) * dpr;
        // An alternative's letter, or an insertion's length (a negative code) where two digits fit.
        const code = f.code[k];
        if (!letters || code === 0 || code < -99 || (code < -9 && r < 5.5) || (size >= 0 && size < SIZE_ONE)) continue;
        alts.begin(SkinToken.surface);
        const t = alts.add(code > 0 ? CHARS[code & 127] : COUNTS[-code]);
        alts.x[t] = at[0];
        alts.y[t] = Math.round(dots.y[i] + lift);
      }
    }
  }
}

/**
 * Points (`stacked` false) or variants into the pass's backend: glyphs, or density past 64 items in view below two
 * device px a cell (`tiny`: density always), or in a wrapped line underlines on the lanes of the line.
 */
export function drawPoints(v: DrawView, f: Features, p: unknown, stacked: boolean, tiny: boolean): void {
  bars.reset();
  dots.reset();
  alts.reset();
  ThemeLocation.reset(v.loc);
  // Streamed (ST3a): region by region, or density, pending regions as placeholders.
  if (v.geometry === AxisGeometry.linear && drawStreamedItems(v, p, true, tiny, stacked ? variant : point, bars)) {
  } else if (v.geometry === AxisGeometry.wrapped) drawLineLanes(v, f, bars);
  else if (tiny || (v.reductions && v.cell * v.dpr < 2 && densityMode(v, itemsInView(v, f, 65), 64)))
    densityBand(v, f, p, true, bars);
  else glyphsAt(v, f, stacked);
  if (bars.n > 0) v.backend.rects(bars);
  if (dots.n > 0) {
    dots.font = v.text.font("ui", 12, 600);
    v.backend.glyphs(dots);
  }
  if (alts.n > 0) v.backend.text(alts);
}

const variant = (v: DrawView, f: Features) => glyphsAt(v, f, true);
const point = (v: DrawView, f: Features) => glyphsAt(v, f, false);

const near = new Int32Array(2);
/**
 * Points' and variants' hit by drawn geometry: the item whose glyph holds the point, the nearest center where glyphs
 * overlap, as `glyphsAt` places and sizes them (`stacked`: a variant at its place in its position's stack, a boundary
 * on its edge); false: none (the column then decides). Linear, held items.
 */
export function hitPoints(v: AxisView, f: Features, out: Hit, stacked = false): boolean {
  const ax = v.axis;
  if (!ax || v.geometry !== AxisGeometry.linear || v.stream !== null) return false;
  glyphFrame(v, stacked);
  const r = frame[0];
  const most = stacked ? r : Math.max(r, (v.height - 2) / 2);
  const px = out.point[0];
  const b = out.block;
  if (b < 0 || (b + 1) * f.lanes >= f.seg.length) return false;
  const dy0 = out.point[1] - (stacked ? frame[3] / v.dpr : v.top + v.height / 2);
  if (dy0 < -most || (!stacked && dy0 > most)) return false;
  // Through module scratch: the point and reach (doubles) would be boxed across the call, every move.
  hitReach[0] = px;
  hitReach[1] = most;
  if (!hitWindowAt(v, b, near)) return false;
  const pitch = frame[1] / v.dpr;
  const sizes = sizeTheme(hitLayer, 0, 1);
  let best = -1;
  let least = Number.POSITIVE_INFINITY;
  for (let lane = 0; lane < f.lanes; lane++) {
    const g = b * f.lanes + lane;
    // One past the window: a boundary after the block's last position.
    segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], near[0], near[1] + 1, win);
    for (let k = win[0]; k < win[1]; k++) {
      if (f.on[k] >= 0 && f.on[k] !== v.row) continue;
      let cx: number;
      if (f.stop[k] <= f.start[k]) {
        // A boundary, on its edge as `glyphAt` places it.
        const e = boundaryEdge(v, b, f.start[k]);
        if (e < 0) continue;
        cx = v.x(b, e as never);
      } else {
        if (!itemCenter(v, b, f.start[k], f.stop[k], ctr)) continue;
        cx = ctr[0];
      }
      const size = stacked ? NO_SIZE : featureSize(hitLayer, f, k, sizes);
      const rk = size < 0 ? r : Math.min(most, Math.max(1, size / SIZE_ONE));
      const dy = stacked ? dy0 - Math.min(f.kind[k], frame[2] - 1) * pitch : dy0;
      const d = (px - cx) * (px - cx) + dy * dy;
      if (d <= (rk + 0.5) * (rk + 0.5) && d < least) {
        least = d;
        best = k;
      }
    }
  }
  return best >= 0 && hitItem(out, ax.set, f.annotation[best], f.row[best]);
}

const hovered = new Int32Array(2);
/**
 * The hovered item: where `drawPoints` draws glyphs (a variant's at its place in the stack), a ring around its glyph
 * inside its cell, else the glyph outlined inset; wrapped, its underline outlined; nothing in density or on the tiny
 * path.
 */
export function outlinePoints(v: DrawView, f: Features, stacked: boolean, tiny: boolean): void {
  if (!hoverEmphasis(v) || !v.axis) return;
  if (v.geometry === AxisGeometry.wrapped) {
    outlineLineLanes(v, f);
    return;
  }
  const x = v.geometry === AxisGeometry.linear && !tiny ? hoveredItem(v, f, hovered) : null;
  if (x === null || (v.reductions && v.cell * v.dpr < 2 && densityShown(v))) return;
  glyphFrame(v, stacked);
  const dpr = v.dpr;
  // Whole device px, inside the drawn glyph (its center is rounded here, not there: half a px under the radius).
  const r = (frame[0] * dpr - 0.5) | 0;
  // The glyph's cell: its row of the stack (the entry's band unstacked) by its positions' cells.
  const tall = (stacked ? frame[1] / 2 : (v.height * dpr) / 2) | 0;
  HoverOutline.reset();
  for (let j = hovered[0]; j < hovered[1]; j++) {
    const k = hoveredAt(v, f, x, j);
    const b = k < 0 ? 0 : f.block[k];
    if (k < 0 || !glyphAt(v, b, f, k)) continue;
    const cy = frame[3] + Math.min(f.kind[k], frame[2] - 1) * frame[1];
    const wide = ((v.cell * v.span[b] * dpr * Math.max(1, f.stop[k] - f.start[k])) / 2) | 0;
    const id = featureGlyph(v, f, k, 0);
    const key = featureHoverKey(v, f, k);
    // The glyph's own radius under a size theme, as the draw places it.
    const size = featureSize(v, f, k, sizeTheme(v, 0, 1));
    const most = Math.max(frame[0], (v.height - 2) / 2);
    const share = Math.max(1, (frame[0] * Math.min(size, SIZE_ONE)) / SIZE_ONE);
    const rk = size === NO_SIZE ? r : ((stacked ? share : Math.min(most, Math.max(1, size / SIZE_ONE))) * dpr) | 0;
    HoverOutline.glyph(
      v,
      key,
      id < 0 ? 0 : id,
      Math.round(at[0]) | 0,
      Math.round(cy) | 0,
      rk,
      Math.min(Math.max(tall, rk), wide),
    );
  }
  HoverOutline.hand(v);
}
