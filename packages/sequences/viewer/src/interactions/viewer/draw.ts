/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fingerprint's draw and `tiny` on batches (seqspec:viewer#extension-points; seqspec:libraries#interactions, as
// KLIFS shows fingerprints), in device pixels. A frame builds the groups of contacted cells in view once (adjacent
// cells in one device column merge: their kinds are the OR of theirs), from prepare's contacted lists: the walk is
// bounded by the contacted positions, not by the positions in view, whatever the zoom. Then per row and group: a gray
// ground and a stripe per kind held, or, in a compact row, one rect in the color of the first kind. Grounds go in a
// batch before the stripes, so no backend's color grouping can draw one over a stripe. Blocks without the property are
// drawn by the `matrix` renderer's tiles. The overlay outlines the hovered stripe, else the hovered cell, inside its
// own pixels (seqspec:interaction#hover).

import {
  AxisGeometry,
  type DrawView,
  drawMatrixTiles,
  hoverEmphasis,
  type MatrixState,
  NO_COLOR,
  outlineItem,
  outlineMatrixCell,
  packAlpha,
  RectBatch,
  SkinToken,
  type Theme,
  ThemeLocation,
  visibleRange,
} from "@molstar/seq-viewer/draw";
import { compactRows, Geo, geo, LEAD, rowFrame, rowTop, stripeTop } from "./geometry.ts";
import type { FingerprintCache, FingerprintState } from "./renderer.ts";

const grounds = new RectBatch(256);
const stripes = new RectBatch(1024);
const range = new Int32Array(2);
/** The entry's top and height around the fallback tiles (a typed slot: a double kept across a call may be boxed). */
const saved = new Float64Array(2);
/** Per stripe, its kind's color and alpha in the layer's theme. */
const stripeColor = new Int32Array(16);
const stripeAlpha = new Uint8Array(16);
/** The frame's groups: device px edges, the range in `gCol` (fingerprint columns) and, per block, its first group. */
let gx0 = new Int32Array(256);
let gx1 = new Int32Array(256);
let gFirst = new Int32Array(256);
let gCount = new Int32Array(256);
let gCol = new Int32Array(1024);
let blockGroups = new Int32Array(16);

const grown = (a: Int32Array<ArrayBuffer>): Int32Array<ArrayBuffer> => {
  const b = new Int32Array(2 * a.length);
  b.set(a);
  return b;
};
function growGroups(): void {
  gx0 = grown(gx0);
  gx1 = grown(gx1);
  gFirst = grown(gFirst);
  gCount = grown(gCount);
}
function growColumns(): void {
  gCol = grown(gCol);
}

/** The kinds' colors through the layer's color theme (a kind without one: the quiet ink), into the stripe tables. */
function stripeColors(v: DrawView, s: FingerprintState): void {
  const t = v.color;
  const alpha = t === null ? undefined : t.alpha;
  const loc = v.loc;
  ThemeLocation.reset(loc);
  loc.layer = v.layer;
  for (let j = 0; j < s.kinds.length; j++) {
    loc.category = s.kinds[j];
    const c = t === null ? NO_COLOR : t.value(loc);
    stripeAlpha[j] = c >= 0 && alpha !== undefined ? alpha(loc) : 255;
    stripeColor[j] = c === NO_COLOR ? SkinToken.inkQuiet : c;
  }
  loc.category = -1;
}

/**
 * Into the module's tables: the contacted cells in view, per fingerprint block, grouped by device column. Returns the
 * group count; `blockGroups[b]..[b + 1]` are block b's.
 */
function groups(v: DrawView, s: FingerprintState): number {
  const ax = v.axis;
  const nb = s.fp.length;
  if (ax === null) return 0;
  if (blockGroups.length < nb + 1) blockGroups = new Int32Array(2 * nb + 2);
  let n = 0;
  let c = 0;
  for (let b = 0; b < nb; b++) {
    blockGroups[b] = n;
    if (s.fp[b] === null) continue;
    visibleRange(v, b, range);
    const lo = range[0];
    const hi = range[1];
    if (lo >= hi) continue;
    const base = ax.base[b];
    const first = n;
    const idx = s.contacted[b];
    const cols = s.contactedCol[b];
    for (let i = 0; i < idx.length; i++) {
      const k = v.shownOf(b, base + idx[i]);
      if (k < lo || k >= hi) continue;
      const x0 = v.xd(b, k);
      const x1 = v.xe(b, k);
      if (c === gCol.length) growColumns();
      gCol[c] = cols[i];
      if (n > first && gx0[n - 1] === x0) {
        if (x1 > gx1[n - 1]) gx1[n - 1] = x1;
        gCount[n - 1]++;
      } else {
        if (n === gx0.length) growGroups();
        gx0[n] = x0;
        gx1[n] = x1;
        gFirst[n] = c;
        gCount[n] = 1;
        n++;
      }
      c++;
    }
  }
  blockGroups[nb] = n;
  return n;
}

/** The blocks without the property as the `matrix` renderer's tiles, in the counts channel, inside `inset` px. */
function tiles(v: DrawView, s: FingerprintState, inset: number): void {
  const held = v.cache as FingerprintCache | null;
  const color = v.color;
  saved[0] = v.top;
  saved[1] = v.height;
  try {
    v.cache = held === null ? null : held.tiles;
    v.color = (v.themes[1] ?? null) as Theme<"color"> | null;
    v.top += inset;
    v.height -= 2 * inset;
    drawMatrixTiles(v, s.fallback as MatrixState);
  } finally {
    v.top = saved[0];
    v.height = saved[1];
    v.cache = held;
    v.color = color;
  }
}

export function draw(v: DrawView, s: FingerprintState, _p: unknown): void {
  const rows = s.matrix.states;
  if (rows === 0 || v.axis === null || v.geometry !== AxisGeometry.linear) return;
  rowFrame(v, rows);
  const a = geo[Geo.a];
  const h = geo[Geo.height];
  const gap = geo[Geo.gap];
  if (h < 1) return;
  if (s.fallback !== null) tiles(v, s, 1);
  const K = s.kinds.length;
  if (K === 0 || groups(v, s) === 0) return;
  stripeColors(v, s);
  const compact = compactRows(K);
  grounds.reset();
  grounds.begin(SkinToken.grid);
  stripes.reset();
  const nb = s.fp.length;
  for (let q = 0; q < rows; q++) {
    const y0 = rowTop(a, h, gap, rows, q);
    const y1 = rowTop(a, h, gap, rows, q + 1) - gap;
    if (y1 <= y0) continue;
    const hh = y1 - y0;
    for (let b = 0; b < nb; b++) {
      if (q >= s.rowsOf[b]) continue;
      const cells = s.cells[b];
      const C = s.width[b];
      for (let g = blockGroups[b], g1 = blockGroups[b + 1]; g < g1; g++) {
        let m = 0;
        for (let i = gFirst[g], e = i + gCount[g]; i < e; i++) m |= cells[q * C + gCol[i]];
        if (m === 0) continue;
        const x0 = gx0[g];
        const w = gx1[g] - x0;
        if (compact) {
          const j = s.stripeOf[LEAD[m]];
          if (j < 0) continue;
          stripes.begin(stripeColor[j], stripeAlpha[j]);
          stripes.push(x0, y0, w < 1 ? 1 : w, hh);
          continue;
        }
        const inset = w >= 4 ? 1 : 0;
        const gw = w - 2 * inset < 1 ? 1 : w - 2 * inset;
        grounds.push(x0 + inset, y0, gw, hh);
        for (let x = m; x !== 0; x &= x - 1) {
          const j = s.stripeOf[31 - Math.clz32(x & -x)];
          if (j < 0) continue;
          const t0 = stripeTop(y0, hh, K, j);
          stripes.begin(stripeColor[j], stripeAlpha[j]);
          stripes.push(x0 + inset, t0, gw, stripeTop(y0, hh, K, j + 1) - t0);
        }
      }
    }
  }
  if (grounds.n > 0) v.backend.rects(grounds);
  if (stripes.n > 0) v.backend.rects(stripes);
}

/** Below the minimum height: one compact row over the whole entry, from the kinds any row holds; no ground. */
export function tiny(v: DrawView, s: FingerprintState, _p: unknown): void {
  if (s.matrix.states === 0 || v.axis === null || v.geometry !== AxisGeometry.linear) return;
  if (s.fallback !== null) tiles(v, s, 0);
  if (s.kinds.length === 0 || groups(v, s) === 0) return;
  stripeColors(v, s);
  const y0 = Math.round(v.top * v.dpr) | 0;
  const y1 = Math.round((v.top + v.height) * v.dpr) | 0;
  const hh = y1 <= y0 ? 1 : y1 - y0;
  stripes.reset();
  for (let b = 0; b < s.fp.length; b++) {
    const union = s.union[b];
    for (let g = blockGroups[b], g1 = blockGroups[b + 1]; g < g1; g++) {
      let m = 0;
      for (let i = gFirst[g], e = i + gCount[g]; i < e; i++) m |= union[gCol[i]];
      const j = m === 0 ? -1 : s.stripeOf[LEAD[m]];
      if (j < 0) continue;
      const w = gx1[g] - gx0[g];
      stripes.begin(stripeColor[j], stripeAlpha[j]);
      stripes.push(gx0[g], y0, w < 1 ? 1 : w, hh);
    }
  }
  if (stripes.n > 0) v.backend.rects(stripes);
}

const el = new Int32Array(4);
/** A stripe's color key, its alpha folded in (as `featureKey`'s), for `hoverInk`. */
const stripeKey = (j: number) => packAlpha(stripeColor[j], stripeAlpha[j]);
/**
 * The hovered cell's emphasis (the layer's own hit, `marks.element`): its stripe (`part`, a kind it holds) outlined
 * inside, else the cell's ground (a compact row: its rect, in its first kind's color); a block without the property, its
 * matrix cell as `matrix` outlines it. Nothing over a cell without contacts.
 */
export function overlay(v: DrawView, s: FingerprintState, _p: unknown): void {
  const rows = s.matrix.states;
  const ax = v.axis;
  if (!hoverEmphasis(v) || rows === 0 || ax === null || v.geometry !== AxisGeometry.linear) return;
  if (v.height < 4 || !v.marks.element(v.layer, el)) return;
  const b = el[0];
  const q = el[2];
  if (b < 0 || b >= s.fp.length || q < 0) return;
  if (s.fp[b] === null) {
    if (s.fallback === null || s.fallback.rows[b] === null) return;
    const color = v.color;
    saved[0] = v.top;
    saved[1] = v.height;
    try {
      v.color = (v.themes[1] ?? null) as Theme<"color"> | null;
      v.top += 1;
      v.height -= 2;
      outlineMatrixCell(v, s.fallback);
    } finally {
      v.top = saved[0];
      v.height = saved[1];
      v.color = color;
    }
    return;
  }
  const j = el[1] - ax.base[b];
  const cols = s.colOf[b];
  const col = q < s.rowsOf[b] && j >= 0 && j < cols.length ? cols[j] : -1;
  const m = col >= 0 ? s.cells[b][q * s.width[b] + col] : 0;
  const k = m !== 0 ? v.shownOf(b, el[1]) : -1;
  if (k < 0) return;
  rowFrame(v, rows);
  const a = geo[Geo.a];
  const h = geo[Geo.height];
  const gap = geo[Geo.gap];
  const y0 = rowTop(a, h, gap, rows, q);
  const hh = rowTop(a, h, gap, rows, q + 1) - gap - y0;
  if (h < 1 || hh <= 0) return;
  stripeColors(v, s);
  const x0 = v.xd(b, k);
  const w = Math.max(1, v.xe(b, k) - x0);
  const K = s.kinds.length;
  if (compactRows(K)) {
    const t = s.stripeOf[LEAD[m]];
    if (t >= 0) outlineItem(v, stripeKey(t), x0, y0, w, hh);
    return;
  }
  const inset = w >= 4 ? 1 : 0;
  const gw = w - 2 * inset < 1 ? 1 : w - 2 * inset;
  const part = el[3];
  const t = part >= 0 && m & (1 << part) ? s.stripeOf[part] : -1;
  if (t < 0) outlineItem(v, SkinToken.grid, x0 + inset, y0, gw, hh);
  else {
    const t0 = stripeTop(y0, hh, K, t);
    outlineItem(v, stripeKey(t), x0 + inset, t0, gw, stripeTop(y0, hh, K, t + 1) - t0);
  }
}
