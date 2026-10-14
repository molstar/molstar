/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The matrix's tiles (seqspec:viewer#extension-points `matrix`; seqspec:sources#levels-and-the-zoom-gate), part of the
// renderer API: per block and level (samples of 2^lx positions, their per-position maximum, at least a device px each),
// tiles of at most 4,096 samples, built through `images.fromPixels` when first in view and kept per prepared state,
// theme, backend and tile height. A frame draws one `image` per tile of each shown run (`shownToRuns`: an unfiltered
// block is one run), never a scan of positions; a reversed piece draws mirrored tiles. The bands are the drawn states
// in their order (`stateOf`, seqspec:format#keyed-sets). The entry's device frame is `matrixFrame`: its top rounded,
// its height ⌊height·dpr⌋ (from the height alone, so a scroll by a fraction of a pixel keeps the tiles). A tile has a
// row per band, which the backend stretches over the frame (device row y shows band ⌊(2y+1)·B/2H⌋, `matrixBandOfRow`),
// or, where bands are thinner than a device pixel (KS2), a row per device pixel: band q in row ⌊q·H/B⌋ (ints, so the
// rows are whole device pixels at any DPR, as the interactions fingerprint's), by the layer's `density`: `majority`,
// the row's hue the plurality of its bands' colors (ties: the color met first), opaque by the share of its bands that
// have a color, as gaps fade an alignment's tiles; `maximum`, the band with the largest value (the first on a tie) in its
// own color, as samples of positions take theirs. The hit reduces the one row again from the data, as the tile did
// (`reduceRow`: at most ⌈B/H⌉ bands × 2^lx positions), and names the band that drew it (the first with the plurality's
// color, or the strongest): the state whose color is drawn (`matrixHitBand`), whichever cache drew last (VL-04). The
// overlay outlines the hovered cell (`outlineMatrixCell`; seqspec:interaction#hover).

import { AxisGeometry, type AxisView } from '../axis.js';
import type { MatrixState } from '../rows.js';
import { type Theme, ThemeLocation } from '../themes.js';
import { type BackendImages, NO_COLOR } from './backend.js';
import { visibleRange } from './cells.js';
import { packAlpha } from './density.js';
import { hoverEmphasis, outlineItem } from './hover.js';
import { hitLayer, type SizeView } from './sizes.js';
import type { DrawView } from './view.js';

/** Samples a tile holds at most (a canvas past about 32k px fails). */
export const MATRIX_TILE = 4096;
/** What a frame did (K4's counters): runs drawn, tiles built, images drawn, positions scanned outside a build. */
export const MatrixVisits = { runs: 0, tiles: 0, images: 0, scans: 0 };

/** One backend's tiles: per (block, level, mirrored) the image ids by tile, -1 until built. */
export interface MatrixTileSet {
  images: BackendImages;
  ids: (Int32Array | undefined)[];
}
/** Per prepared state (the engine's `cache`, else the module's): the theme the tiles were made in and a set per backend. */
export class MatrixTiles {
  theme: unknown = null;
  /** The tiles' rows: the bands, or the entry's device pixels where those are fewer (a reduced tile). */
  height = 0;
  /** The prepared state drawn. */
  state: MatrixState | null = null;
  readonly sets: MatrixTileSet[] = [];
  dispose(): void {
    for (const t of this.sets)
      for (const ids of t.ids) for (let i = 0; ids && i < ids.length; i++) if (ids[i] >= 0) t.images.dispose(ids[i]);
    this.sets.length = 0;
  }
}
/** A draw with no engine cache keeps its tiles here, by prepared state. */
const own = new WeakMap<MatrixState, MatrixTiles>();
/** Edge-rate, apart from the draw. */
function newTiles(s: MatrixState): MatrixTiles {
  const c = new MatrixTiles();
  own.set(s, c);
  return c;
}
function newSet(c: MatrixTiles, images: BackendImages): MatrixTileSet {
  const t = { images, ids: [] };
  c.sets.push(t);
  return t;
}
function tiles(v: DrawView, s: MatrixState, height: number): MatrixTiles {
  const c = v.cache instanceof MatrixTiles ? v.cache : (own.get(s) ?? newTiles(s));
  const fresh = c.state !== s;
  if (fresh) c.state = s;
  // Another prepared state (the data, the order or the `density`), theme or height: the tiles are made again.
  if (fresh || c.theme !== v.color || c.height !== height) {
    c.dispose();
    c.theme = v.color;
    c.height = height;
  }
  return c;
}
function tileSet(c: MatrixTiles, images: BackendImages): MatrixTileSet {
  for (let i = 0; i < c.sets.length; i++) if (c.sets[i].images === images) return c.sets[i];
  return newSet(c, images);
}

/**
 * The entry's device frame into `out`: its top (rounded) and its height, ⌊height·dpr⌋ (from the height alone: a frame
 * scrolled by a fraction keeps it). Ints out; the draw, the hit and the outline share it.
 */
export function matrixFrame(v: AxisView, out: Int32Array): void {
  out[0] = Math.round(v.top * v.dpr) | 0;
  out[1] = (v.height * v.dpr) | 0;
}
/** The band a stretched tile (B rows over H device rows, sampled at pixel centers) shows in device row y. */
export function matrixBandOfRow(y: number, rows: number, bands: number): number {
  return (((2 * y + 1) * bands) / (2 * rows)) | 0;
}
/** The first device row that shows band q of a stretched tile: the least y with ⌊(2y+1)·B/2H⌋ ≥ q. */
export function matrixRowOfBand(q: number, rows: number, bands: number): number {
  const n = 2 * q * rows - bands;
  return n <= 0 ? 0 : ((n + 2 * bands - 1) / (2 * bands)) | 0;
}

const LEVELS = 32;
let px = new Uint32Array(MATRIX_TILE * 4);
/** A reduced row's bands' colors and their votes (edge rate, grown only). */
let votes = new Int32Array(64);
/** The first band of reduced row y of H rows over B bands: the least q with ⌊q·H/B⌋ ≥ y. Ints in and out. */
export function matrixFirstBand(y: number, rows: number, bands: number): number {
  return ((y * bands + rows - 1) / rows) | 0;
}
/**
 * Block b's tile t at level lx (`mirror`: samples right to left), built on first use (edge rate): per sample and row
 * the maximum over its positions through the color theme, as RGBA pixels, transparent where no position has a value;
 * a reduced row (`set`'s height under the bands), the plurality of its bands' colors.
 */
function tile(
  v: DrawView,
  s: MatrixState,
  set: MatrixTileSet,
  b: number,
  lx: number,
  t: number,
  mirror: number,
  h: number,
) {
  const key = (b * LEVELS + lx) * 2 + mirror;
  const samples = Math.ceil((s.rows[b] as Int32Array).length / 2 ** lx);
  let ids = set.ids[key];
  if (!ids) {
    ids = new Int32Array(Math.ceil(samples / MATRIX_TILE)).fill(-1);
    set.ids[key] = ids;
  }
  if (ids[t] < 0) ids[t] = build(v, s, set.images, b, lx, t, mirror, samples, h);
  return ids[t];
}
/** The quiet gray drawn without a theme, as a pixel (an int32: the pixel's bits). */
const QUIET = (255 << 24) | (0x98 << 16) | (0x8f << 8) | 0x8a;
/** `cellMax`'s result (a typed slot: a double returned is boxed). */
const cell = new Float64Array(1);
/** Band q's largest value at block b's positions [j0, j1) into `cell[0]`; NaN where none has one. */
function cellMax(s: MatrixState, b: number, q: number, j0: number, j1: number): void {
  const rows = s.rows[b] as Int32Array;
  const values = s.ann[b]?.matrix?.values;
  const stride = s.stride[b];
  const st = s.stateOf[q];
  let max = Number.NaN;
  for (let j = j0; values && st < stride && j < j1; j++) {
    const val = rows[j] >= 0 ? values[rows[j] * stride + st] : Number.NaN;
    if (!Number.isNaN(val) && !(val <= max)) max = val;
  }
  cell[0] = max;
}
/**
 * Band q's color at block b's positions [j0, j1) as a pixel's bits in an int32 (a double past 2^31 would be boxed
 * across the call): `theme` at their maximum; 0 where none has a value or the theme gives no color. The draw's view or
 * the hit's layer (`hitLayer`).
 */
function cellColor(
  v: SizeView,
  theme: Theme<'color'> | null,
  s: MatrixState,
  b: number,
  q: number,
  j0: number,
  j1: number,
  base: number,
): number {
  cellMax(s, b, q, j0, j1);
  if (Number.isNaN(cell[0])) return 0;
  const ax = v.axis;
  if (!theme || !ax) return QUIET;
  const loc = ThemeLocation.at(v.loc, ax, v.layer, b, base + j0);
  loc.value = cell[0];
  loc.state = s.stateOf[q];
  const c = theme.value(loc);
  const a = c >= 0 && theme.alpha ? theme.alpha(loc) : 255;
  loc.value = Number.NaN;
  loc.state = -1;
  // ImageData's byte order, R first in memory.
  return c < 0 ? 0 : (a << 24) | ((c & 255) << 16) | (c & 0xff00) | ((c >> 16) & 255);
}
/** The band of rows [q0, q1) with the largest value at [j0, j1) (the first on a tie); -1: none has one. */
function strongest(s: MatrixState, b: number, q0: number, q1: number, j0: number, j1: number): number {
  let best = -1;
  let max = Number.NEGATIVE_INFINITY;
  for (let q = q0; q < q1; q++) {
    cellMax(s, b, q, j0, j1);
    if (!Number.isNaN(cell[0]) && (best < 0 || cell[0] > max)) {
      max = cell[0];
      best = q;
    }
  }
  return best;
}
/** `reduceRow`'s pixel (an int32: its bits) and band offset. */
const reduced = new Int32Array(2);
/**
 * Reduced row [q0, q1) of block b's positions [j0, j1): the pixel into `reduced[0]` (0: none) and the offset in the row
 * of the band that drew it into `reduced[1]` (-1: none). `maximum`: the strongest band in its own color; `majority`: the
 * plurality of the bands' colors (the color met first on a tie), opaque by the share of the bands with a color, drawn
 * by its first band of that color. Shared by the tile and the hit, so the hit names what the tile drew.
 */
function reduceRow(
  v: SizeView,
  theme: Theme<'color'> | null,
  s: MatrixState,
  b: number,
  q0: number,
  q1: number,
  j0: number,
  j1: number,
  base: number,
): void {
  reduced[0] = 0;
  reduced[1] = -1;
  if (s.maximum) {
    const q = strongest(s, b, q0, q1, j0, j1);
    const pixel = q < 0 ? 0 : cellColor(v, theme, s, b, q, j0, j1, base);
    if (pixel === 0) return;
    reduced[0] = pixel;
    reduced[1] = q - q0;
    return;
  }
  if (votes.length < 2 * (q1 - q0 + 1)) votes = new Int32Array(4 * (q1 - q0 + 1));
  let n = 0;
  let colored = 0;
  for (let q = q0; q < q1; q++) {
    const pixel = cellColor(v, theme, s, b, q, j0, j1, base);
    if (pixel === 0) continue;
    colored++;
    let k = 0;
    while (k < n && votes[2 * k] !== pixel) k++;
    if (k === n) {
      votes[2 * n] = pixel;
      votes[2 * n + 1] = 0;
      n++;
    }
    votes[2 * k + 1]++;
  }
  if (n === 0) return;
  let best = 0;
  for (let k = 1; k < n; k++) if (votes[2 * k + 1] > votes[2 * best + 1]) best = k;
  const m = q1 - q0;
  const hue = votes[2 * best];
  const alpha = (((hue >>> 24) * colored + (m >> 1)) / m) | 0;
  reduced[0] = (hue & 0xffffff) | (alpha << 24);
  for (let q = q0; q < q1; q++)
    if (cellColor(v, theme, s, b, q, j0, j1, base) === hue) {
      reduced[1] = q - q0;
      return;
    }
}
function build(
  v: DrawView,
  s: MatrixState,
  images: BackendImages,
  b: number,
  lx: number,
  t: number,
  mirror: number,
  samples: number,
  H: number,
) {
  const rows = s.rows[b] as Int32Array;
  const B = s.states;
  const W = Math.min(MATRIX_TILE, samples - t * MATRIX_TILE);
  if (px.length < W * H) px = new Uint32Array(W * H);
  px.fill(0, 0, W * H);
  const base = v.axis ? v.axis.base[b] : 0;
  const theme = v.color;
  MatrixVisits.tiles++;
  for (let x = 0; s.ann[b]?.matrix && x < W; x++) {
    const j0 = (t * MATRIX_TILE + x) * 2 ** lx;
    const j1 = Math.min(rows.length, j0 + 2 ** lx);
    const col = mirror ? W - 1 - x : x;
    for (let y = 0; y < H; y++) {
      if (H >= B) {
        const pixel = cellColor(v, theme, s, b, y, j0, j1, base);
        if (pixel !== 0) px[y * W + col] = pixel;
        continue;
      }
      reduceRow(v, theme, s, b, matrixFirstBand(y, H, B), matrixFirstBand(y + 1, H, B), j0, j1, base);
      if (reduced[0] !== 0) px[y * W + col] = reduced[0];
    }
  }
  return images.fromPixels(W, H, px);
}

/**
 * The band the hit names at block b's block index i and device row y of the frame (H rows): a stretched tile's band
 * there; a reduced row's band whose color the tile draws, reduced again from the data at the level the frame draws
 * (`reduceRow`, with the layer's themes from `hitLayer`); its first band where nothing is drawn. -1: outside. Ints in
 * and out.
 */
export function matrixHitBand(v: AxisView, s: MatrixState, b: number, i: number, y: number, H: number): number {
  const B = s.states;
  if (y < 0 || y >= H || B === 0) return -1;
  if (H >= B) return matrixBandOfRow(y, H, B);
  const q0 = matrixFirstBand(y, H, B);
  const rows = s.rows[b];
  if (!rows || !s.ann[b]?.matrix || i < 0 || i >= rows.length) return q0;
  const lx = levelOf(v, b);
  const j0 = ((i / 2 ** lx) | 0) * 2 ** lx;
  const j1 = Math.min(rows.length, j0 + 2 ** lx);
  const base = v.axis ? v.axis.base[b] : 0;
  reduceRow(hitLayer, colorOf(hitLayer.themes), s, b, q0, matrixFirstBand(y + 1, H, B), j0, j1, base);
  return reduced[1] < 0 ? q0 : q0 + reduced[1];
}
/** A layer's first color theme, as the draw's `color` (`DrawView.setThemes`). */
function colorOf(themes: readonly (Theme | null)[]): Theme<'color'> | null {
  for (let k = 0; k < themes.length; k++) if (themes[k]?.kind === 'color') return themes[k] as Theme<'color'>;
  return null;
}
/** The level of block b's device columns: samples of 2^lx positions, each at least a device pixel. */
function levelOf(v: AxisView, b: number): number {
  const per = 1 / (v.cell * v.span[b] * v.dpr);
  return per <= 1 ? 0 : Math.min(LEVELS - 1, Math.ceil(Math.log2(per)));
}

const r = new Float64Array(8);
const frame = new Int32Array(2);
const range = new Int32Array(2);
let runs = new Int32Array(32);
/**
 * Block b's run [i0, i1) of block indices, shown over device x [x0, x1): each tile of level lx it meets, one image
 * (fractional source edges, so the run's positions land on exactly its pixels).
 */
function drawRun(
  v: DrawView,
  s: MatrixState,
  set: MatrixTileSet,
  b: number,
  lx: number,
  i0: number,
  i1: number,
  x0: number,
  x1: number,
  mirror: number,
  h: number,
) {
  const base = v.axis ? v.axis.base[b] : 0;
  const u0 = (i0 - base) / 2 ** lx;
  const u1 = (i1 - base) / 2 ** lx;
  const scale = (x1 - x0) / (u1 - u0);
  const samples = Math.ceil((s.rows[b] as Int32Array).length / 2 ** lx);
  MatrixVisits.runs++;
  for (let t = Math.floor(u0 / MATRIX_TILE); t * MATRIX_TILE < u1; t++) {
    const t0 = t * MATRIX_TILE;
    const W = Math.min(MATRIX_TILE, samples - t0);
    const a = Math.max(u0, t0);
    const e = Math.min(u1, t0 + W);
    if (e <= a) continue;
    r[0] = mirror ? t0 + W - e : a - t0;
    r[2] = e - a;
    r[4] = mirror ? x0 + (u1 - e) * scale : x0 + (a - u0) * scale;
    r[6] = (e - a) * scale;
    MatrixVisits.images++;
    v.backend.image(tile(v, s, set, b, lx, t, mirror, h), r, false);
  }
}

/**
 * Each block's shown runs in view as tiles at the level of its device columns, bands top to bottom over the entry: a
 * row each, or reduced to a row per device pixel where they are more.
 */
export function drawMatrixTiles(v: DrawView, s: MatrixState): void {
  const ax = v.axis;
  if (!ax || v.geometry !== AxisGeometry.linear || s.states === 0) return;
  matrixFrame(v, frame);
  const dh = frame[1];
  if (dh <= 0) return;
  const h = dh < s.states ? dh : s.states;
  const c = tiles(v, s, h);
  const set = tileSet(c, v.backend.images);
  r[1] = 0;
  r[3] = h;
  r[5] = frame[0];
  r[7] = dh;
  for (let b = 0; b < s.rows.length; b++) {
    visibleRange(v, b, range);
    if (!s.rows[b] || range[0] >= range[1]) continue;
    const lx = levelOf(v, b);
    let n = ax.shownToRuns(b, range[0], range[1], runs);
    if (2 * n > runs.length) {
      runs = new Int32Array(4 * n);
      n = ax.shownToRuns(b, range[0], range[1], runs);
    }
    // Runs in shown order, each of consecutive block indices: a run's shown length is its length.
    for (let j = 0, k = range[0]; j < n; j++) {
      const i0 = runs[2 * j];
      const i1 = runs[2 * j + 1];
      const mirror = i1 - i0 > 1 && ax.indexOf(b, k as never) !== i0 ? 1 : 0;
      drawRun(v, s, set, b, lx, i0, i1, v.xd(b, k), v.xd(b, k + i1 - i0), mirror, h);
      k += i1 - i0;
    }
  }
}

const el = new Int32Array(4);
/**
 * The hovered cell outlined: the position's columns (a device px at least) by the state's band as the tiles place it
 * (a reduced tile: its device row), in the hover ink of its color (the theme at the cell's value: a reduced row's hit
 * names the state of the drawn color). Only for the layer's own hit (`marks.element`).
 */
export function outlineMatrixCell(v: DrawView, s: MatrixState): void {
  const ax = v.axis;
  if (!hoverEmphasis(v) || !ax || v.geometry !== AxisGeometry.linear || !v.marks.element(v.layer, el)) return;
  const b = el[0];
  const q = el[2];
  const band = q >= 0 && q < s.bandOf.length ? s.bandOf[q] : -1;
  const rows = b < s.rows.length ? s.rows[b] : null;
  const k = rows !== null && band >= 0 ? v.shownOf(b, el[1]) : -1;
  if (k < 0 || rows === null) return;
  const i = el[1] - ax.base[b];
  const row = i >= 0 && i < rows.length ? rows[i] : -1;
  const values = s.ann[b]?.matrix?.values;
  const stride = s.stride[b];
  const value = row >= 0 && values && q < stride ? values[row * stride + q] : Number.NaN;
  let key = NO_COLOR;
  const theme = v.color;
  if (!Number.isNaN(value) && theme) {
    const loc = ThemeLocation.at(v.loc, ax, v.layer, b, el[1]);
    loc.value = value;
    loc.state = q;
    const c = theme.value(loc);
    const a = c >= 0 && theme.alpha ? theme.alpha(loc) : 255;
    loc.value = Number.NaN;
    loc.state = -1;
    key = c < 0 ? c : packAlpha(c, a);
  } else if (!Number.isNaN(value)) key = 0x8a8f98;
  const x0 = v.xd(b, k);
  const x1 = Math.max(x0 + 1, v.xe(b, k));
  matrixBandRows(v, s, band, frame);
  outlineItem(v, key, x0, frame[0], x1 - x0, Math.max(1, frame[1] - frame[0]));
}
/**
 * Band q's device rows [out[0], out[1]) in the entry's frame: a reduced tile's band is its device row
 * (`matrixFirstBand`'s inverse); else the rows its stretched band shows in.
 */
export function matrixBandRows(v: AxisView, s: MatrixState, q: number, out: Int32Array): void {
  matrixFrame(v, out);
  const y = out[0];
  const h = out[1];
  const B = s.states;
  out[0] = y + (h < B ? ((q * h) / B) | 0 : matrixRowOfBand(q, h, B));
  out[1] = h < B ? out[0] + 1 : y + matrixRowOfBand(q + 1, h, B);
}
