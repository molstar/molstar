/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// (its `imageOf`, `paint` and `cells`; the rest of the rows draw is renderers/rows.ts')
// The `rows` layer's tile images (seqspec:viewer#extension-points; seqspec:viewer#the-render-backend,
// seqspec:sources#levels-and-the-zoom-gate): 256 × 256 images, a sample a pixel, drawn magnified onto the rows through
// `images.fromPixels`, never smoothed (a sample is at least a device pixel). An alignment in memory at the finest level
// is themed per cell (row, column, letter); a reduced tile (the region cache's, `DrawView.stream`) per sample's
// majority residue, opaque by its share of the sample's cells (gaps fade), always with the theme's alpha (case 3's
// conservation shading). Each tile drawn at its whole scale, in parts; the coarse level in the place of a fine tile not
// arrived, and there alone. Images are
// made on a miss (edge rate: `builds` counts them) and kept per (level, tile) in an LRU of slots, all dropped when the
// theme or the backend's images change; a frame whose images are made allocates nothing.

import { TILE } from "@molstar/seq-core/format";
import {
  type Axis,
  type BackendImages,
  type DrawView,
  NO_COLOR,
  OPAQUE,
  RectBatch,
  type RowsState,
  rowAt,
  type Theme,
  ThemeLocation,
} from "../draw/api.js";

/** Levels per axis (the region cache's); the overview's images (`ov`) are a set of levels of their own. */
const L = 16;
const OVERVIEW = L * L;
const SLOTS = 320;
/** A reduced sample of gaps only: a quiet gray (an all-gap region reads as such, never as nothing). */
const GAPS = 0x9aa0a8;
/** The least opacity of a reduced sample on the overview's map, out of 255 (a coarse majority is small: 45%). */
const OVERVIEW_SHARE = 115;

/** Per prepared state and skin (the rows entry's `cache`): tile images by level and tile, in slots with LRU stamps. */
export class RowTiles {
  /** Images made, ever: a frame with its images made adds none. */
  builds = 0;
  theme: Theme | null = null;
  images: BackendImages | null = null;
  clock = 0;
  /** Per level (ly·16 + lx, the overview's past 256): tile index → slot. */
  readonly levels: (Map<number, number> | undefined)[] = [];
  /** Per slot: its image (-1: free), level, tile index and the frame it was last drawn in. */
  img = new Int32Array(SLOTS).fill(-1);
  level = new Int32Array(SLOTS);
  index = new Int32Array(SLOTS);
  stamp = new Int32Array(SLOTS);
  used = 0;
  /** Every image disposed (a theme or backend change, the cache's end). */
  dispose(): void {
    for (let k = 0; k < this.used; k++) if (this.img[k] >= 0) this.images?.dispose(this.img[k]);
    this.img.fill(-1);
    this.used = 0;
    this.levels.length = 0;
  }
}

/** A draw's start: the images dropped when the theme or the backend's images changed; a new frame for the LRU. */
export function beginTiles(v: DrawView, c: RowTiles): void {
  if (c.theme !== v.color || c.images !== v.backend.images) {
    c.dispose();
    c.theme = v.color;
    c.images = v.backend.images;
  }
  c.clock++;
}

/** One tile's pixels, shared by every image made (ImageData's byte order: R, G, B, A in memory). */
const px = new Uint32Array(TILE * TILE);
/** 0xRRGGBB and an alpha byte as a little-endian RGBA pixel, an int32 (a Smi-sized return; the array stores it). */
const pixel = (c: number, a: number) => (a << 24) | ((c & 255) << 16) | (c & 0xff00) | ((c >> 16) & 255);

/**
 * Tile (ty, tx) of level (ly, lx) as an image, block b's columns: kept, or made from the held codes (finest level) or
 * the region cache's tile; -1: not arrived. `ov`: the overview's (a reduction's majority at least `OVERVIEW_SHARE`).
 */
export function tileImage(
  v: DrawView,
  s: RowsState,
  c: RowTiles,
  b: number,
  ly: number,
  lx: number,
  ty: number,
  tx: number,
  ov: boolean,
): number {
  const level = ly * L + lx + (ov ? OVERVIEW : 0);
  const index = (ty * Math.ceil(s.columns / (TILE << lx)) + tx) | 0;
  const m = c.levels[level];
  const held = m === undefined ? undefined : m.get(index);
  if (held !== undefined) {
    c.stamp[held] = c.clock;
    return c.img[held];
  }
  const exact = s.codes !== null && s.slot !== null && ly === 0 && lx === 0;
  const tile = exact || !v.stream ? null : v.stream.tile(ly, lx, ty, tx);
  if ((!exact && tile === null) || !v.axis) return -1;
  const theme = v.color;
  const ax = v.axis;
  const codes = s.codes as Uint8Array;
  const slots = s.slot as Int32Array;
  const C = s.columns;
  px.fill(0);
  for (let y = 0; y < TILE; y++) {
    const k = (ty * TILE + y) << ly;
    if (k >= s.lanes) break;
    const slot = exact ? slots[k] : -1;
    const row = ly === 0 ? rowAt(s, k) : -1;
    for (let x = 0; x < TILE; x++) {
      const col = (tx * TILE + x) << lx;
      if (col >= C) break;
      const cls = exact ? (slot >= 0 ? codes[slot * C + col] : 0) : (tile as { cls: Uint8Array }).cls[y * TILE + x];
      if (cls === 0) continue;
      const s0 = exact ? 255 : (tile as { share: Uint8Array }).share[y * TILE + x];
      const share = ov ? (OVERVIEW_SHARE + ((255 - OVERVIEW_SHARE) * s0 + 127) / 255) | 0 : s0;
      if ((cls === 45 || cls === 46) && (ly > 0 || lx > 0)) {
        px[y * TILE + x] = pixel(GAPS, ov ? share : share >> 1);
        continue;
      }
      const loc = ThemeLocation.at(v.loc, ax, v.layer, b, col);
      loc.row = row;
      loc.letter = s.map !== null ? s.map[cls & 127] : cls;
      const color = theme ? theme.value(loc) : GAPS;
      if (color < 0) continue;
      const alpha = theme?.alpha ? theme.alpha(loc) : OPAQUE;
      px[y * TILE + x] = pixel(color, ((share * alpha + 127) / 255) | 0);
    }
  }
  return keep(c, level, index, (v.backend.images as BackendImages).fromPixels(TILE, TILE, px));
}

/** A made image into a slot: a free one, else the least recently drawn one not drawn this frame, else a new one. */
function keep(c: RowTiles, level: number, index: number, img: number): number {
  let k = c.used;
  if (k === c.img.length) {
    k = -1;
    for (let j = 0; j < c.used; j++) if (c.stamp[j] < c.clock && (k < 0 || c.stamp[j] < c.stamp[k])) k = j;
    if (k >= 0) {
      c.levels[c.level[k]]?.delete(c.index[k]);
      c.images?.dispose(c.img[k]);
    } else {
      k = c.used;
      grow(c);
    }
  }
  if (k === c.used) c.used++;
  c.img[k] = img;
  c.level[k] = level;
  c.index[k] = index;
  c.stamp[k] = c.clock;
  let m = c.levels[level];
  if (m === undefined) c.levels[level] = m = new Map();
  m.set(index, k);
  c.builds++;
  return img;
}
/** More slots than `SLOTS`: a frame drew more tiles than the LRU holds (edge rate). */
function grow(c: RowTiles): void {
  const n = c.img.length * 2;
  const g = (a: Int32Array) => {
    const b = new Int32Array(n);
    b.set(a);
    return b;
  };
  c.img = g(c.img).fill(-1, c.used);
  c.level = g(c.level);
  c.index = g(c.index);
  c.stamp = g(c.stamp);
}

/**
 * Where `paintTiles` draws, CSS px (a double passed to a call is boxed): lane `k0`'s top (0), a lane's height (1), and
 * a single row's inset top and bottom (2; 0 for dense lanes).
 */
export const rowGeo = new Float64Array(3);
const r = new Float64Array(8);
let runs = new Int32Array(64);
const pend = new RectBatch(16);

/**
 * Lanes [k0, k0 + n) at level (ly, lx) over block b's shown columns [lo, hi) (linear), lane k0's top at `rowGeo[0]`.
 * Each tile is drawn at its whole tile's scale (its columns placed through `v.xd` from its first, its lanes from k0's
 * top), clipped to the run and the lanes, so a column's pixels never depend on where the view starts (a copied band
 * equals a redraw at any zoom). A tile not arrived is counted; with `pending` a color (the overview's) it is filled
 * with it, and with `under` a coarser level (ly << 4 | lx) that level is drawn in its place: a tile's pixels depend on
 * its own state alone, so a region redrawn on an arrival equals a full draw.
 * `open`: the lanes above k0 (1) and below k0 + n (2) are outside the caller's clip, as the columns before a `lo` above
 * 0 and after a `hi` below the block's end are: a part reaching past such an edge is not clipped there (a clip per part
 * costs the canvas a saved state). Returns the tiles not arrived.
 */
export function paintTiles(
  v: DrawView,
  s: RowsState,
  c: RowTiles,
  b: number,
  ly: number,
  lx: number,
  k0: number,
  n: number,
  lo: number,
  hi: number,
  pending: number,
  under = -1,
  open = 0,
): number {
  const d = v.dpr;
  const end = (v.axis as Axis).shown[b];
  const inset = n === 1 && rowGeo[2] > 0;
  const ya = Math.round((rowGeo[0] + (inset ? rowGeo[2] : 0)) * d) | 0;
  const yb = Math.round((rowGeo[0] + (inset ? rowGeo[1] - rowGeo[2] : n * rowGeo[1])) * d) | 0;
  let missing = 0;
  pend.reset();
  pend.begin(pending);
  let count = v.shownToRuns(b, lo, hi, runs);
  if (2 * count > runs.length) {
    runs = new Int32Array(4 * count);
    count = v.shownToRuns(b, lo, hi, runs);
  }
  for (let q = 0, ks = lo; q < count; q++) {
    const i0 = runs[2 * q];
    const len = runs[2 * q + 1] - i0;
    // A run read in descending order (a reversed piece) goes column by column, a sample over its column.
    const up = len === 1 || v.indexOf(b, ks as never) === i0;
    for (let j = 0; j < (up ? 1 : len); j++) {
      const k = ks + j;
      const m = up ? len : 1;
      const ia = up ? i0 : v.indexOf(b, k as never);
      const x0 = v.xd(b, k);
      const x1 = v.xd(b, k + m);
      const free = (k === lo && lo > 0 ? 1 : 0) | (k + m === hi && hi < end ? 2 : 0) | (open << 2);
      missing += tiles(v, s, c, b, ly, lx, k0, n, k, ia, m, up, x0, x1, ya, yb, pending, under, free);
    }
    ks += len;
  }
  if (pend.n > 0) v.backend.rects(pend);
  return missing;
}

/**
 * The tiles at (ly, lx) over lanes [k0, k0 + n) and columns [ia, ia + m) (shown k on), within device px [x0, x1) ×
 * [ya, yb): each at its whole tile's scale, clipped to its part but at the `free` edges (left, right, top, bottom).
 */
function tiles(
  v: DrawView,
  s: RowsState,
  c: RowTiles,
  b: number,
  ly: number,
  lx: number,
  k0: number,
  n: number,
  k: number,
  ia: number,
  m: number,
  up: boolean,
  x0: number,
  x1: number,
  ya: number,
  yb: number,
  pending: number,
  under: number,
  free: number,
): number {
  const d = v.dpr;
  // A single row's sample is drawn over the row (its inset), not at a tile's scale.
  const one = n === 1 && ly === 0 && rowGeo[2] > 0;
  let missing = 0;
  for (let ty = k0 >> (ly + 8); ty <= (k0 + n - 1) >> (ly + 8); ty++) {
    // Integers (`| 0`): Math.max and min may hand a double, boxed when passed to `parts`.
    const cy0 = one ? ya : Math.max(ya, Math.round((rowGeo[0] + (((ty * TILE) << ly) - k0) * rowGeo[1]) * d)) | 0;
    const cy1 = one ? yb : Math.min(yb, Math.round((rowGeo[0] + ((((ty + 1) * TILE) << ly) - k0) * rowGeo[1]) * d)) | 0;
    if (cy1 <= cy0) continue;
    for (let tx = ia >> (lx + 8); tx <= (ia + m - 1) >> (lx + 8); tx++) {
      const cx0 = up ? Math.max(x0, v.xd(b, k + (((tx * TILE) << lx) - ia))) | 0 : x0;
      const cx1 = up ? Math.min(x1, v.xd(b, k + ((((tx + 1) * TILE) << lx) - ia))) | 0 : x1;
      if (cx1 <= cx0) continue;
      const img = tileImage(v, s, c, b, ly, lx, ty, tx, pending !== NO_COLOR);
      if (img < 0) {
        missing++;
        if (pending !== NO_COLOR) pend.push(cx0, cy0, cx1 - cx0, cy1 - cy0);
        if (under >= 0)
          tiles(v, s, c, b, under >> 4, under & 15, k0, n, k, ia, m, up, cx0, cx1, cy0, cy1, NO_COLOR, -1, 0);
        continue;
      }
      if (pending === NO_COLOR) parts(v, img, b, ly, lx, k0, k, ia, up, one, ty, tx, cx0, cx1, cy0, cy1, free);
      else {
        // The overview's map is drawn whole, never copied: the run's samples of the tile in one image, as before.
        const w0 = Math.max(k0 >> ly, ty * TILE) - ty * TILE;
        const w1 = Math.min((k0 + n - 1) >> ly, ty * TILE + TILE - 1) - ty * TILE + 1;
        const u0 = Math.max(ia >> lx, tx * TILE) - tx * TILE;
        const u1 = Math.min((ia + m - 1) >> lx, tx * TILE + TILE - 1) - tx * TILE + 1;
        const xa = v.xd(b, k + (((tx * TILE + u0) << lx) - ia));
        const ya = Math.round((rowGeo[0] + (((ty * TILE + w0) << ly) - k0) * rowGeo[1]) * d) | 0;
        r[0] = u0;
        r[1] = w0;
        r[2] = u1 - u0;
        r[3] = w1 - w0;
        r[4] = xa;
        r[5] = ya;
        r[6] = v.xd(b, k + (((tx * TILE + u1) << lx) - ia)) - xa;
        r[7] = (Math.round((rowGeo[0] + (((ty * TILE + w1) << ly) - k0) * rowGeo[1]) * d) | 0) - ya;
        v.backend.image(img, r, false);
      }
    }
  }
  return missing;
}

/** About the most device px a tile's part covers (`parts`). */
const PART = 512;

/**
 * Tile (ty, tx)'s image within [cx0, cx1) × [cy0, cy1), device px, in parts of an odd count of samples, each about
 * `PART` device px at most, its edges its samples' (columns through `v.xd`, lanes from k0's top: the tile's scale
 * whatever the view's first column or row). Nearest sampling maps a pixel's center into a part: with an odd count over
 * a whole px it never lands on a sample's edge, and a small part keeps it further from one than float rounding reaches,
 * so the sample a pixel shows does not depend on where the view starts. A part over the clip's edge is clipped, but
 * at a `free` edge (one the caller's clip cuts already).
 */
function parts(
  v: DrawView,
  img: number,
  b: number,
  ly: number,
  lx: number,
  k0: number,
  k: number,
  ia: number,
  up: boolean,
  one: boolean,
  ty: number,
  tx: number,
  cx0: number,
  cx1: number,
  cy0: number,
  cy1: number,
  free: number,
): void {
  const d = v.dpr;
  // Samples a part: odd, at least one; a single row's or a reversed column's sample alone.
  const fx = Math.floor(PART / (v.cell * v.span[b] * d * (1 << lx))) | 0;
  const fy = Math.floor(PART / (rowGeo[1] * d * (1 << ly))) | 0;
  const px = !up || fx < 1 ? 1 : fx >= TILE ? TILE - 1 : fx - 1 + (fx & 1);
  const py = one || fy < 1 ? 1 : fy >= TILE ? TILE - 1 : fy - 1 + (fy & 1);
  const r0 = one ? k0 - ty * TILE : 0;
  for (let a = r0, h = 0; a < (one ? r0 + 1 : TILE); a += h) {
    h = Math.min(py, TILE - a) | 0;
    if ((h & 1) === 0) h--;
    const ya = one ? cy0 : Math.round((rowGeo[0] + (((ty * TILE + a) << ly) - k0) * rowGeo[1]) * d) | 0;
    const yb = one ? cy1 : Math.round((rowGeo[0] + (((ty * TILE + a + h) << ly) - k0) * rowGeo[1]) * d) | 0;
    if (yb <= cy0) continue;
    if (ya >= cy1) break;
    const c0 = up ? 0 : (ia >> lx) - tx * TILE;
    for (let u = c0, w = 0; u < (up ? TILE : c0 + 1); u += w) {
      w = Math.min(px, TILE - u) | 0;
      if ((w & 1) === 0) w--;
      const xa = up ? v.xd(b, k + (((tx * TILE + u) << lx) - ia)) : cx0;
      const xb = up ? v.xd(b, k + (((tx * TILE + u + w) << lx) - ia)) : cx1;
      if (xb <= cx0) continue;
      if (xa >= cx1) break;
      r[0] = u;
      r[1] = a;
      r[2] = w;
      r[3] = h;
      r[4] = xa;
      r[5] = ya;
      r[6] = xb - xa;
      r[7] = yb - ya;
      const clip =
        (xa < cx0 && !(free & 1)) ||
        (xb > cx1 && !(free & 2)) ||
        (ya < cy0 && !(free & 4)) ||
        (yb > cy1 && !(free & 8));
      if (clip) v.backend.pushClip(cx0, cy0, cx1 - cx0, cy1 - cy0);
      v.backend.image(img, r, false);
      if (clip) v.backend.popClip();
    }
  }
}

/**
 * Lanes [k0, k0 + n) of block b at the levels `lv` (`blockLevels`: fine ly, lx, coarse ly, lx): the fine level, and the
 * coarse one in place of each fine tile in view not arrived.
 */
export function paintLevels(
  v: DrawView,
  s: RowsState,
  c: RowTiles,
  b: number,
  lv: Int32Array,
  k0: number,
  n: number,
  lo: number,
  hi: number,
  open = 0,
): void {
  const fine = s.codes !== null && lv[0] === 0 && lv[1] === 0;
  const coarse = lv[2] !== lv[0] || lv[3] !== lv[1];
  paintTiles(v, s, c, b, lv[0], lv[1], k0, n, lo, hi, NO_COLOR, !fine && coarse ? (lv[2] << 4) | lv[3] : -1, open);
}
