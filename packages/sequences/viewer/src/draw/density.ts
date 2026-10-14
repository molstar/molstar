/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Density composed per device column (seqspec:viewer#extension-points `density`;
// seqspec:sources#levels-and-the-zoom-gate), part of the renderer API: below a pixel `ranges`, `points`, `variants`,
// `emphasis` and pair ends draw, per device column of the data area, their items' counts per color slot: `stacked` each
// color's share of one bar, in palette order from the bottom, the bar's height by the column's total; `majority` the
// plurality color at its share of that height. Never one color's bars over another's. Items or density is one decision
// per layer per frame, for every block, with hysteresis (`densityMode`); in density a block whose count pyramid serves
// its columns adds their bins (`columnEnd`, `columnBins`), any other block (filtered, too few positions a column, a
// rows track) its items, into the same columns: one band, no item rects. Allocation-free once warm: module scratch,
// grown only.

import { type Axis, AxisGeometry } from "../axis.ts";
import type { Features, StreamFields } from "../prepare.ts";
import { binCounts, segWindow } from "../pyramids.ts";
import { densityOn, ITEM_BASE, ItemDensity, SLOTS } from "../reductions.ts";
import { type Theme, ThemeLocation } from "../themes.ts";
import { NO_COLOR, type RectBatch } from "./backend.ts";
import { columnBins, columnEnd, columnLevel, shownHull, visibleRange } from "./cells.ts";
import type { SizeView } from "./sizes.ts";
import type { DrawView } from "./view.ts";

/** The layer's `density` param (seqspec:format#layer-types): `stacked` (the default) or `majority`. */
export const DensityKind = { stacked: 0, majority: 1 } as const;
export const densityKindOf = (p: unknown): number =>
  (p as { density?: unknown } | null)?.density === "majority" ? DensityKind.majority : DensityKind.stacked;

/**
 * Feature k's color key: its color theme's at its location, its alpha folded in where it fades (bits 24 up hold 255 −
 * alpha, so an opaque color is its own key and NO_COLOR stays −1); NO_COLOR without a theme. Leaves `v.loc` at k.
 */
export function featureKey(v: DrawView, f: Features, k: number): number {
  const t = v.color;
  if (!t || !v.axis) return NO_COLOR;
  const loc = featureLoc(v, f, k);
  // A streamed item has no row in the set: a field theme reads its region's table, its values kept (ST3b).
  const s = f.fields;
  if (s != null && t.raw !== undefined) return s.theme === t ? (s.keys as Int32Array)[k] : fieldValues(f, s, t, k);
  const c = t.value(loc);
  const a = c >= 0 && t.alpha ? t.alpha(loc) : 255;
  return packAlpha(c, a);
}
/**
 * A color with its alpha folded in as a key (bits 24 up hold 255 − alpha; opaque: the color itself), every renderer's
 * one rule (LC-06). Alpha 0 is kept at 1: 255 in the top byte over a near-white color would land on a skin token's key
 * or NO_COLOR (VL-21).
 */
export function packAlpha(c: number, a: number): number {
  return a === 255 ? c : c | ((255 - (a > 0 ? a : 1)) << 24);
}
/**
 * Feature k's glyph: its glyph theme's at its location, a streamed item's from its region's table (`Theme.raw`, ST4b);
 * `none` without a glyph theme. Leaves `v.loc` at k.
 */
export function featureGlyph(v: DrawView, f: Features, k: number, none: number): number {
  const t = v.glyph;
  if (!t || !v.axis) return none;
  const loc = featureLoc(v, f, k);
  const s = f.fields;
  if (s != null && t.raw !== undefined)
    return s.glyphTheme === t ? (s.glyphs as Int32Array)[k] : glyphValues(f, s, t, k);
  return t.value(loc);
}
/** Edge rate: a region's items' glyphs by `t`; feature k's. */
function glyphValues(f: Features, s: StreamFields, t: Theme, k: number): number {
  const raw = t.raw as NonNullable<Theme["raw"]>;
  s.glyphs = Int32Array.from({ length: f.count }, (_, q) => raw(s.table, f.row[q]));
  s.glyphTheme = t;
  return s.glyphs[k];
}
/** Edge rate (a region's first draw under a theme): its items' values of `t` (`Theme.raw`); feature k's. */
function fieldValues(f: Features, s: StreamFields, t: Theme, k: number): number {
  const raw = t.raw as NonNullable<Theme["raw"]>;
  s.keys = Int32Array.from({ length: f.count }, (_, q) => raw(s.table, f.row[q]));
  s.theme = t;
  return s.keys[k];
}
/** Feature k's location in `v.loc` (its first position, annotation, row and letter), for its themes. */
export function featureLoc(v: SizeView, f: Features, k: number): ThemeLocation {
  const loc = ThemeLocation.at(v.loc, v.axis as Axis, v.layer, f.block[k], f.start[k]);
  loc.annotation = f.annotation[k];
  loc.row = f.row[k];
  loc.letter = f.code[k] > 0 ? f.code[k] : -1;
  return loc;
}
/** A key's color (keys from −64 to −1: skin tokens and NO_COLOR, as they are). */
export const keyColor = (key: number): number => (key >= -64 && key < 0 ? key : key & 0xffffff);
export const keyAlpha = (key: number): number => (key >= -64 && key < 0 ? 255 : 255 - (key >>> 24));

/** Per (prepared state): its `ItemDensity` for a theme and each slot's palette rank (the theme's legend order). */
interface Slots {
  theme: unknown;
  point: boolean;
  d: ItemDensity;
  rank: Uint8Array;
}
const slotsOf = new WeakMap<object, Slots>();
/** Edge-rate, apart from the draw: the closures here would make their contexts on every frame. */
function newSlots(v: DrawView, f: Features, owner: object, point: boolean): Slots {
  const d = new ItemDensity(f, (k) => featureKey(v, f, k), point);
  const legend: number[] = [];
  v.color?.legend?.({ item: (_l, color) => legend.push(color), ramp: () => {}, size: () => {} });
  const at = (s: number) => {
    const i = legend.indexOf(keyColor(d.color[s]));
    return i < 0 ? legend.length + s : i;
  };
  const order = Array.from({ length: d.slots }, (_, s) => s).sort((a, b) => at(a) - at(b));
  const x = { theme: v.color, point, d, rank: Uint8Array.from(order) };
  slotsOf.set(owner, x);
  return x;
}

// The columns: per device column of the data area and slot, a count; raw coverage as differences, summed at the draw.
let counts = new Float64Array(SLOTS * 1024);
let cover = new Float64Array(SLOTS * 1025);
const keys = new Int32Array(SLOTS);
/** Slots in palette order, bottom up. */
const rank = new Uint8Array(SLOTS);
const column = new Float64Array(SLOTS);
const state = new Int32Array(4);
const LEFT = 0;
const N = 1;
const USED = 2;
const COVER = 3;

/** Starts a band over the view's data area: its columns cleared, no slots (`densitySlot` assigns them as seen). */
export function densityStart(v: DrawView): void {
  const left = Math.round(v.left * v.dpr) | 0;
  const n = Math.max(1, Math.round(v.right * v.dpr) - left) | 0;
  if (counts.length < n * SLOTS) {
    counts = new Float64Array(2 * n * SLOTS);
    cover = new Float64Array((2 * n + 1) * SLOTS);
  }
  counts.fill(0, 0, n * SLOTS);
  state[LEFT] = left;
  state[N] = n;
  state[USED] = state[COVER] = 0;
}
/** The data area's device columns: `ranges`' threshold (density above one item a column). */
export const densityColumns = (v: DrawView): number =>
  Math.max(1, Math.round(v.right * v.dpr) - Math.round(v.left * v.dpr)) | 0;
/** A key's slot, assigned in order of first sight (colors past 8 share the last). */
export function densitySlot(key: number): number {
  let s = 0;
  while (s < state[USED] && keys[s] !== key) s++;
  if (s < state[USED]) return s;
  if (s === SLOTS) return SLOTS - 1;
  keys[s] = key;
  rank[s] = s;
  state[USED]++;
  return s;
}
/** `n` items (an int) of slot s in the column at device x (dropped outside the data area). */
export function densityAdd(x: number, s: number, n: number): void {
  const c = x - state[LEFT];
  if (c >= 0 && c < state[N]) counts[c * SLOTS + s] += n;
}
/** One item of slot s covering device columns [x0, x1] (coverage: each column it spans counts it once). */
export function densityCover(x0: number, x1: number, s: number): void {
  const n = state[N];
  const c0 = Math.max(0, x0 - state[LEFT]);
  const c1 = Math.min(n - 1, x1 - state[LEFT]);
  if (c1 < c0) return;
  if (state[COVER] === 0) cover.fill(0, 0, (n + 1) * SLOTS);
  state[COVER] = 1;
  cover[c0 * SLOTS + s]++;
  cover[(c1 + 1) * SLOTS + s]--;
}

let shown = new Int32Array(8);
const hull = new Int32Array(2);
const win = new Int32Array(2);
/** Block b's items over shown [k0, k1) counted raw into the columns: a point once at its center, a range per column. */
function addItems(v: DrawView, f: Features, slot: Uint8Array, point: boolean, b: number, k0: number, k1: number) {
  shownHull(v, b, k0, k1, hull);
  for (let lane = 0; lane < f.lanes; lane++) {
    const g = b * f.lanes + lane;
    segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], hull[0], hull[1], win);
    for (let k = win[0]; k < win[1]; k++) {
      const stop = Math.max(f.stop[k], f.start[k] + 1);
      if (stop <= hull[0] || (f.on[k] >= 0 && f.on[k] !== v.row)) continue;
      if (point) {
        const q = v.shownOf(b, (f.start[k] + stop - 1) >> 1);
        if (q >= k0 && q < k1) densityAdd(v.xd(b, q), slot[k], 1);
        continue;
      }
      let m = v.runToShown(b, f.start[k], stop, shown);
      if (2 * m > shown.length) {
        shown = new Int32Array(4 * m);
        m = v.runToShown(b, f.start[k], stop, shown);
      }
      for (let r = 0; r < m; r++) {
        const a = Math.max(k0, shown[2 * r]);
        const e = Math.min(k1, shown[2 * r + 1]);
        if (a < e) densityCover(v.xd(b, a), v.xd(b, e - 1), slot[k]);
      }
    }
  }
}

const range = new Int32Array(2);
const bins = new Int32Array(2);
/**
 * A feature layer's band into `rects`, between device rows y0 and y1: every block's items counted per device column
 * and color slot (`point`: once at its center, summed; else coverage, a column its bins' maximum), from the count
 * pyramid where its level serves the block's columns, else its items; drawn by `kind` (`DensityKind`). `owner` keys the
 * slots (the prepared state); a slot with no color draws `none` (emphasis: the accent; NO_COLOR: nothing).
 */
export function drawDensity(
  v: DrawView,
  f: Features,
  owner: object,
  point: boolean,
  kind: number,
  none: number,
  rects: RectBatch,
  y0: number,
  y1: number,
): void {
  const ax = v.axis;
  if (!ax) return;
  const x = slotsOf.get(owner);
  const sl =
    x?.theme === v.color && x.point === point && x.d.slot.length === f.count ? x : newSlots(v, f, owner, point);
  const d = sl.d;
  densityStart(v);
  for (let s = 0; s < d.slots; s++) {
    keys[s] = d.color[s];
    rank[s] = sl.rank[s];
  }
  state[USED] = d.slots;
  for (let b = 0; b < ax.count && (b + 1) * f.lanes < f.seg.length; b++) {
    visibleRange(v, b, range);
    if (range[0] >= range[1]) continue;
    const pyr =
      v.positions[b] === null && v.geometry === AxisGeometry.linear ? d.pyramid(b, ax.base[b], ax.length[b]) : null;
    const level = pyr ? columnLevel(v, b, pyr.counts.length, ITEM_BASE) : -1;
    if (!pyr || level < 0) {
      addItems(v, f, d.slot, point, b, range[0], range[1]);
      continue;
    }
    for (let k = range[0]; k < range[1]; ) {
      const k1 = columnEnd(v, b, k);
      if (columnBins(v, b, k, k1, level, bins)) {
        binCounts(pyr, level, bins[0], bins[1], !point, column);
        // Inline: a count is a double, and a double passed to a call V8 does not inline is boxed.
        const c = v.xd(b, k) - state[LEFT];
        if (c >= 0 && c < state[N]) for (let s = 0; s < d.slots; s++) counts[c * SLOTS + s] += column[s];
      } else addItems(v, f, d.slot, point, b, k, k1);
      k = k1;
    }
  }
  densityDraw(rects, y0, y1, kind, none);
}

/**
 * The band of the columns counted since `densityStart` into `rects`, between device rows y0 and y1: per column a bar
 * by the square root of its total over the frame's largest, `stacked` in palette order from the bottom (every color
 * present at least a device px) or `majority` (the plurality color, ties to the earlier in the palette, at its share).
 */
export function densityDraw(rects: RectBatch, y0: number, y1: number, kind: number, none: number): void {
  const n = state[N];
  const S = state[USED];
  if (state[COVER] === 1)
    for (let s = 0; s < S; s++)
      for (let c = 0, run = 0; c < n; c++) {
        run += cover[c * SLOTS + s];
        counts[c * SLOTS + s] += run;
      }
  let max = 0;
  for (let c = 0; c < n; c++) {
    let t = 0;
    for (let s = 0; s < S; s++) t += counts[c * SLOTS + s];
    if (t > max) max = t;
  }
  const h = y1 - y0;
  for (let c = 0; max > 0 && h > 0 && c < n; c++) {
    let t = 0;
    let nz = 0;
    let best = -1;
    for (let r = 0; r < S; r++) {
      const s = rank[r];
      const m = counts[c * SLOTS + s];
      if (m === 0) continue;
      t += m;
      nz++;
      if (best < 0 || m > counts[c * SLOTS + best]) best = s;
    }
    if (t === 0) continue;
    const bar = Math.min(h, Math.max(2, nz, Math.round(h * Math.sqrt(t / max))));
    const x = state[LEFT] + c;
    if (kind === DensityKind.majority) {
      const top = Math.max(1, Math.round((bar * counts[c * SLOTS + best]) / t));
      slotColor(rects, best, none);
      rects.push(x, y1 - top, 1, top);
      continue;
    }
    let below = 0;
    let sum = 0;
    for (let r = 0, j = 0; r < S; r++) {
      const s = rank[r];
      if (counts[c * SLOTS + s] === 0) continue;
      sum += counts[c * SLOTS + s];
      // At most the bar less a row for each color still to come (so at most the bar).
      const top = Math.min(bar - (nz - ++j), Math.max(below + 1, Math.round((bar * sum) / t)));
      // More colors than the band has rows: those left no row are not drawn, never out of the band (VL-13).
      if (top <= below) continue;
      slotColor(rects, s, none);
      rects.push(x, y1 - top, 1, top - below);
      below = top;
    }
  }
}
function slotColor(rects: RectBatch, s: number, none: number) {
  const key = keys[s];
  rects.begin(key === NO_COLOR ? none : keyColor(key), keyAlpha(key));
}

/**
 * The items of a feature layer in view (every block drawn), counted to `limit`: what decides items or density (ranges
 * over the data area's device columns, points and variants over 64).
 */
export function itemsInView(v: DrawView, f: Features, limit: number): number {
  const ax = v.axis;
  let n = 0;
  // The whole view's, in a pan's strip too (`viewLo`): one decision for the frame, whatever part of it is drawn.
  const lo = v.visLo;
  const hi = v.visHi;
  v.visLo = v.viewLo;
  v.visHi = v.viewHi;
  for (let b = 0; ax && b < ax.count && (b + 1) * f.lanes < f.seg.length && n < limit; b++) {
    visibleRange(v, b, range);
    if (range[0] >= range[1]) continue;
    shownHull(v, b, range[0], range[1], hull);
    for (let lane = 0; lane < f.lanes && n < limit; lane++) {
      const g = b * f.lanes + lane;
      segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], hull[0], hull[1], win);
      if (!f.hasRows) n += win[1] - win[0];
      else for (let k = win[0]; k < win[1] && n < limit; k++) if (f.on[k] < 0 || f.on[k] === v.row) n++;
    }
  }
  v.visLo = lo;
  v.visHi = hi;
  return n < limit ? n : limit;
}

/**
 * Items or density for the layer drawn, once per frame for every block: density past `threshold` items in view, items
 * again only below 0.8 of it, kept per layer and view so a resize near it flips once.
 */
export function densityMode(v: DrawView, count: number, threshold: number): boolean {
  const m = v.densityModes;
  const key = v.layer * 65536 + v.onlyBlock + 1;
  const on = densityOn(count, threshold, m.get(key) === 1);
  m.set(key, on ? 1 : 0);
  return on;
}
/** Whether the layer drawn was last drawn as density (`densityMode`'s decision, per view; an overlay reads its main pass's). */
export const densityShown = (v: DrawView): boolean => v.densityModes.get(v.layer * 65536 + v.onlyBlock + 1) === 1;
