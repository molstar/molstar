/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The size channels in the draws (seqspec:format#layer-types, seqspec:format#themes): a layer's size themes are the
// last of its themes (the type's size channels, in the order the kernel table lists them), each value an integer in
// 1/`SIZE_ONE` of the channel's unit (px, or a share of the lane, the row or the cell), `NO_SIZE` the viewer's own
// size. Part of the renderer API.

import type { Axis, AxisView } from "../axis.ts";
import type { Features } from "../prepare.ts";
import { NO_SIZE, SIZE_ONE, type Theme, ThemeLocation } from "../themes.ts";
import { shownHull } from "./cells.ts";
import { featureLoc } from "./density.ts";

/**
 * The viewer's own dot, CSS px: the radius of a glyph no size theme sizes (a pair's end, a point at a narrow cell; 9 px
 * across, as Nightingale draws them). Explicit sizes of a document's size theme are never scaled.
 */
export const DOT_R = 4.5;
/** The largest radius a point's glyph grows to with its cell (a row's height may hold it less). */
export const DOT_R_MAX = 6;
/** A variant's mark in a stack of them (7 px across; the stack's rows keep fitting its band) and the most its cell grows it to. */
export const DOT_STACK_R = 3.5;
export const DOT_STACK_R_MAX = 8;

/** What sizes are read through: a `DrawView` in the draw, `hitLayer` in a hit. */
export interface SizeView {
  axis: Axis | null;
  layer: number;
  themes: readonly (Theme | null)[];
  readonly loc: ThemeLocation;
}
/**
 * The layer a renderer's `hit` runs for: hit.ts sets its axis, index and themes before each layer's `hit`, so a hit
 * sizes a glyph or a band as the draw does. Not re-entrant.
 */
export const hitLayer: SizeView = { axis: null, layer: -1, themes: [], loc: ThemeLocation.create() };

/**
 * Into `out` the block-index hull [lo, hi) of block b's shown positions within `reach` px of x (a hit's point): the
 * items whose glyphs may hold it. False: none.
 */
export function hitWindow(v: AxisView, b: number, x: number, reach: number, out: Int32Array): boolean {
  hitReach[0] = x;
  hitReach[1] = reach;
  return hitWindowAt(v, b, out);
}
/** `hitWindow`'s point and reach (CSS px), for a caller whose doubles must not cross a call boxed (a per-move hit). */
export const hitReach = new Float64Array(2);
/** `hitWindow` at `hitReach`'s point and reach. */
export function hitWindowAt(v: AxisView, b: number, out: Int32Array): boolean {
  const ax = v.axis;
  if (!ax) return false;
  const w = ax.span[b] * v.cell;
  const k0 = Math.max(0, Math.floor((hitReach[0] - hitReach[1] - v.blockX[b]) / w));
  const k1 = Math.min(ax.shown[b], Math.floor((hitReach[0] + hitReach[1] - v.blockX[b]) / w) + 1);
  shownHull(v, b, k0, k1, out);
  return out[1] > out[0];
}

/** The layer's size theme `k` of its type's `n` size channels (the last in its channel order); null: the viewer's size. */
export function sizeTheme(v: SizeView, k: number, n: number): Theme | null {
  const t = v.themes[v.themes.length - n + k];
  return t?.kind === "size" ? t : null;
}
/** Feature k's size by `t` (an int, see above); `NO_SIZE` without a theme. Leaves `v.loc` at k. */
export function featureSize(v: SizeView, f: Features, k: number, t: Theme | null): number {
  return t && v.axis ? t.value(featureLoc(v, f, k)) : NO_SIZE;
}
/** A size in CSS px, clamped to [least, most]; `none` where the theme gives none. */
export function sizePx(size: number, none: number, least: number, most: number): number {
  return size < 0 ? none : Math.min(most, Math.max(least, size / SIZE_ONE));
}
/** The bar `[y0, y1)` (device px) shrunk to a band share of its lane and centered, into `bar[0]` and `bar[1]`. */
export const bar = new Int32Array(2);
export function barBand(size: number, y0: number, y1: number): void {
  bar[0] = y0;
  bar[1] = y1;
  if (size < 0 || size >= SIZE_ONE) return;
  const h = Math.max(1, ((y1 - y0) * size + (SIZE_ONE >> 1)) >> 8);
  bar[0] = (((y0 + y1) >> 1) - (h >> 1)) | 0;
  bar[1] = bar[0] + h;
}
