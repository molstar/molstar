/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The window of the navigator, the vertical overview and the scrollbar (seqspec:viewer#behaviors,
// seqspec:interaction#selection): where a point falls on it. One classification for the gesture a press starts there
// and the cursor shown before the press,
// so what the cursor says is what a drag does. Headless and allocation-free: a frame is computed from the section's
// axis view and layout, a part from a frame. The hit region's kind, the lit part and the rail's grab zone read plain
// records the caller passes (`ui.window.current`, `ui.drag.value`, the pointer; IN2).

import type { AxisView } from './axis.js';
import { HitRegion } from './hit.js';

/** `pan` drags the window; the ends resize it (zoom); `page` jumps (the scrollbar pages). */
export const WindowPart = { none: 0, pan: 1, left: 2, right: 3, top: 4, bottom: 5, page: 6 } as const;
/** `windowFrame`'s kinds. */
export const WindowKind = { navigator: 0, overview: 1, scrollbar: 2, ring: 3 } as const;
/** Px of an end's handle: a tab just outside the window's end, and at most this (and a quarter of it) inside. */
export const WINDOW_EDGE = 5;
/**
 * A window or thumb is drawn at least this long (at most half its span, so it still travels): its true share, the
 * zoom limit's share (which the view cannot pass) or this floor, whichever is longest. A floored window is a classic
 * scrollbar thumb: its position maps the scroll over the span less the thumb, so a drag covers the whole range.
 */
export const WINDOW_FLOOR = 16;
/** A window is hit over at least this many px, centered on it (inside its span); its ends' handles just outside. */
export const WINDOW_HIT = 12;
/** The numbers in a window's frame (`windowFrame`). */
export const WINDOW_FRAME = 18;

/** What `windowFrame` reads of the section's layout (LY1's `SectionLayout` has these). */
export interface WindowLayout {
  readonly width: number;
  readonly height: number;
  readonly scrolledTop: number;
  readonly vnavWidth: number;
}

/**
 * One axis of a window into `g` at `o` from its span at `s` and fractions at `f`: its drawn ends in px (at least
 * `WINDOW_FLOOR` long, its position proportional over the travel), the fraction one px of travel moves (at `per`), and
 * the ends it is hit between (at `hit`: at least `WINDOW_HIT` px, centered on it, inside the span). Indices only: a
 * double passed to a call V8 does not inline is boxed.
 */
function span(g: Float64Array, o: number, per: number, hit: number, s: number, f: number) {
  const lo = g[s];
  const len = g[s + 1];
  const f0 = g[f];
  const frac = Math.min(1, Math.max(0, g[f + 1] - f0));
  const w = Math.min(len, Math.max(Math.min(WINDOW_FLOOR, len / 2), frac * len));
  const travel = len - w;
  const p0 = lo + (frac < 1 ? (Math.max(0, f0) / (1 - frac)) * travel : 0);
  g[o] = Math.min(p0, lo + travel);
  g[o + 1] = g[o] + w;
  g[per] = travel > 0 ? (1 - frac) / travel : 0;
  const hw = Math.min(len, Math.max(WINDOW_HIT, w));
  g[hit] = Math.max(lo, Math.min(lo + len - hw, g[o] + (w - hw) / 2));
  g[hit + 1] = g[hit] + hw;
}

/**
 * The window's frame for a section into `g` (`WINDOW_FRAME` numbers): span x0, width, y0, height; window fractions
 * x0, x1, y0, y1; the window as drawn, px (x0, x1, y0, y1: at least `WINDOW_FLOOR` long); fractions per px of travel along x and y
 * (a drag of the body moves the view by them); the window as hit, px (x0, x1, y0, y1: at least `WINDOW_HIT` long).
 * `kind` (`WindowKind`) 0: the navigator (x over its map of the axis), 1: the vertical overview (x and y over its
 * column), 2: the scrollbar (y).
 */
export function windowFrame(v: AxisView, l: WindowLayout, kind: number, g: Float64Array): void {
  const dw = v.navWidth > 0 ? v.navWidth : v.right - v.left;
  const rh = l.height - l.scrolledTop;
  g[0] = kind === WindowKind.overview ? l.width - l.vnavWidth : v.left;
  g[1] = kind === WindowKind.overview ? l.vnavWidth : dw;
  g[2] = l.scrolledTop;
  g[3] = rh;
  // No data area (the header and data columns fill the width) shows no columns: a share of none, floored at the start
  //, never 0 / 0's NaN, which drew the overview's window as a line and washed the rows in view.
  g[4] = dw > 0 ? (v.navX0 - v.left) / dw : 0;
  g[5] = dw > 0 ? (v.navX1 - v.left) / dw : 0;
  g[6] = rh > 0 ? (v.vnavY0 - l.scrolledTop) / rh : 0;
  g[7] = rh > 0 ? (v.vnavY1 - l.scrolledTop) / rh : 1;
  windowRect(g);
}
/** From a frame's span and fractions (0..7), the window as drawn (8..11), its travel (12, 13) and as hit (14..17). */
export function windowRect(g: Float64Array) {
  span(g, 8, 12, 14, 0, 4);
  span(g, 10, 13, 16, 2, 6);
}

/**
 * An end of the window dragged `total` px from the press, over the press's frame `g` (along y with `vertical`; the
 * start's end with `start`): only that end moves, by the pointer's travel (no snap to the pointer, so nothing jumps),
 * the other held; it stops at the span's end and where the drawn window would be shorter than `least` of the span
 * (the zoom limit) or than the floor (`WINDOW_FLOOR`), never flipping or pushing the other end. The window's new
 * fractions into `io` (0, 1), from its drawn ends (above the floor, drawn and real ends agree; a floored thumb's
 * end dragged out zooms out from there); false: the ends as at the press (a floored thumb dragged in keeps the view).
 * `total` and `least` are read from `io[0]` and `io[1]`: a window drag calls this per move, and a double passed to a
 * call V8 does not inline is boxed.
 */
export function windowEnds(g: Float64Array, vertical: boolean, start: boolean, io: Float64Array): boolean {
  const total = io[0];
  const least = io[1];
  const lo = g[vertical ? 2 : 0];
  const len = g[vertical ? 3 : 1];
  const a0 = g[vertical ? 10 : 8];
  const b0 = g[vertical ? 11 : 9];
  const min = Math.min(b0 - a0, Math.max(Math.min(WINDOW_FLOOR, len / 2), least * len));
  const a = start ? Math.max(lo, Math.min(b0 - min, a0 + total)) : a0;
  const b = start ? b0 : Math.min(lo + len, Math.max(a0 + min, b0 + total));
  if (len <= 0 || (Math.abs(a - a0) < 1e-6 && Math.abs(b - b0) < 1e-6)) return false;
  io[0] = (a - lo) / len;
  io[1] = (b - lo) / len;
  return true;
}

/** The end zone of the point along x or y (frame `g`, its hit ends): -1 the start's, 1 the end's, 0 neither. */
function end(g: Float64Array, p: { readonly x: number; readonly y: number }, vertical: boolean): number {
  const at = vertical ? p.y : p.x;
  const a = g[vertical ? 16 : 14];
  const b = g[vertical ? 17 : 15];
  const lo = g[vertical ? 2 : 0];
  const hi = lo + g[vertical ? 3 : 1];
  const inner = Math.min(WINDOW_EDGE, (b - a) / 4);
  if (at >= Math.max(lo, a - WINDOW_EDGE) && at < a + inner) return -1;
  if (at > b - inner && at <= Math.min(hi, b + WINDOW_EDGE)) return 1;
  return 0;
}
/**
 * The `WindowPart` under a point (a hit or an input record: a double passed to a call is boxed, a record's field
 * is not) on a window of `kind`, from its frame `g`. The handles are tabs just outside the window's ends and at most a
 * quarter of it inside, so a press in its middle always moves it. Row zoom is a control only with `rows` (the
 * section's `verticalZoom`): else the window only moves.
 */
export function windowPart(
  kind: number,
  at: { readonly x: number; readonly y: number },
  g: Float64Array,
  rows: boolean,
): number {
  if (kind < 0) return WindowPart.none;
  const inX = at.x >= g[14] && at.x <= g[15];
  const inY = at.y >= g[16] && at.y <= g[17];
  // The navigator has no vertical ends (its frame's y is the rows' span, not its row).
  const ey = rows && (kind === WindowKind.scrollbar || (kind === WindowKind.overview && inX)) ? end(g, at, true) : 0;
  if (kind === WindowKind.scrollbar) {
    return ey < 0 ? WindowPart.top : ey > 0 ? WindowPart.bottom : inY ? WindowPart.pan : WindowPart.page;
  }
  const ex = end(g, at, false);
  if (ex !== 0 && (kind === WindowKind.navigator || inY || ey === 0))
    return ex < 0 ? WindowPart.left : WindowPart.right;
  if (ey !== 0) return ey < 0 ? WindowPart.top : WindowPart.bottom;
  return inX && (kind === WindowKind.navigator || inY) ? WindowPart.pan : WindowPart.page;
}

/**
 * The grips a lit part lights, as bits: 1 the start's, 2 the end's, of the ends along y with `vertical`, else along x.
 * The body (`pan`) lights both; outside the window (`page`) and the other axis's ends, none.
 */
export function gripsLit(part: number, vertical: boolean): number {
  if (part === WindowPart.pan) return 3;
  if (part === (vertical ? WindowPart.top : WindowPart.left)) return 1;
  if (part === (vertical ? WindowPart.bottom : WindowPart.right)) return 2;
  return 0;
}

/**
 * Whether a horizontal window has nowhere to go (the navigator's, with the rail's): the view shows the whole axis
 * (the content fits the data area) and zoom-in is at its limit, so it can neither zoom in or out nor pan. A fixed cell
 * never zooms or scrolls. Plain records in (the viewport's state and its axis), no allocation; it is false again as soon
 * as a narrower view or more columns make zoom possible.
 */
export function windowIdle(
  s: { readonly cell: number; readonly maxCell: number; readonly contentWidth: number; readonly width: number },
  ax: { readonly insetX: number; readonly fixedCell: number } | null,
): boolean {
  if (!ax || s.width <= 0) return false;
  return s.cell >= s.maxCell * (1 - 1e-9) && (ax.fixedCell > 0 || s.contentWidth <= s.width - ax.insetX + 1e-6);
}

/** The window kind (`WindowKind`) of a hit region; -1: none. */
export const windowKind = (region: number) =>
  region === HitRegion.navigator
    ? WindowKind.navigator
    : region === HitRegion.vnav
      ? WindowKind.overview
      : region === HitRegion.scrollbar
        ? WindowKind.scrollbar
        : -1;

/** What `railGrown` reads of the layout (LY1's `SectionLayout` has these). */
export interface RailLayout {
  readonly railEntry: number;
  readonly railTop: number;
  readonly dataLeft: number;
  readonly dataRight: number;
  readonly entryTop: ArrayLike<number>;
  readonly entryHeight: ArrayLike<number>;
}
/**
 * The rail navigator drawn grown (thicker): while its drag runs, wherever the pointer goes, and else while the pointer
 * is over its grab zone (the `rulerRail`'s, layout.ts `railTop`). The same rule as the cursor staying `grabbing`.
 * Plain records in, so the surface reads it per frame: `ui.drag.value`, the section's layout, `ctx.pointer.current`.
 */
export function railGrown(
  drag: { readonly section: number; readonly owner: string } | null,
  l: RailLayout | null,
  p: { readonly section: number; readonly x: number; readonly y: number; readonly inside: boolean },
  section: number,
): boolean {
  if (drag && drag.section === section) return drag.owner === 'navigator';
  if (!l || l.railEntry < 0 || !p.inside || p.section !== section) return false;
  const bottom = l.entryTop[l.railEntry] + l.entryHeight[l.railEntry];
  return p.y >= l.railTop && p.y < bottom && p.x >= l.dataLeft && p.x < l.dataRight;
}

/** The part lit on a section's window of `kind` (`ui.window.current`, set with the cursor); `none` where another is lit. */
export function litPart(
  w: { readonly section: number; readonly kind: number; readonly part: number },
  section: number,
  kind: number,
): number {
  return w.section === section && w.kind === kind ? w.part : WindowPart.none;
}
