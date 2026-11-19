/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The window mark (_spec/README.md "Render", Window): what a ring shows of a linear detail, and a navigator's window,
// drawn by a layout's chrome in the `window` pass (the first of the emphasis, under the highlights) in a strip's frame,
// through its warp (on a ring a sector over the rings). One geometry and one hit for every window: its body in
// `seq-window`, a grip as a tab just outside each end in `seq-grip`, lit per part; a floor of 16 px; a hit of at least
// 12 px (24 by touch); idle (nothing follows it) drawn quiet and not a target.

import type { PassWriter } from './draw-lists.js';
import { NoColor } from './draw-lists.js';
import type { WindowPart } from './layout/types.js';
import type { SeqTokens } from './tokens.js';

/** A window's look and hit, CSS px (fixed, as other hit boxes are). */
export const WindowMetrics = {
  /** The least width a window is drawn at: a narrower one grows about its middle. */
  floor: 16,
  /** A grip's width (a tab just outside an end) and its share of the window's height. */
  grip: 4,
  gripShare: 0.5,
  /** The least hit width of a part, pointer and touch. */
  hit: 12,
  touchHit: 24,
  /** The body's fill opacity, and lit (hovered or dragged). */
  fill: 0.12,
  litFill: 0.2,
} as const;

/** The window's ends after the floor: into `out[0]`, `out[1]`. */
export function windowExtent(x0: number, x1: number, out: Float64Array): void {
  const w = x1 - x0;
  if (w >= WindowMetrics.floor) {
    out[0] = x0;
    out[1] = x1;
    return;
  }
  const c = (x0 + x1) / 2;
  out[0] = c - WindowMetrics.floor / 2;
  out[1] = c + WindowMetrics.floor / 2;
}

const Extent = new Float64Array(2);

// the metrics a draw reads, as numbers read once: a draw reads no property of an imported object
const Grip: number = WindowMetrics.grip;
const GripShare: number = WindowMetrics.gripShare;
// an imported binding read in a draw is a property read under a module runner (it boxed one number a frame here)
const None: number = NoColor;

/**
 * Draws a window over `[x0, x1)` across `y`, `h` (a strip's frame): idle a quiet outline; active its body tinted and
 * outlined in `seq-window` and its grips in `seq-grip`; the `lit` part drawn grown (its body brighter, its grip wider).
 */
export function drawWindow(
  out: PassWriter,
  x0: number,
  x1: number,
  y: number,
  h: number,
  active: boolean,
  lit: WindowPart | undefined,
  tokens: SeqTokens,
): void {
  const g = out.args;
  g[0] = x0;
  g[1] = x1;
  g[2] = y;
  g[3] = h;
  drawWindowArgs(out, active, lit, tokens);
}

/**
 * `drawWindow` with its ends and band in `out.args` (x0, x1, y, h): what an emphasis pass calls, since a number passed
 * to a call that is not inlined is boxed (nothing allocated once warm).
 */
export function drawWindowArgs(out: PassWriter, active: boolean, lit: WindowPart | undefined, tokens: SeqTokens): void {
  const g = out.args;
  const y = g[2];
  const h = g[3];
  // `windowExtent` inline: its numbers stay unboxed
  let a = g[0];
  let b = g[1];
  if (b - a < WindowMetrics.floor) {
    const c = (a + b) / 2;
    a = c - WindowMetrics.floor / 2;
    b = c + WindowMetrics.floor / 2;
  }
  const ink = tokens.colors['seq-window'];
  g[0] = a;
  g[1] = y;
  g[2] = b - a;
  g[3] = h;
  // the inks at their opacities inline: a color passed to a call is boxed (nothing allocated once warm)
  const fill = lit === 'body' ? WindowMetrics.litFill : WindowMetrics.fill;
  g[4] = active ? ((ink & 0xffffff00) | Math.round(fill * 255)) >>> 0 : None;
  g[5] = active ? ink : ((ink & 0xffffff00) | Math.round(0.4 * 255)) >>> 0;
  g[6] = active && lit !== undefined ? 2 : 1;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
  if (!active) return;
  const grip = tokens.colors['seq-grip'];
  const gh = h * GripShare;
  for (let k = 0; k < 2; k++) {
    const lit1 = (k === 0 && lit === 'start') || (k === 1 && lit === 'end');
    const gw = lit1 ? Grip + 2 : Grip;
    g[0] = k === 0 ? a - gw : b;
    g[1] = y + (h - gh) / 2;
    g[2] = gw;
    g[3] = gh;
    g[4] = grip;
    g[5] = None;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
}

/**
 * The part of an active window under a point (a strip's frame), or undefined: an end within half the hit width of
 * its edge (its grip outside it included), else the body. An end wins over the body near it, so a narrow window keeps
 * both grips.
 */
export function windowPart(
  x: number,
  y: number,
  x0: number,
  x1: number,
  top: number,
  h: number,
  touch = false,
): WindowPart | undefined {
  if (y < top || y >= top + h) return undefined;
  windowExtent(x0, x1, Extent);
  const a = Extent[0];
  const b = Extent[1];
  const half = (touch ? WindowMetrics.touchHit : WindowMetrics.hit) / 2;
  const da = Math.abs(x - a);
  const db = Math.abs(x - b);
  if (da <= half || db <= half) return da <= db ? 'start' : 'end';
  return x > a && x < b ? 'body' : undefined;
}
