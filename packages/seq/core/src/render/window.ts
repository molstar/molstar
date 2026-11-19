/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The window mark (_spec/README.md "Render", Window; Navigators): what a ring shows of a linear detail, a navigator's
// window and a ruler rail's, drawn in the `window` pass (the first of the emphasis, under the highlights): a ring's by
// its layout's chrome in a strip's frame, through its warp (on a ring a sector over the rings); a navigator's and a
// rail's on the page. One geometry and one hit for every window: its body in `seq-window`, a grip as a tab just outside
// each end in `seq-grip`, lit per part; a floor of 16 px; a hit of at least 12 px (24 by touch); idle (nothing can
// move) drawn quiet and not a target. A navigator's window adds what is outside it (`outside`: a click centres, a drag
// sweeps), the wash over what is not in view, and a hairline at the true extent where the floor drew it wider.

import type { RenderConfig } from './config.js';
import type { PassWriter } from './draw-lists.js';
import { NoColor } from './draw-lists.js';
import type { WindowPart } from './layout/types.js';
import type { SeqTokens } from './tokens.js';

/** A window's look and hit, CSS px (fixed, as other hit boxes are). */
export const WindowMetrics = {
  /** The least width a window is drawn at: a narrower one grows about its middle. */
  floor: 16,
  /** A grip's width (a tab just outside an end) and its share of the window's height. */
  grip: 5,
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

// ---- a navigator's window, along one axis or two (_spec/README.md "Render", Navigators) ----

type WindowConfig = RenderConfig['navigator']['window'];

/**
 * A window along one axis after the floor (`navigator.window`): its true ends `[a, b)` within the span `[lo, hi)`, drawn
 * at least `minPx` wide (at most half the span), `centred` on the truth and clamped inside the span, or placed as a
 * scrollbar's `thumb` (its position proportional to the scroll). Into `out[i]`, `out[i + 1]`.
 */
export function floorExtent(
  a: number,
  b: number,
  lo: number,
  hi: number,
  config: WindowConfig,
  out: Float64Array,
  i = 0,
): void {
  const span = hi - lo;
  const min = Math.min(config.minPx, span / 2);
  if (b - a >= min) {
    out[i] = a;
    out[i + 1] = b;
    return;
  }
  let s: number;
  if (config.floor === 'thumb') {
    const room = span - (b - a);
    s = room > 0 ? lo + ((a - lo) * (span - min)) / room : lo;
  } else s = (a + b) / 2 - min / 2;
  s = Math.max(lo, Math.min(hi - min, s));
  out[i] = s;
  out[i + 1] = s + min;
}

/**
 * The part of a navigator's window under `x` along one axis (its drawn ends `[a, b)`, its span `[lo, hi)`): an end's
 * grip (a tab `handlePx` just outside it, wider by touch, and at most a quarter of the window inside, so a press in the
 * middle pans), the body (hit over at least `hitPx`, `touchHitPx` by touch, centred and inside the span), what is
 * `outside` it within the span; undefined past the span.
 */
export function navigatorPart(
  x: number,
  a: number,
  b: number,
  lo: number,
  hi: number,
  config: WindowConfig,
  touch = false,
): WindowPart | undefined {
  const handle = touch ? Math.max(config.handlePx, config.touchHitPx / 4) : config.handlePx;
  if (x < lo - handle || x >= hi + handle) return undefined;
  const w = b - a;
  const inside = Math.min(config.handlePx, w / 4);
  if (x >= a - handle && x < a + inside) return 'start';
  if (x >= b - inside && x < b + handle) return 'end';
  const hit = Math.min(hi - lo, Math.max(w, touch ? config.touchHitPx : config.hitPx));
  const ha = Math.max(lo, Math.min(hi - hit, (a + b) / 2 - hit / 2));
  if (x >= ha && x < ha + hit) return 'body';
  return x >= lo && x < hi ? 'outside' : undefined;
}

/**
 * What `drawNavigatorWindow` reads, set by the caller (numbers in an array are never boxed: nothing allocated once
 * warm), at the indices of `NavWin`: the drawn ends along x and y, the true ends, the span the wash covers.
 */
export const NavWindow = new Float64Array(12);
export const NavWin = {
  x0: 0,
  x1: 1,
  y0: 2,
  y1: 3,
  tx0: 4,
  tx1: 5,
  ty0: 6,
  ty1: 7,
  lx: 8,
  hx: 9,
  ly: 10,
  hy: 11,
} as const;

/** How a navigator's window is drawn: the axes it has grips on (it moves along them), an export's outline, idle. */
export const NavWindowFlags = { x: 1, y: 2, outline: 4, idle: 8 } as const;

/**
 * Draws a navigator's window from `NavWindow`: the wash over the span outside it (`seq-wash`), its body tinted and
 * outlined in `seq-window`, a hairline at the true extent along an axis the floor widened (when at least a device
 * pixel), the grips of the axes it moves along (`flags`, lit by `lit`: over an end its grip, over the body both); idle
 * a quiet outline, no wash, no grip; an export's an outline with its wash, no grips.
 */
export function drawNavigatorWindow(
  out: PassWriter,
  flags: number,
  lit: WindowPart | undefined,
  dpr: number,
  config: WindowConfig,
  tokens: SeqTokens,
): void {
  const W = NavWindow;
  const g = out.args;
  const ink = tokens.colors['seq-window'];
  const idle = (flags & 8) !== 0;
  const outline = (flags & 4) !== 0;
  const x0 = W[0];
  const x1 = W[1];
  const y0 = W[2];
  const y1 = W[3];
  if (!idle) {
    const wash = tokens.colors['seq-wash'];
    washRect(out, W[8], W[10], x0 - W[8], W[11] - W[10], wash);
    washRect(out, x1, W[10], W[9] - x1, W[11] - W[10], wash);
    washRect(out, x0, W[10], x1 - x0, y0 - W[10], wash);
    washRect(out, x0, y1, x1 - x0, W[11] - y1, wash);
  }
  const fill = lit === 'body' ? WindowMetrics.litFill : WindowMetrics.fill;
  g[0] = x0;
  g[1] = y0;
  g[2] = x1 - x0;
  g[3] = y1 - y0;
  g[4] = idle || outline ? None : ((ink & 0xffffff00) | Math.round(fill * 255)) >>> 0;
  g[5] = idle ? ((ink & 0xffffff00) | Math.round(0.4 * 255)) >>> 0 : ink;
  g[6] = !idle && lit !== undefined && lit !== 'outside' ? 2 : 1;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
  if (idle) return;
  // the true extent where the floor drew it wider: a hairline, when it is at least a device pixel
  const px = 1 / dpr;
  if (W[5] - W[4] < x1 - x0 - px && W[5] - W[4] >= px) out.fill((W[4] + W[5]) / 2 - px / 2, y0, px, y1 - y0, ink);
  if (W[7] - W[6] < y1 - y0 - px && W[7] - W[6] >= px) out.fill(x0, (W[6] + W[7]) / 2 - px / 2, x1 - x0, px, ink);
  if (outline) return;
  const grip = tokens.colors['seq-grip'];
  const hw = config.handlePx;
  if (flags & 1) {
    const gh = (y1 - y0) * GripShare;
    for (let k = 0; k < 2; k++) {
      const on = lit === 'body' || (k === 0 ? lit === 'start' : lit === 'end');
      const w = on ? hw + 1 : hw;
      g[0] = k === 0 ? x0 - w : x1;
      g[1] = y0 + (y1 - y0 - gh) / 2;
      g[2] = w;
      g[3] = gh;
      g[4] = on ? ink : grip;
      g[5] = None;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
  }
  if (flags & 2) {
    const gw = (x1 - x0) * GripShare;
    for (let k = 0; k < 2; k++) {
      const on = lit === 'body' || (k === 0 ? lit === 'start' : lit === 'end');
      const h = on ? hw + 1 : hw;
      g[0] = x0 + (x1 - x0 - gw) / 2;
      g[1] = k === 0 ? y0 - h : y1;
      g[2] = gw;
      g[3] = h;
      g[4] = on ? ink : grip;
      g[5] = None;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
  }
}

function washRect(out: PassWriter, x: number, y: number, w: number, h: number, color: number): void {
  if (w > 0 && h > 0) out.fill(x, y, w, h, color);
}
