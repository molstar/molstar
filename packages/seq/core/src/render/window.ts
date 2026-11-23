/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The window mark (_spec/README.md "Render", Window): what a ring shows of a linear detail, a navigator's window and a
// ruler rail's, drawn in the `window` pass (the first of the emphasis, under the highlights): a ring's by its layout's
// chrome in a strip's frame, through its warp; a navigator's and a rail's on the page. One geometry and one hit for
// every window, read from `navigator.window`: its floor (`floorExtent`), its parts (`navigatorPart`: a grip as a tab
// just outside each end, the body, what is `outside` it), its drawing (`drawNavigatorWindow`: the body in
// `seq-window`, the grips in `seq-grip`, lit per part, the wash over what is not in view, a hairline at the true extent
// where the floor drew it wider; idle a quiet outline, not a target).

import { DefaultRenderConfig, type RenderConfig } from './config.js';
import type { PassWriter } from './draw-lists.js';
import { NoColor } from './draw-lists.js';
import type { WindowPart } from './layout/types.js';
import type { SeqTokens } from './tokens.js';

type WindowConfig = RenderConfig['navigator']['window'];

/**
 * The zoom strip (_spec/README.md "Render", Navigators, The zoom strip): a reachable span narrower than `minPx` on the
 * navigator gets a zoom strip under it, `heightPx` more (its connectors' gap, `ZoomGapPx`, within it). The viewer's
 * `navigator.zoomStrip` sets them; these are the defaults (an export's).
 */
export const DefaultZoomStrip: {
  readonly minPx: number;
  readonly heightPx: number;
  readonly overBlocks: ZoomOverBlocks;
} = { minPx: 48, heightPx: 16, overBlocks: 'always' };

/**
 * When a navigator over blocks (a genome) has the zoom strip (_spec/README.md "Render", Navigators): `always`, the
 * chromosome the control row under the genome whatever the window's width (the default: a rule a reader can see);
 * `narrow`, only while the window is narrower than `minPx` (with hysteresis: the earlier look, for a host short of
 * height); `never`, one row with a floored window.
 */
export type ZoomOverBlocks = 'always' | 'narrow' | 'never';

/**
 * The window's look, fixed as other marks' are (a look, not a reader's choice): its body's fill opacity, and lit
 * (hovered or dragged); an idle outline's opacity; a grip's share of the window's thickness.
 */
export const WindowFill = 0.12;
const LitFill = 0.2;
const GripShare = 0.5;
// an imported binding read in a draw is a property read under a module runner (it boxed one number a frame here)
const None: number = NoColor;

/**
 * How a navigator's window is drawn: the axes it has grips on (it moves along them), an export's outline, idle; and a
 * grip left out (`noStart`, `noEnd`: a piece of a window over a cyclic origin has the grip of its true end only).
 */
export const NavWindowFlags = { x: 1, y: 2, outline: 4, idle: 8, noStart: 16, noEnd: 32, quiet: 64 } as const;
// the flags as numbers read once: a draw reads no property of an object
const FlagX: number = NavWindowFlags.x;
const FlagY: number = NavWindowFlags.y;
const FlagOutline: number = NavWindowFlags.outline;
const FlagIdle: number = NavWindowFlags.idle;
const FlagNoStart: number = NavWindowFlags.noStart;
const FlagNoEnd: number = NavWindowFlags.noEnd;
const FlagQuiet: number = NavWindowFlags.quiet;

/**
 * A window along one axis after the floor (`navigator.window`): its true ends `[a, b)` within the span `[lo, hi)` read
 * from `e[i]` to `e[i + 3]` (numbers in an array are never boxed: nothing allocated once warm), drawn at least `minPx`
 * wide (at most half the span), `centred` on the truth and clamped inside the span, or placed as a scrollbar's
 * `thumb` (its position proportional to the scroll). The drawn ends into `e[i]`, `e[i + 1]`. An unbounded span (a
 * ring's frame) only grows the window about its middle.
 */
export function floorExtent(e: Float64Array, i: number, config: WindowConfig): void {
  const a = e[i];
  const b = e[i + 1];
  const lo = e[i + 2];
  const hi = e[i + 3];
  const span = hi - lo;
  const min = Math.min(config.minPx, span / 2);
  if (b - a >= min) return;
  let s: number;
  if (config.floor === 'thumb' && span < Number.POSITIVE_INFINITY) {
    const room = span - (b - a);
    s = room > 0 ? lo + ((a - lo) * (span - min)) / room : lo;
  } else s = (a + b) / 2 - min / 2;
  s = Math.max(lo, Math.min(hi - min, s));
  e[i] = s;
  e[i + 1] = s + min;
}

/**
 * The part of a window under `x` along one axis (its drawn ends `[a, b)`, its span `[lo, hi)`): an end's grip (a tab
 * just outside it, hit over half of `hitPx` there, of `touchHitPx` by touch, at least `handlePx`; and at most a quarter
 * of the window inside it, so a press in the middle pans), the body (hit over at least `hitPx`, centred and inside the
 * span), what is `outside` it within the span; undefined past the span.
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
  const hit = touch ? config.touchHitPx : config.hitPx;
  const out = Math.max(config.handlePx, hit / 2);
  if (x < lo - out || x > hi + out) return undefined;
  const w = b - a;
  const inside = Math.min(out, w / 4);
  if (x >= a - out && x < a + inside) return 'start';
  if (x > b - inside && x <= b + out) return 'end';
  const h = Math.min(hi - lo, Math.max(w, hit));
  const ha = Math.max(lo, Math.min(hi - h, (a + b) / 2 - h / 2));
  if (x >= ha && x < ha + h) return 'body';
  return x >= lo && x < hi ? 'outside' : undefined;
}

const Extent = new Float64Array(4);

/** Floors `Extent[0]`, `Extent[1]` with no span (a ring's frame has no ends). */
function floorFree(config: WindowConfig): void {
  Extent[2] = Number.NEGATIVE_INFINITY;
  Extent[3] = Number.POSITIVE_INFINITY;
  floorExtent(Extent, 0, config);
}

/**
 * A ring's window's part under a point of its strip's frame (`[x0, x1)` across `top`, `h`), or undefined off it: the
 * shared parts over the floored window, with no span.
 */
export function windowPart(
  x: number,
  y: number,
  x0: number,
  x1: number,
  top: number,
  h: number,
  touch = false,
  config: WindowConfig = DefaultRenderConfig.navigator.window,
): WindowPart | undefined {
  if (y < top || y >= top + h) return undefined;
  Extent[0] = x0;
  Extent[1] = x1;
  floorFree(config);
  const part = navigatorPart(x, Extent[0], Extent[1], Extent[2], Extent[3], config, touch);
  return part === 'outside' ? undefined : part;
}

/**
 * A ring's window with its ends and band in `out.args` (x0, x1, y, h, in a strip's frame): the shared drawing over
 * the floored window, its grips along x, no wash (a ring washes its sector itself); idle a quiet outline.
 */
export function drawWindowArgs(
  out: PassWriter,
  active: boolean,
  lit: WindowPart | undefined,
  tokens: SeqTokens,
  config: WindowConfig = DefaultRenderConfig.navigator.window,
): void {
  const g = out.args;
  const W = NavWindow;
  W[4] = g[0];
  W[5] = g[1];
  W[6] = g[2];
  W[7] = g[2] + g[3];
  // the ends through an array: a number passed to a call is boxed (nothing allocated once warm)
  Extent[0] = g[0];
  Extent[1] = g[1];
  floorFree(config);
  W[0] = Extent[0];
  W[1] = Extent[1];
  W[2] = W[6];
  W[3] = W[7];
  // the span is the window itself: no wash
  W[8] = W[0];
  W[9] = W[1];
  W[10] = W[2];
  W[11] = W[3];
  drawNavigatorWindow(out, active ? FlagX : FlagIdle, active ? lit : undefined, 1, config, tokens);
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

/**
 * The rect in `out.args` (x, y, width, height, fill, stroke, stroke width), written unless empty: the values reach the
 * pass through the array (a number passed to a call that is not inlined is boxed).
 */
function emit(out: PassWriter): void {
  const g = out.args;
  if (!(g[2] > 0) || !(g[3] > 0)) return;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * Draws a window from `NavWindow`: the wash over the span outside it (`seq-wash`), its body tinted and outlined in
 * `seq-window`, a hairline at the true extent along an axis the floor widened (when at least a device pixel), the
 * grips of the axes it moves along (`flags`, lit by `lit`: over an end its grip, over the body both); idle a quiet
 * outline, no wash, no grip; an export's an outline with its wash, no grips.
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
  const idle = (flags & FlagIdle) !== 0;
  const outline = (flags & FlagOutline) !== 0;
  const x0 = W[0];
  const x1 = W[1];
  const y0 = W[2];
  const y1 = W[3];
  if (!idle) {
    const wash = tokens.colors['seq-wash'];
    // the four sides; an empty one is skipped
    for (let k = 0; k < 4; k++) {
      g[0] = k === 0 ? W[8] : k === 1 ? x1 : x0;
      g[1] = k < 3 ? W[10] : y1;
      g[2] = k === 0 ? x0 - W[8] : k === 1 ? W[9] - x1 : x1 - x0;
      g[3] = k < 2 ? W[11] - W[10] : k === 2 ? y0 - W[10] : W[11] - y1;
      g[4] = wash;
      g[5] = 0;
      g[6] = 0;
      emit(out);
    }
  }
  const fill = lit === 'body' ? LitFill : WindowFill;
  g[0] = x0;
  g[1] = y0;
  g[2] = x1 - x0;
  g[3] = y1 - y0;
  const quiet = idle || (flags & FlagQuiet) !== 0;
  g[4] = idle || outline ? None : ((ink & 0xffffff00) | Math.round(fill * 255)) >>> 0;
  // two stores, not a conditional between two numbers (one is boxed a frame)
  g[5] = ink;
  if (quiet) g[5] = ((ink & 0xffffff00) | Math.round(config.idleAlpha * 255)) >>> 0;
  g[6] = config.strokePx;
  emit(out);
  if (idle) return;
  // the true extent where the floor drew it wider: a hairline, when it is at least a device pixel
  const px = 1 / dpr;
  if (W[5] - W[4] < x1 - x0 - px && W[5] - W[4] >= px) {
    g[0] = (W[4] + W[5]) / 2 - px / 2;
    g[1] = y0;
    g[2] = px;
    g[3] = y1 - y0;
    g[4] = ink;
    g[5] = 0;
    g[6] = 0;
    emit(out);
  }
  if (W[7] - W[6] < y1 - y0 - px && W[7] - W[6] >= px) {
    g[0] = x0;
    g[1] = (W[6] + W[7]) / 2 - px / 2;
    g[2] = x1 - x0;
    g[3] = px;
    g[4] = ink;
    g[5] = 0;
    g[6] = 0;
    emit(out);
  }
  if (outline) return;
  const grip = tokens.colors['seq-grip'];
  const hw = config.handlePx;
  if (flags & FlagX) {
    const gh = (y1 - y0) * GripShare;
    for (let k = 0; k < 2; k++) {
      if (flags & (k === 0 ? FlagNoStart : FlagNoEnd)) continue;
      const on = lit === 'body' || (k === 0 ? lit === 'start' : lit === 'end');
      const w = on ? hw + 1 : hw;
      g[0] = k === 0 ? x0 - w : x1;
      g[1] = y0 + (y1 - y0 - gh) / 2;
      g[2] = w;
      g[3] = gh;
      g[4] = on ? ink : grip;
      g[5] = None;
      g[6] = 0;
      emit(out);
    }
  }
  if (flags & FlagY) {
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
      emit(out);
    }
  }
}
