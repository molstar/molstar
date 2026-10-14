/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ring's window (seqspec:interaction#linked-sections, CG3): a section on rings linked by `follow` to a linear one
// shows that section's view as a sector of its rings, the linear navigator's window drawn round the circle. The
// window is the other section's viewport read through the link (its center element's image and the px a position
// there, `followImage`), so it follows that section whoever moves it; dragging it pans that section, a grip zooms it
// with the other end held, both through the same link back (`centerOn`). Across a whole cyclic unit's origin the
// sector wraps. Headless and allocation-free: frames go in a reused Float64Array, points in plain records.

import type { ResolvedView } from "@molstar/seq-core/format";
import { type Axis, AxisGeometry, type AxisView, type ShownIndex } from "./axis.ts";
import type { SeqViewerContext } from "./context.ts";
import { defineBehavior } from "./default-behaviors.ts";
import { followImage } from "./follow.ts";
import { DragPhase, PointerType } from "./input.ts";
import type { SectionLayout } from "./layout.ts";
import { Ring } from "./ring.ts";
import type { Origin } from "./state.ts";
import { Viewport } from "./viewport.ts";
import { WINDOW_EDGE, WindowPart } from "./windows.ts";

/**
 * A ring window's frame (`ringWindowFrame`), unrolled CSS px along the outer ring: the window's ends x0, x1 (x0 in the
 * block, x1 past its end across a cyclic origin), the ring block's extent lo, hi, cyclic (1: a whole cyclic unit, the
 * window wraps), the ring block, px a shown position there, the px a ring position takes in the linear section's view,
 * that view's data width, its section, the link from the ring to it, the least and most window (its zoom limits), and
 * the image of its content's ends on the ring (cl, ch, cl < ch, unwrapped past the block's end; ±Infinity: the content
 * covers a whole cyclic block). Then the set's input (x0, length), the pointer's unrolled x (`ringWindowAt`) and a
 * drag's state (the pointer's last unrolled x, its travel since the press).
 */
export const RING_WINDOW = 21;
/** A frame's numbers the drawing reads: the window's ends, the block's extent, the block cyclic (1). */
export const RW_X0 = 0;
export const RW_X1 = 1;
export const RW_LO = 2;
export const RW_HI = 3;
export const RW_CYCLIC = 4;
const BLOCK = 5;
const PX = 6;
const PER = 7;
const AREA = 8;
const TARGET = 9;
const BACK = 10;
const LEAST = 11;
const MOST = 12;
const CL = 13;
const CH = 14;
const SET_X = 15;
const SET_LEN = 16;
const AT_X = 17;
const LAST = 18;
const TRAVEL = 19;
/** Scratch for the content's ends (`contentImage`). */
const CX = 20;

// Enum-like members read once (a read per event allocated under the test runner's module wrapper).
const NONE = WindowPart.none;
const PAN = WindowPart.pan;
const LEFT = WindowPart.left;
const RIGHT = WindowPart.right;
const DOWN = PointerType.down;
const END = DragPhase.end;
const LINEAR = AxisGeometry.linear;
const CIRCULAR = AxisGeometry.circular;

const io = new Float64Array(2);
const link = new Int32Array(2);
const polar = new Float64Array(2);
/** The rim's and the last ring's unrolled y (`ringBounds`). */
const bounds = new Float64Array(2);
const ring = new Ring();

/**
 * The linear section a ring's window shows: the first linked to `section` by `follow`, either way, into `link[0]` (its
 * section) and `link[1]` (the link from it to the ring, `2 * pair + side`); false: none.
 */
function linked(ctx: SeqViewerContext, section: number): boolean {
  const view = ctx.graph.assembled.value?.view;
  if (!view) return false;
  const v = view.resolved;
  const f = v.follows;
  const me = ctx.layout.sectionRef(section);
  for (let p = 0; me !== null && p < f.count; p++) {
    const a = view.set.strings.get(v.sections.ref[f.a[p]]);
    const b = view.set.strings.get(v.sections.ref[f.b[p]]);
    const other = a === me ? b : b === me ? a : null;
    if (other === null || other === me || ctx.layout.get(other)?.current.geometry !== LINEAR) continue;
    link[0] = ctx.layout.sectionIndex(other);
    link[1] = 2 * p + (other === b ? 1 : 0);
    return link[0] >= 0;
  }
  return false;
}

/** The viewport a ring section's window shows (edge rate: a surface follows it); null: none. */
export function ringWindowTarget(ctx: SeqViewerContext, section: number): Viewport | null {
  const ref = ctx.layout.sectionRef(section);
  if (ref === null || ctx.layout.get(ref)?.current.geometry !== CIRCULAR || !linked(ctx, section)) return null;
  const vp = ctx.viewports.get(ctx.layout.sectionRef(link[0]) ?? "");
  return vp instanceof Viewport ? vp : null;
}

/** Whether block b of the ring's axis is a whole cyclic unit (its window wraps across the origin). */
function wholeCyclic(ax: NonNullable<Viewport["axis"]>, b: number): boolean {
  const u = ax.unit[b];
  const m = u >= 0 ? ax.set.molecules[ax.set.units.molecule[u]] : null;
  return !!m?.cyclic && ax.base[b] === 0 && ax.length[b] === m.length;
}

/** The window of ring section `section` into `g` (`RING_WINDOW`); false: no ring, no linked linear section or view. */
export function ringWindowFrame(ctx: SeqViewerContext, section: number, g: Float64Array): boolean {
  const vp = ringWindowTarget(ctx, section);
  const ref = ctx.layout.sectionRef(section);
  const av = ref === null ? null : ctx.layout.axis(ref);
  const ax = av?.axis;
  const model = ctx.graph.assembled.value?.view;
  if (!vp?.axis || !av || !ax || !model) return false;
  const s = vp.live.current;
  const area = s.width - vp.axis.insetX;
  const b = followImage(model.resolved, link[1], vp.centerElement(io), io, ax);
  if (b < 0 || b >= ax.count || area <= 0 || !(io[1] > 0) || ax.shown[b] === 0) return false;
  const i = Math.floor(io[0]);
  let k: number = ax.shownOf(b, i);
  let fr = io[0] - i;
  if (k < 0) {
    k = ax.nearestShown(b, i);
    fr = 0.5;
  }
  if (k < 0) return false;
  const px = av.span[b] * av.cell;
  const lo = av.blockX[b];
  const hi = lo + ax.shown[b] * px;
  const len = Math.min(hi - lo, (area / io[1]) * px);
  const cyclic = wholeCyclic(ax, b);
  let x0 = lo + (ax.dir[b] > 0 ? k + fr : k + 1 - fr) * px - len / 2;
  if (!cyclic) x0 = Math.max(lo, Math.min(hi - len, x0));
  else if (x0 < lo) x0 += hi - lo;
  else if (x0 >= hi) x0 -= hi - lo;
  g[RW_X0] = x0;
  g[RW_X1] = x0 + len;
  g[RW_LO] = lo;
  g[RW_HI] = hi;
  g[RW_CYCLIC] = cyclic ? 1 : 0;
  g[BLOCK] = b;
  g[PX] = px;
  g[PER] = io[1];
  g[AREA] = area;
  g[TARGET] = link[0];
  g[BACK] = link[1] ^ 1;
  // At the zoom limits a ring position takes `maxCell / cell` (`minCell / cell`) times the px it takes now.
  g[LEAST] = Math.min(hi - lo, (area / (io[1] * (s.maxCell / s.cell))) * px);
  g[MOST] = Math.min(hi - lo, (area / (io[1] * (s.minCell / s.cell))) * px);
  contentImage(model.resolved, vp.axis, ax, av, g);
  return true;
}

/**
 * Ring block g[BLOCK]'s index io[0] (with its fraction) as unrolled x into g[CX]; false: another block or not shown.
 */
function xOnRing(ax: Axis, av: AxisView, b: number, g: Float64Array): boolean {
  const i = Math.floor(io[0]);
  const k: number = ax.shownOf(b, i);
  if (k < 0) return false;
  const fr = io[0] - i;
  g[CX] = av.blockX[b] + (ax.dir[b] > 0 ? k + fr : k + 1 - fr) * g[PX];
  return true;
}
/**
 * The image on the ring of the linear section's content, from its first shown position's start to its last's end,
 * into g[CL], g[CH] (unwrapped across a cyclic origin); ±Infinity: the whole of a cyclic block, or no image.
 */
function contentImage(v: ResolvedView, tx: Axis, ax: Axis, av: AxisView, g: Float64Array): void {
  g[CL] = Number.NEGATIVE_INFINITY;
  g[CH] = Number.POSITIVE_INFINITY;
  const b = g[BLOCK];
  const last = tx.count - 1;
  if (tx.count === 0 || tx.shown[0] === 0 || tx.shown[last] === 0) return;
  // Just inside each end (an image past the content is the nearest element's middle).
  let i = tx.indexOf(0, 0 as ShownIndex);
  io[0] = tx.dir[0] > 0 ? i + 1e-6 : i + 1 - 1e-6;
  io[1] = 1;
  if (followImage(v, g[BACK] ^ 1, 0, io, ax) !== b || !xOnRing(ax, av, b, g)) return;
  const c0 = g[CX];
  i = tx.indexOf(last, (tx.shown[last] - 1) as ShownIndex);
  io[0] = tx.dir[last] > 0 ? i + 1 - 1e-6 : i + 1e-6;
  io[1] = 1;
  if (followImage(v, g[BACK] ^ 1, last, io, ax) !== b || !xOnRing(ax, av, b, g)) return;
  let lo = Math.min(c0, g[CX]);
  let hi = Math.max(c0, g[CX]);
  const span = g[RW_HI] - g[RW_LO];
  if (g[RW_CYCLIC] === 1) {
    // In shown order from the first's image; the whole block (less a position) is no bound.
    lo = c0;
    hi = g[CX] > c0 ? g[CX] : g[CX] + span;
    if (hi - lo >= span - 1.5 * g[PX]) return;
  }
  g[CL] = lo;
  g[CH] = hi;
}

/**
 * The linked section moved so the window spans `g[SET_X]` for `g[SET_LEN]` px (unrolled): its center element through
 * the link, at the px a position that shows that span across its data width (clamped there to its zoom limits), with
 * the origin of what moved it (a drag's `pointer`, the keys' `keyboard`: seqspec:viewer#state).
 */
export function ringWindowSet(
  ctx: SeqViewerContext,
  section: number,
  g: Float64Array,
  origin: Origin = "pointer",
): void {
  const ref = ctx.layout.sectionRef(section);
  const ax = ref === null ? null : ctx.layout.axis(ref)?.axis;
  const vp = ctx.viewports.get(ctx.layout.sectionRef(g[TARGET]) ?? "");
  const model = ctx.graph.assembled.value?.view;
  if (!ax || !(vp instanceof Viewport) || !model) return;
  const b = g[BLOCK];
  const n = ax.shown[b];
  let pos = (g[SET_X] + g[SET_LEN] / 2 - g[RW_LO]) / g[PX];
  pos = g[RW_CYCLIC] === 1 ? ((pos % n) + n) % n : Math.max(0, Math.min(n - 1e-6, pos));
  const k = Math.floor(pos);
  const i = ax.indexOf(b, k as ShownIndex);
  if (i < 0) return;
  io[0] = ax.dir[b] > 0 ? i + pos - k : i + 1 - (pos - k);
  io[1] = (g[AREA] * g[PX]) / g[SET_LEN];
  const t = followImage(model.resolved, g[BACK], b, io, vp.axis);
  if (t >= 0) vp.centerOn(t, io, origin);
}

/**
 * The window's part (`WindowPart`: pan, left, right or none) under surface point `at` of ring section `section`, its
 * frame into `g` and the point's unrolled x (unwrapped past the window's start) at `g[AT_X]`: over the rings or the
 * rim just outside them, within the window or on an end's handle (`WINDOW_EDGE` along the rim, at most a quarter of the
 * window inside it), as the navigator's (windows.ts).
 */
export function ringWindowAt(
  ctx: SeqViewerContext,
  section: number,
  at: { readonly x: number; readonly y: number },
  g: Float64Array,
): number {
  const ref = ctx.layout.sectionRef(section);
  const l = ref === null ? null : ctx.layout.get(ref)?.current;
  const av = ref === null ? null : ctx.layout.axis(ref);
  if (!l || !av || l.geometry !== CIRCULAR || !ringWindowFrame(ctx, section, g)) return NONE;
  ring.from(l, av, 1).toUnrolledAt(at, polar);
  // Radially: from the last ring's inner edge to a grip's half past the rim.
  ringBounds(l, bounds);
  const rim = ring.r - (bounds[0] - ring.top) * ring.scale;
  const depth = ring.r - (polar[1] - ring.top) * ring.scale;
  if (depth > rim + 7 || depth < ring.r - (bounds[1] - ring.top) * ring.scale) return NONE;
  g[AT_X] = polar[0];
  // The handles' px along the rim, as unrolled px (the outer ring's).
  const e = (WINDOW_EDGE * ring.r) / Math.max(1, rim);
  // Off the window's block (another block of the ring), nothing; on a cyclic one, the copy of x nearest the window.
  let x = polar[0];
  if (x < g[RW_LO] - e || x > g[RW_HI] + e) return NONE;
  // Both copies computed on every call: an arithmetic branch first taken after warm-up deoptimized this (measured).
  const span = g[RW_HI] - g[RW_LO];
  const up = x + span;
  const down = x - span;
  if (g[RW_CYCLIC] === 1 && x < g[RW_X0] - e) x = up;
  else if (g[RW_CYCLIC] === 1 && x > g[RW_X1] + e && down >= g[RW_X0] - e) x = down;
  const inner = Math.min(e, (g[RW_X1] - g[RW_X0]) / 4);
  if (x >= g[RW_X0] - e && x < g[RW_X0] + inner) return LEFT;
  if (x > g[RW_X1] - inner && x <= g[RW_X1] + e) return RIGHT;
  return x >= g[RW_X0] && x <= g[RW_X1] ? PAN : NONE;
}

/**
 * The window's section zoomed by `by[0]` (above 1: in) about the window's center, within its zoom limits and its
 * content (a ring's `+` and `-`: seqspec:interaction#the-keyboard), with origin `keyboard` (GG-04). False: no window.
 */
export function ringWindowZoom(
  ctx: SeqViewerContext,
  section: number,
  by: Float64Array,
  origin: Origin = "keyboard",
): boolean {
  const g = zg;
  if (!ringWindowFrame(ctx, section, g)) return false;
  const len = g[RW_X1] - g[RW_X0];
  const n = Math.max(g[LEAST], Math.min(g[MOST], len / by[0]));
  let x0 = g[RW_X0] + (len - n) / 2;
  if (g[RW_CYCLIC] === 1 && g[CL] > Number.NEGATIVE_INFINITY && x0 < g[CL]) x0 += g[RW_HI] - g[RW_LO];
  g[SET_X] = Math.max(g[CL], Math.min(g[CH] - n, x0));
  g[SET_LEN] = n;
  ringWindowSet(ctx, section, g, origin);
  return true;
}
const zg = new Float64Array(RING_WINDOW);

/**
 * The rim the window's grips sit on (the outermost ruler's inner edge, its baseline, else the rings' top) and the last
 * ring's inner edge, unrolled y, into out[0], out[1].
 */
export function ringBounds(l: SectionLayout, out: Float64Array): void {
  const k = l.pinnedCount;
  out[0] =
    k < l.entryCount && l.ringTop > l.scrolledTop + l.entryTop[k] + 0.5
      ? l.scrolledTop + l.entryTop[k] + l.entryHeight[k]
      : l.ringTop;
  const e = l.entryCount - 1;
  out[1] = e >= 0 ? l.scrolledTop + l.entryTop[e] + l.entryHeight[e] : l.scrolledTop;
}

export const RingWindowBehavior = defineBehavior(
  "ringWindow",
  "interaction",
  100,
  {},
  "On rings, the window of the linear section linked by `follow`: dragging it pans that section, dragging a grip " +
    "zooms it with the other end held (a plain press of the primary button; a click still selects).",
  ({ ctx, h, on }) => {
    let mode: number = NONE;
    /**
     * The window's frame at the press, and the drag's state in its last slots (doubles in a Float64Array: a captured
     * `let` holding one boxes at every store); the pointer's surface point in a drag.
     */
    const g0 = new Float64Array(RING_WINDOW);
    const pt = { x: 0.5, y: 0.5 };
    on(ctx.interaction.pointer, (e) => {
      const p = e.input;
      if (p.type !== DOWN || p.button !== 0 || p.modifiers !== 0 || e.owner !== "") return;
      if (!h.forSection(e.hit.section)) return;
      const part = ringWindowAt(ctx, e.hit.section, p, g0);
      if (part === NONE) return;
      e.claim("ringWindow");
      mode = part;
      g0[LAST] = g0[AT_X];
      g0[TRAVEL] = 0;
    });
    on(ctx.interaction.drag, (e) => {
      if (e.owner !== "ringWindow") return;
      const sec = e.start.section;
      const vp = ctx.viewports.get(ctx.layout.sectionRef(g0[TARGET]) ?? "");
      if (e.phase >= END) return vp?.release();
      const ref = ctx.layout.sectionRef(sec);
      const l = ref === null ? null : ctx.layout.get(ref)?.current;
      const av = ref === null ? null : ctx.layout.axis(ref);
      if (!l || !av || !vp) return;
      pt.x = e.start.x + e.totalX;
      pt.y = e.start.y + e.totalY;
      ring.from(l, av, 1).toUnrolledAt(pt, polar);
      dragTo(g0, polar, ring.w, mode);
      ringWindowSet(ctx, sec, g0);
    });
  },
);

/**
 * A window drag's step: the pointer's unrolled x (`at[0]`) taken the shorter way round from the last (`w` the
 * circumference) into the travel, and the window it asks for into `g[SET_X]`, `g[SET_LEN]`: moved (`PAN`) or one end
 * moved with the other held, within its zoom limits, the block and the linear section's content (`cl`, `ch`; a cyclic
 * block's window first taken to the content's copy).
 */
function dragTo(g: Float64Array, at: Float64Array, w: number, mode: number): void {
  let d = at[0] - g[LAST];
  if (d > w / 2) d -= w;
  else if (d < -w / 2) d += w;
  g[LAST] = at[0];
  g[TRAVEL] += d;
  const t = g[TRAVEL];
  const span = g[RW_HI] - g[RW_LO];
  const cyclic = g[RW_CYCLIC] === 1;
  let x0 = g[RW_X0];
  const len = g[RW_X1] - x0;
  if (cyclic && g[CL] > Number.NEGATIVE_INFINITY && x0 < g[CL]) x0 += span;
  const x1 = x0 + len;
  // Where the window may go: the content's image (on a whole cyclic block, anywhere), within a linear block.
  const lo = cyclic ? g[CL] : Math.max(g[RW_LO], g[CL]);
  const hi = cyclic ? g[CH] : Math.min(g[RW_HI], g[CH]);
  if (mode === PAN) {
    g[SET_X] = Math.max(lo, Math.min(hi - len, x0 + t));
    g[SET_LEN] = len;
  } else if (mode === LEFT) {
    const n = Math.max(g[LEAST], Math.min(g[MOST], x1 - lo, len - t));
    g[SET_X] = x1 - n;
    g[SET_LEN] = n;
  } else {
    const n = Math.max(g[LEAST], Math.min(g[MOST], hi - x0, len + t));
    g[SET_X] = x0;
    g[SET_LEN] = n;
  }
}
