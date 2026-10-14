/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The circular geometry's map (seqspec:viewer#the-render-backend): a section drawn on rings is laid out and drawn as
// its axis unrolled along the outer ring (`SectionAxisView.unrolled`: x along the circumference, the entries' y as
// depth from the outer edge), and this maps that frame to the circle: x to an angle clockwise from `theta0` (the top),
// y to a radius. Draws go through `RingBackend` (draw/ring.ts), which maps their batches; hits come back through
// `toUnrolled`. Headless and allocation-free.

import type { AxisView, ShownIndex } from "./axis.js";
import type { SectionLayout } from "./layout.js";

const TAU = 2 * Math.PI;

/** `RingView` (draw/view.ts) of a section's rings, and the polar map both ways. Fill with `from` per frame. */
export class Ring {
  /** The center, device px; the angle of unrolled x `left` (−π/2: the top), radians. */
  cx = 0;
  cy = 0;
  theta0 = -Math.PI / 2;
  /** The view drawing an entry (its `top` and `height`): `rOuter` and `rInner` are its ring, device px. */
  drawn: { readonly top: number; readonly height: number } | null = null;
  /** Angular columns: the outer circumference in device px. */
  columns = 0;
  /** The outermost ring's inner radius, device px: text anchored outside it is a label of the ring (the ruler's). */
  edge = 0;
  /** The frame, CSS px: the outer radius, the radius per px of depth, the unrolled x and y origins, the circumference. */
  r = 0;
  scale = 1;
  left = 0;
  top = 0;
  w = 0;
  dpr = 1;
  /** The axis view the unrolled x is read from. */
  view: AxisView | null = null;

  from(l: SectionLayout, v: AxisView, dpr: number): this {
    this.dpr = dpr;
    this.cx = l.ringCx * dpr;
    this.cy = l.ringCy * dpr;
    this.r = l.ringR;
    this.scale = l.ringScale;
    this.left = l.dataLeft;
    this.top = l.scrolledTop;
    this.w = l.ringW > 0 ? l.ringW : 1;
    this.columns = Math.ceil(l.ringW * dpr);
    const k = l.pinnedCount;
    this.edge = k < l.entryCount ? (this.r - (l.entryTop[k] + l.entryHeight[k]) * this.scale) * dpr : 0;
    this.view = v;
    return this;
  }
  get rOuter(): number {
    const v = this.drawn;
    return v ? (this.r - (v.top - this.top) * this.scale) * this.dpr : 0;
  }
  get rInner(): number {
    const v = this.drawn;
    return v ? (this.r - (v.top + v.height - this.top) * this.scale) * this.dpr : 0;
  }
  /** Unrolled x (device px) → radians. */
  angleAt(xd: number): number {
    return this.theta0 + (TAU * (xd / this.dpr - this.left)) / this.w;
  }
  /** Surface y (device px) → radius (device px). */
  radiusAt(yd: number): number {
    return (this.r - (yd / this.dpr - this.top) * this.scale) * this.dpr;
  }
  angleOf(block: number, shown: ShownIndex): number {
    const v = this.view;
    return v ? this.angleAt(v.x(block, shown) * this.dpr) : this.theta0;
  }
  /** [block, shown] under an angle into `out`; false: a block gap. */
  shownAt(angle: number, out: Int32Array): boolean {
    const v = this.view;
    const ax = v?.axis;
    if (!v || !ax) return false;
    const x = this.left + (((((angle - this.theta0) % TAU) + TAU) % TAU) / TAU) * this.w;
    for (let b = 0; b < ax.count; b++) {
      const k = Math.floor((x - v.blockX[b]) / (v.span[b] * v.cell));
      if (k < 0 || k >= ax.shown[b]) continue;
      out[0] = b;
      out[1] = k;
      return true;
    }
    return false;
  }
  /**
   * A surface point (CSS px) into the unrolled frame (`out[0]` x, `out[1]` y, CSS px): the angle along the axis, the
   * depth from the outer ring down. Points off the rings land outside the entries (y < top, or past the last).
   */
  toUnrolled(x: number, y: number, out: Float64Array): void {
    const dx = x - this.cx / this.dpr;
    const dy = y - this.cy / this.dpr;
    const a = Math.atan2(dy, dx);
    out[0] = this.left + (((((a - this.theta0) % TAU) + TAU) % TAU) / TAU) * this.w;
    out[1] = this.top + (this.r - Math.sqrt(dx * dx + dy * dy)) / this.scale;
  }
  /**
   * `toUnrolled` of a record's point (`p.x`, `p.y`): read from its fields, so no double crosses a call (a pointer-rate
   * caller's; CG3).
   */
  toUnrolledAt(p: { readonly x: number; readonly y: number }, out: Float64Array): void {
    const dx = p.x - this.cx / this.dpr;
    const dy = p.y - this.cy / this.dpr;
    const a = Math.atan2(dy, dx);
    out[0] = this.left + (((((a - this.theta0) % TAU) + TAU) % TAU) / TAU) * this.w;
    out[1] = this.top + (this.r - Math.sqrt(dx * dx + dy * dy)) / this.scale;
  }
  /** A point of the unrolled frame (CSS px) onto the surface (`out[0]` x, `out[1]` y, CSS px): `toUnrolled`'s inverse. */
  toSurface(x: number, y: number, out: Float64Array): void {
    const a = this.theta0 + (TAU * (x - this.left)) / this.w;
    const r = Math.max(0, this.r - (y - this.top) * this.scale);
    out[0] = this.cx / this.dpr + r * Math.cos(a);
    out[1] = this.cy / this.dpr + r * Math.sin(a);
  }
}
