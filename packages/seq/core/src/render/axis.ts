/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The axis transform (_spec/README.md "Render", Tracks layout): a section's positions to x. The positions shown are
// pieces in display order (`[start, end, …]`, an axis's shape: two when a cyclic range wraps), laid end to end as
// display units; the viewport is a window of display units drawn across the data box. Cell edges snap to the device
// grid.

/** A window of display units: from `start` to `end` (fractional), along a section's pieces end to end. */
export interface Viewport {
  readonly start: number;
  readonly end: number;
}

/** Which ends of a span the view or a piece's edge cut (they are drawn square). */
export const Cut = { left: 1, right: 2 } as const;

export class AxisTransform {
  /** Each piece's first display unit. */
  readonly offsets: number[] = [];
  readonly total: number;
  readonly u0: number;
  readonly u1: number;
  /** CSS px per position. */
  readonly scale: number;

  constructor(
    readonly system: string,
    /** `[start, end, …]` in display order. */
    readonly pieces: readonly number[],
    readonly x0: number,
    readonly width: number,
    readonly dpr: number,
    viewport?: Viewport,
  ) {
    let total = 0;
    for (let i = 0; i + 1 < pieces.length; i += 2) {
      this.offsets.push(total);
      total += pieces[i + 1] - pieces[i];
    }
    this.total = total;
    const u0 = viewport ? Math.max(0, Math.min(viewport.start, total)) : 0;
    const u1 = viewport ? Math.max(u0, Math.min(viewport.end, total)) : total;
    this.u0 = u0;
    this.u1 = u1 > u0 ? u1 : u0 + 1;
    this.scale = width / (this.u1 - this.u0);
  }

  /** Device pixels per position: what the density steps read. */
  get devicePx(): number {
    return this.scale * this.dpr;
  }

  /** The data box's device columns. */
  get columns(): number {
    return Math.max(1, Math.ceil(this.width * this.dpr));
  }

  /** x of a display unit. */
  xOfUnit(u: number): number {
    return this.x0 + (u - this.u0) * this.scale;
  }

  /** x snapped to the device grid. */
  snap(x: number): number {
    return Math.round(x * this.dpr) / this.dpr;
  }

  /** The device column of x (clamped to the box). */
  column(x: number): number {
    const c = Math.floor((x - this.x0) * this.dpr);
    return c < 0 ? 0 : c >= this.columns ? this.columns - 1 : c;
  }

  /** The piece holding position `p` (`[s, e)`); with `edge`, a position at a piece's end counts as in it. */
  pieceOf(p: number, edge = false): number {
    const n = this.pieces.length;
    for (let k = 0; k + 1 < n; k += 2) {
      const s = this.pieces[k];
      const e = this.pieces[k + 1];
      if (p >= s && (p < e || (edge && p === e))) return k >> 1;
    }
    return -1;
  }

  /** x of the left edge of position `p`, or NaN where it is not shown. */
  x(p: number): number {
    const k = this.pieceOf(p, true);
    if (k < 0) return Number.NaN;
    return this.xOfUnit(this.offsets[k] + p - this.pieces[2 * k]);
  }

  /**
   * The cell of position `p`: whether both its edges are shown, its left and right x (`x(p)`, `x(p + 1)`, unsnapped) in
   * `sx0` and `sx1`. No number is returned: the hover's draw calls it.
   */
  cell(p: number): boolean {
    const k0 = this.pieceOf(p, true);
    if (k0 < 0) return false;
    const k1 = this.pieceOf(p + 1, true);
    if (k1 < 0) return false;
    this.sx0 = this.xOfUnit(this.offsets[k0] + p - this.pieces[2 * k0]);
    this.sx1 = this.xOfUnit(this.offsets[k1] + (p + 1) - this.pieces[2 * k1]);
    return Number.isFinite(this.sx0) && Number.isFinite(this.sx1);
  }

  /** The position under x (an integer), or −1 where there is none. */
  positionAt(x: number): number {
    if (x < this.x0 || x >= this.x0 + this.width) return -1;
    const u = this.u0 + (x - this.x0) / this.scale;
    for (let k = 0; k < this.offsets.length; k++) {
      const s = this.pieces[2 * k];
      const len = this.pieces[2 * k + 1] - s;
      if (u >= this.offsets[k] && u < this.offsets[k] + len) return s + Math.floor(u - this.offsets[k]);
    }
    return -1;
  }

  /** The last `spanIn`'s x extents (unsnapped, clipped to the data box) and which ends were cut: read right after it. */
  sx0 = 0.5;
  sx1 = 0.5;
  scut = 0;

  /**
   * Positions `[s, e)` in piece `k`: whether any of it is in view there, its x extents and cut in `sx0`, `sx1`, `scut`.
   * What `spans` makes a triple of, with no array: the layout's per-position loops call it directly.
   */
  spanIn(k: number, s: number, e: number): boolean {
    const ps = this.pieces[2 * k];
    const pe = this.pieces[2 * k + 1];
    const a = Math.max(s, ps);
    const b = Math.min(e, pe);
    // a boundary (zero length) inside a piece is a hairline at its place
    if (b < a || (b === a && (a === pe || s !== e))) return false;
    const left = this.x0;
    const right = this.x0 + this.width;
    let xa = this.xOfUnit(this.offsets[k] + a - ps);
    let xb = this.xOfUnit(this.offsets[k] + b - ps);
    // a span that only touches the view's edge is not in view; a boundary on the edge is
    if (s === e ? xb < left || xa > right : xb <= left || xa >= right) return false;
    let cut = 0;
    if (a > s || xa < left) cut |= Cut.left;
    if (b < e || xb > right) cut |= Cut.right;
    if (xa < left) xa = left;
    if (xb > right) xb = right;
    this.sx0 = xa;
    this.sx1 = xb;
    this.scut = cut;
    return true;
  }

  /**
   * The x extents of positions `[s, e)` as triples `x0, x1, cut` (unsnapped, clipped to the data box), a triple per
   * piece it touches in view; `cut` says which ends the view or a piece's edge cut. Appended to `out`.
   */
  spans(s: number, e: number, out: number[]): number[] {
    for (let k = 0; k < this.offsets.length; k++) if (this.spanIn(k, s, e)) out.push(this.sx0, this.sx1, this.scut);
    return out;
  }

  /** The positions in view, `[start, end, …]` per piece (whole cells: a partly shown one counts), appended to `out`. */
  visible(out: number[]): number[] {
    for (let k = 0; k < this.offsets.length; k++) {
      const ps = this.pieces[2 * k];
      const off = this.offsets[k];
      const len = this.pieces[2 * k + 1] - ps;
      const a = Math.max(this.u0, off);
      const b = Math.min(this.u1, off + len);
      if (b <= a) continue;
      out.push(ps + Math.floor(a - off), ps + Math.ceil(b - off));
    }
    return out;
  }
}
