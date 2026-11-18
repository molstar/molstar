/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The axis transform (_spec/README.md "Render", Tracks layout): a section's positions to x and back. The positions
// shown are pieces in display order (`[start, end, …]`, start below end, an axis's shape), laid end to end as display
// units in the axis's direction: on a reversed axis a piece `[s, e)` shows `e − 1 − k` at its offset `k` (EC-32). The
// viewport is a window of display units drawn across the data box. A break between two pieces is zero units wide. Cell
// edges snap to the device grid. Lookups are binary searches over the pieces (a reversed axis of 150,000 pieces, EC-33,
// pans with no scan of them).

import { axisBreaks, monotonic, pieceAtOffset, pieceOffsets, pieceOfPosition } from '@molstar/seqx-core/view';

/** A window of display units: from `start` to `end` (fractional), along a section's pieces end to end. */
export interface Viewport {
  readonly start: number;
  readonly end: number;
}

/** Which ends of a span the view or a piece's edge cut, on screen (they are drawn square). */
export const Cut = { left: 1, right: 2 } as const;

/** What a transform reads of its axis besides its pieces: its direction, and for its breaks, its system's. */
export interface AxisShapeOptions {
  readonly dir?: 1 | -1;
  readonly length?: number;
  readonly cyclic?: boolean;
}

export class AxisTransform {
  /** Each piece's first display unit. */
  readonly offsets: Float64Array;
  readonly total: number;
  readonly u0: number;
  readonly u1: number;
  /** CSS px per position. */
  readonly scale: number;
  /** 1: positions count up left to right; −1: down. */
  readonly dir: 1 | -1;
  readonly length?: number;
  readonly cyclic: boolean;
  /** The first and last piece in view (display order). */
  readonly k0: number;
  readonly k1: number;
  /** Whether the pieces run one way in display order: a position's piece is a binary search. */
  private readonly sorted: boolean;

  constructor(
    readonly system: string,
    /** `[start, end, …]` in display order. */
    readonly pieces: readonly number[],
    readonly x0: number,
    readonly width: number,
    readonly dpr: number,
    viewport?: Viewport,
    shape: AxisShapeOptions = {},
  ) {
    const n = pieces.length >> 1;
    const all = pieceOffsets(pieces);
    this.offsets = all.subarray(0, n);
    this.dir = shape.dir ?? 1;
    if (shape.length !== undefined) this.length = shape.length;
    this.cyclic = shape.cyclic ?? false;
    this.sorted = monotonic(pieces, this.dir);
    const total = all[n];
    this.total = total;
    const u0 = viewport ? Math.max(0, Math.min(viewport.start, total)) : 0;
    const u1 = viewport ? Math.max(u0, Math.min(viewport.end, total)) : total;
    this.u0 = u0;
    this.u1 = u1 > u0 ? u1 : u0 + 1;
    this.scale = width / (this.u1 - this.u0);
    this.k0 = n ? pieceAtOffset(all, u0) : 0;
    let k1 = n ? pieceAtOffset(all, this.u1) : -1;
    // a piece starting at the window's end is not in it
    while (k1 > this.k0 && all[k1] >= this.u1) k1--;
    this.k1 = k1;
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

  /** The display unit of the boundary before position `b` in piece `k` (its left edge forward, its right reversed). */
  unitOf(k: number, b: number): number {
    return this.dir === -1 ? this.offsets[k] + this.pieces[2 * k + 1] - b : this.offsets[k] + b - this.pieces[2 * k];
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
    return pieceOfPosition(this.pieces, this.dir, p, edge);
  }

  /** x of the boundary before position `p` (the left edge of its cell forward, the right reversed), or NaN. */
  x(p: number): number {
    const k = this.pieceOf(p, true);
    if (k < 0) return Number.NaN;
    return this.xOfUnit(this.unitOf(k, p));
  }

  /**
   * The cell of position `p`: whether it is shown, its left and right x on screen (unsnapped) in `sx0` and `sx1`. No
   * number is returned: the hover's draw calls it.
   */
  cell(p: number): boolean {
    const k = this.pieceOf(p, false);
    if (k < 0) return false;
    const a = this.xOfUnit(this.unitOf(k, p));
    const b = this.xOfUnit(this.unitOf(k, p + 1));
    this.sx0 = a < b ? a : b;
    this.sx1 = a < b ? b : a;
    return Number.isFinite(this.sx0) && Number.isFinite(this.sx1);
  }

  /** The display unit under x (fractional), or −1 outside the data box. */
  unitAt(x: number): number {
    if (x < this.x0 || x >= this.x0 + this.width) return -1;
    return this.u0 + (x - this.x0) / this.scale;
  }

  /** The position under x (an integer), or −1 where there is none. */
  positionAt(x: number): number {
    const u = this.unitAt(x);
    if (u < 0) return -1;
    return this.elementAtUnit(u);
  }

  /** The position whose cell holds display unit `u`, or −1 past the pieces. */
  elementAtUnit(u: number): number {
    const n = this.offsets.length;
    if (!n || u >= this.total) return -1;
    let k = pieceAtOffset(pieceOffsets(this.pieces), u);
    // past pieces of no length at the same offset
    while (k + 1 < n && u >= this.offsets[k + 1]) k++;
    const d = Math.floor(u - this.offsets[k]);
    return this.dir === -1 ? this.pieces[2 * k + 1] - 1 - d : this.pieces[2 * k] + d;
  }

  /** The last `spanIn`'s x extents (unsnapped, clipped to the data box) and which ends were cut: read right after it. */
  sx0 = 0.5;
  sx1 = 0.5;
  scut = 0;
  /** The last `piecesIn`'s first and last piece (display order): read right after it. */
  ka = 0;
  kb = -1;

  /**
   * The pieces in view that positions `[s, e)` may touch, `ka` to `kb` (display order; none: `kb < ka`): a binary
   * search where the pieces run one way, else every piece in view. What a per-item loop over `spanIn` iterates.
   */
  piecesIn(s: number, e: number): void {
    this.ka = this.k0;
    this.kb = this.k1;
    if (!this.sorted || this.k1 < this.k0) return;
    const r = this.pieces;
    const n = r.length >> 1;
    // genomic order: the first piece ending at or after `s`, the last starting at or before `e`
    let lo = 0;
    let hi = n;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      const k = this.dir === -1 ? n - 1 - mid : mid;
      if (r[2 * k + 1] < s) lo = mid + 1;
      else hi = mid;
    }
    const g0 = lo;
    hi = n;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      const k = this.dir === -1 ? n - 1 - mid : mid;
      if (r[2 * k] <= e) lo = mid + 1;
      else hi = mid;
    }
    const g1 = lo - 1;
    if (g1 < g0) {
      this.kb = this.ka - 1;
      return;
    }
    const a = this.dir === -1 ? n - 1 - g1 : g0;
    const b = this.dir === -1 ? n - 1 - g0 : g1;
    this.ka = a > this.k0 ? a : this.k0;
    this.kb = b < this.k1 ? b : this.k1;
  }

  /**
   * Positions `[s, e)` in piece `k`: whether any of it is in view there, its x extents on screen and cut in `sx0`,
   * `sx1`, `scut` (the screen's sides: on a reversed piece the start of the range is on the right). What `spans` makes
   * a triple of, with no array: the layout's per-position loops call it directly.
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
    const ua = this.xOfUnit(this.unitOf(k, a));
    const ub = this.xOfUnit(this.unitOf(k, b));
    let xa = ua < ub ? ua : ub;
    let xb = ua < ub ? ub : ua;
    // a span that only touches the view's edge is not in view; a boundary on the edge is
    if (s === e ? xb < left || xa > right : xb <= left || xa >= right) return false;
    let cut = 0;
    // the range cut at its start and its end by the piece's edges, on the screen's sides; an end on a piece's edge at
    // a break is cut there too (a band projected across an intron is two ranges, square where the break cuts them)
    const cutStart = a > s || (a === ps && this.breakBelow(k));
    const cutEnd = b < e || (b === pe && this.breakAbove(k));
    if (this.dir === -1 ? cutEnd : cutStart) cut |= Cut.left;
    if (this.dir === -1 ? cutStart : cutEnd) cut |= Cut.right;
    if (xa < left) cut |= Cut.left;
    if (xb > right) cut |= Cut.right;
    if (xa < left) xa = left;
    if (xb > right) xb = right;
    this.sx0 = xa;
    this.sx1 = xb;
    this.scut = cut;
    return true;
  }

  /** Whether a break borders piece `k`'s lowest position (the piece below it in genomic order ends short of it). */
  private breakBelow(k: number): boolean {
    if (!this.sorted) return false;
    const r = this.pieces;
    const j = this.dir === -1 ? k + 1 : k - 1;
    return j >= 0 && 2 * j + 1 < r.length && r[2 * j + 1] !== r[2 * k];
  }

  /** Whether a break borders piece `k`'s highest edge (the piece above it in genomic order starts past it). */
  private breakAbove(k: number): boolean {
    if (!this.sorted) return false;
    const r = this.pieces;
    const j = this.dir === -1 ? k - 1 : k + 1;
    return j >= 0 && 2 * j < r.length && r[2 * j] !== r[2 * k + 1];
  }

  /**
   * The x extents of positions `[s, e)` as triples `x0, x1, cut` (unsnapped, clipped to the data box), a triple per
   * piece it touches in view, in display order; `cut` says which ends the view or a piece's edge cut. Appended to `out`.
   */
  spans(s: number, e: number, out: number[]): number[] {
    this.piecesIn(s, e);
    const ka = this.ka;
    const kb = this.kb;
    for (let k = ka; k <= kb; k++) if (this.spanIn(k, s, e)) out.push(this.sx0, this.sx1, this.scut);
    return out;
  }

  /** The positions in view, `[start, end, …]` per piece (whole cells: a partly shown one counts), appended to `out`. */
  visible(out: number[]): number[] {
    for (let k = this.k0; k <= this.k1; k++) {
      const ps = this.pieces[2 * k];
      const pe = this.pieces[2 * k + 1];
      const off = this.offsets[k];
      const a = Math.max(this.u0, off);
      const b = Math.min(this.u1, off + pe - ps);
      if (b <= a) continue;
      if (this.dir === -1) out.push(pe - Math.ceil(b - off), pe - Math.floor(a - off));
      else out.push(ps + Math.floor(a - off), ps + Math.ceil(b - off));
    }
    return out;
  }

  /** The axis's pieces as an axis shape's `ranges` (what the view's lookups read). */
  get ranges(): readonly number[] {
    return this.pieces;
  }

  /**
   * The breaks in view: `[x, hiddenStart, hiddenEnd, …]` triples (x unsnapped, in display order), appended to `out`. A
   * break on the view's edge counts.
   */
  breaks(out: number[]): number[] {
    const b = axisBreaks(this);
    const m = b.length / 3;
    if (!m) return out;
    // the first break at or after the window's start: a binary search over their offsets
    let lo = 0;
    let hi = m;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (b[3 * mid] < this.u0) lo = mid + 1;
      else hi = mid;
    }
    for (let i = lo; i < m && b[3 * i] <= this.u1; i++) out.push(this.xOfUnit(b[3 * i]), b[3 * i + 1], b[3 * i + 2]);
    return out;
  }
}
