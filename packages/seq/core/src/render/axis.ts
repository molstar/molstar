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
// pans with no scan of them). Groups (_spec/README.md "Render", Section layouts: text) put a gap after every `size`
// display units, counted from the axis's first: a cell skips the gaps, a span bridges them.

import { axisBreaks, monotonic, pieceAtOffset, pieceOffsets, pieceOfPosition } from '@molstar/seqx-core/view';

/**
 * A section's viewport (_spec/README.md "Render", Section layouts): `start` and `end` in display units (fractional),
 * read by the section's layout (tracks: the window; text: `start` the anchor, `end` derived; rings: the window), its
 * vertical scroll `top`, `scale` (text: px a position), and per band its own scroll and row zoom (`bands`, by track
 * address) and a grid's second axis (`y`).
 */
export interface Viewport {
  readonly start: number;
  readonly end: number;
  readonly top?: number;
  readonly scale?: number;
  readonly bands?: Readonly<Record<string, BandViewport>>;
  readonly y?: { readonly start: number; readonly end: number };
}

/** Which ends of a span the view or a piece's edge cut, on screen (they are drawn square). */
export const Cut = { left: 1, right: 2 } as const;

/** What a transform reads of its axis besides its pieces: its direction, and for its breaks, its system's. */
export interface AxisShapeOptions {
  readonly dir?: 1 | -1;
  readonly length?: number;
  readonly cyclic?: boolean;
  /** A gap of `gap` px after every `size` display units (text's groups of ten); none by default. */
  readonly groups?: { readonly size: number; readonly gap: number; readonly origin?: number };
}

/**
 * A band's own viewport (by track address): its rows' scroll (`bandTop`: the band's top in display rows below its
 * pinned band, fractional, as a step's view declares it; the layout turns it into px), row zoom (`rowScale`, px a row)
 * and the clades its rows layer folds (`folds`, by name: session state, compared by identity; the layer's `collapsed`
 * when absent).
 */
export interface BandViewport {
  readonly bandTop: number;
  readonly rowScale?: number;
  readonly folds?: readonly string[];
}

export class AxisTransform {
  /**
   * A transform on a resolved axis (its system, direction, length and cyclic flag) over `pieces` (its own, a block's,
   * a follow's image) between `x0` and `x0 + width`.
   */
  static on(
    axis: { readonly system: string; readonly dir: 1 | -1; readonly length?: number; readonly cyclic: boolean },
    pieces: readonly number[],
    x0: number,
    width: number,
    dpr: number,
    viewport?: Viewport,
  ): AxisTransform {
    return new AxisTransform(axis.system, pieces, x0, width, dpr, viewport, {
      dir: axis.dir,
      ...(axis.length !== undefined ? { length: axis.length } : {}),
      cyclic: axis.cyclic,
    });
  }

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
  /** Groups: their size in display units (0: none), the gap after each (px), the group the window starts in. */
  readonly gsize: number;
  readonly ggap: number;
  private readonly g0: number;
  /** The display unit groups are counted from (0; a text line's first unit when a piece edge starts the line). */
  private readonly gorigin: number;

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
    const g = shape.groups;
    this.gsize = g && g.size > 0 ? g.size : 0;
    this.ggap = this.gsize ? g!.gap : 0;
    this.gorigin = this.gsize ? (g!.origin ?? 0) : 0;
    this.g0 = this.gsize ? Math.floor((u0 - this.gorigin) / this.gsize) : 0;
    // the gaps inside the window take their px from the width
    const inner = this.gsize ? Math.max(0, Math.ceil((this.u1 - this.gorigin) / this.gsize) - 1 - this.g0) : 0;
    this.scale = Math.max(width - inner * this.ggap, 1e-9) / (this.u1 - this.u0);
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

  /** x of a display unit (a cell's left edge: past the gap before a group's first unit). */
  xOfUnit(u: number): number {
    if (this.gsize === 0) return this.x0 + (u - this.u0) * this.scale;
    return this.x0 + (u - this.u0) * this.scale + (Math.floor((u - this.gorigin) / this.gsize) - this.g0) * this.ggap;
  }

  /**
   * The x of units `va` to `vb` (either order) into `sx0` (`xOfUnit` of the lower) and `sx1` (`xOfEnd` of the higher;
   * the same x when equal). Inline, no number returned: once the hover and emphasis draws have met sections of pieces
   * and of groups, a double returned from a call that is not inlined is boxed.
   */
  private edges(va: number, vb: number): void {
    const lo = va < vb ? va : vb;
    const hi = va < vb ? vb : va;
    if (this.gsize === 0) {
      this.sx0 = this.x0 + (lo - this.u0) * this.scale;
      this.sx1 = this.x0 + (hi - this.u0) * this.scale;
      return;
    }
    this.sx0 =
      this.x0 + (lo - this.u0) * this.scale + (Math.floor((lo - this.gorigin) / this.gsize) - this.g0) * this.ggap;
    if (lo === hi) {
      this.sx1 = this.sx0;
      return;
    }
    const k = Math.ceil((hi - this.gorigin) / this.gsize) - 1 - this.g0;
    this.sx1 = this.x0 + (hi - this.u0) * this.scale + (k > 0 ? k : 0) * this.ggap;
  }

  /** x of a display unit as a right edge: before the gap after a group's last unit. */
  xOfEnd(u: number): number {
    if (this.gsize === 0) return this.x0 + (u - this.u0) * this.scale;
    const k = Math.ceil((u - this.gorigin) / this.gsize) - 1 - this.g0;
    return this.x0 + (u - this.u0) * this.scale + (k > 0 ? k : 0) * this.ggap;
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
    this.edges(this.unitOf(k, p), this.unitOf(k, p + 1));
    return Number.isFinite(this.sx0) && Number.isFinite(this.sx1);
  }

  /** The display unit under x (fractional), or −1 outside the data box; in a group's gap, the next group's first. */
  unitAt(x: number): number {
    if (x < this.x0 || x >= this.x0 + this.width) return -1;
    if (this.gsize === 0) return this.u0 + (x - this.x0) / this.scale;
    const d = x - this.x0;
    // px to the window's first group boundary, then whole groups with their gaps
    const first = (this.gorigin + (this.g0 + 1) * this.gsize - this.u0) * this.scale;
    if (d < first) return this.u0 + d / this.scale;
    const span = this.gsize * this.scale + this.ggap;
    const k = Math.floor((d - first) / span);
    const r = d - first - k * span - this.ggap;
    const u = this.gorigin + (this.g0 + 1 + k) * this.gsize;
    return r < 0 ? u : u + r / this.scale;
  }

  /**
   * Whether display unit `u` ends a group (its last unit): where a text's group label stands; `gapped` only when a gap
   * follows it in the window (an inline label stands in the gap).
   */
  endsGroup(u: number, gapped = false): boolean {
    if (this.gsize === 0) return false;
    const e = Math.floor(u) + 1;
    return (gapped ? e < this.u1 : e <= this.u1) && (e - this.gorigin) % this.gsize === 0;
  }

  /** Whether x falls in a group's gap (no position there). */
  inGap(x: number): boolean {
    if (this.gsize === 0 || x < this.x0 || x >= this.x0 + this.width) return false;
    const first = (this.gorigin + (this.g0 + 1) * this.gsize - this.u0) * this.scale;
    const d = x - this.x0 - first;
    if (d < 0) return false;
    const span = this.gsize * this.scale + this.ggap;
    return d - Math.floor(d / span) * span < this.ggap;
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
    // with groups the left edge past a gap, the right before one, so a span bridges the gaps inside it only
    this.edges(this.unitOf(k, a), this.unitOf(k, b));
    let xa = this.sx0;
    let xb = this.sx1;
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
