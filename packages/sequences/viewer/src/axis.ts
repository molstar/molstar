/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A section's axis (seqspec:model#pieces-and-layout-maps): per block its piece (`first`, `dir`, the layout map), shown
// count, span, offset and key, side by side or stacked, and the block gap. Every index names its coordinate system: `k`
// the layout's shown index (`ShownIndex`), `i` the block coordinate system's. The run converters, `runToShown` and
// `shownToRuns`, are the only code that reads `dir` and a cyclic origin; `indexOf`, `shownOf` and `breakAt` read the
// same per-block segment table, built once from `blocks.layout` and the filter's `positions`. `AxisView` is the axis
// under the current viewport, one per section, which hits and draws share.

import { Arrangement, LigandForm, type ResolvedView } from '@molstar/seq-core/format';
import { type Brand, OrderedSet, type SeqSet, type SortedArray } from '@molstar/seq-core/model';
import type { RowModel } from './layout.js';
import { type PreparedRun, runXInto } from './runs.js';
import type { StreamView } from './streams.js';
import type { ViewportState } from './viewport.js';

/** An index in a section's layout coordinate system: a block's shown position, from 0. */
export type ShownIndex = Brand<number, 'shown'>;

export const clamp = (v: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, v));
export const sameList = <T>(a: readonly T[], b: readonly T[]) => a.length === b.length && a.every((x, i) => x === b[i]);

/** The scratch a block's segment walk fills: at most two segments (a cyclic origin), in shown order. */
const SEGS = 2;

/** A section's axis as the viewport reads it: per block its piece, shown count, span and key. Built at edge rate. */
export class Axis {
  readonly set: SeqSet;
  readonly count: number;
  readonly stacked: boolean;
  readonly unit: Int32Array;
  readonly alignment: Int32Array;
  readonly members: readonly (SortedArray | null)[];
  /** Per block: the filter's block indices; null: the whole piece. */
  readonly positions: readonly (OrderedSet | null)[];
  /** Per block: the domain's start (per-position arrays index `i - base`), the index at shown 0, the direction. */
  readonly base: Int32Array;
  /** Per block: its coordinate system's domain length (per-position arrays and pyramids span it). */
  readonly length: Int32Array;
  readonly first: Int32Array;
  readonly dir: Int8Array;
  readonly shown: Int32Array;
  readonly span: Float64Array;
  /** Cells before each block, side by side. */
  readonly offset: Float64Array;
  /** Side by side: every block's; stacked: the widest. */
  readonly cells: number;
  /** The unit's or alignment's key: committed blocks name them. */
  readonly names: readonly string[];
  /** What a block shows by key (names, members, piece, span), so equal axes compare across sets and contexts. */
  readonly keys: readonly string[];
  /**
   * Px between side-by-side blocks: `blockGap` up to 8 blocks, then shrinking so the gaps together take 7 of them (at
   * least 4 px each, never more than `blockGap`; 24 chromosomes: 9 px); 0 stacked, where every block starts at 0.
   */
  readonly gap: number;
  /** Set by the layout (LY1): the data area's left, and the surface px outside it across (header, columns, gutters, vnav). */
  left = 0;
  insetX = 0;
  /**
   * Set by the layout: the data area's left with the first k of the section's data columns shown, at `lefts[k]`, and
   * the px right of the data area at most (the gutter and the overview); `fitWidth` picks k for a width.
   */
  lefts: Float64Array | null = null;
  reserve = 0;
  /** Set by the layout: the least data area the columns leave (the section's `dataFloor`); 0: they never yield. */
  floor = 0;
  /** The data columns shown at the last `fitWidth`. */
  columnsShown = 0;
  /** Wrapped: the cell, fixed (no zoom); 0: linear. */
  fixedCell = 0;
  /** Set by the layout (LY1): the entries pinned above the scrolled ones, and the scrolled entries over [y0, y1) into [first, last]. */
  pinnedCount = 0;
  entriesIn: ((range: Float64Array) => void) | null = null;
  /** Builds the layout if stale: a vertical move needs the rows' extent (set by the layout manager). */
  ensure: (() => void) | null = null;
  /** The section's rows layer's display rows (set by the layout; a tree's header hit reads them). */
  rows: RowModel | null = null;
  /**
   * Per layer of the view: where the nearest track of residues (`letters`, `rows`) of its section is, 1 above its
   * track, -1 below, 0 none (a pair layer's `side: null` faces them; seqspec:format#layer-types).
   */
  readonly facing: Int8Array;
  /**
   * Per layer of the view: 1 when a `letters` layer of its own track is drawn under it (before it in the track), so a
   * mark it draws between positions keeps off their text (`ranges`' boundaries); else 0.
   */
  readonly lettered: Uint8Array;
  /**
   * Per block, `SEGS` slots of its layout map's segments in shown order: the lowest block index and length, the
   * direction, the first shown index and shown count, and the filter's ranks [lb, ub) in `positions` (unfiltered:
   * [0, length)). `segs` counts the used slots.
   */
  private readonly segs: Uint8Array;
  private readonly sI0: Int32Array;
  private readonly sLen: Int32Array;
  private readonly sDir: Int8Array;
  private readonly sK0: Int32Array;
  private readonly sN: Int32Array;
  private readonly sLb: Int32Array;
  private readonly sUb: Int32Array;
  /** CSR per block: the shown indices that follow a break, ascending. */
  private readonly breakStart: Int32Array;
  private readonly breaks: Int32Array;

  constructor(view: ResolvedView, set: SeqSet, s: number, gap: number, ligandCells: number) {
    const bv = view.blocks;
    const b0 = view.sections.blockStart[s];
    const n = view.sections.blockCount[s];
    const ligand = view.sections.ligandForm[s] === LigandForm.position;
    const key = (u: number) => set.strings.get(set.units.key[u]);
    this.set = set;
    this.count = n;
    this.stacked = view.sections.arrangement[s] === Arrangement.stacked;
    this.unit = bv.unit.slice(b0, b0 + n);
    this.alignment = bv.alignment.slice(b0, b0 + n);
    this.members = bv.members.slice(b0, b0 + n);
    this.positions = bv.positions.slice(b0, b0 + n);
    this.base = bv.base.slice(b0, b0 + n);
    this.length = bv.length.slice(b0, b0 + n);
    this.first = bv.first.slice(b0, b0 + n);
    this.dir = bv.dir.slice(b0, b0 + n);
    this.shown = new Int32Array(n);
    this.span = new Float64Array(n);
    this.offset = new Float64Array(n);
    this.segs = new Uint8Array(n);
    const slots = () => new Int32Array(SEGS * n);
    [this.sI0, this.sLen, this.sK0, this.sN, this.sLb, this.sUb] = [
      slots(),
      slots(),
      slots(),
      slots(),
      slots(),
      slots(),
    ];
    this.sDir = new Int8Array(SEGS * n);
    const breaks: number[] = [];
    this.breakStart = new Int32Array(n + 1);
    const names: string[] = [];
    const keys: string[] = [];
    let cells = 0;
    for (let b = 0; b < n; b++) {
      this.segments(b, bv.layout[b0 + b], breaks);
      this.breakStart[b + 1] = breaks.length;
      this.span[b] = ligand && bv.ligand[b0 + b] ? ligandCells : 1;
      this.offset[b] = this.stacked ? 0 : cells;
      const w = this.shown[b] * this.span[b];
      cells = this.stacked ? Math.max(cells, w) : cells + w;
      const a = set.alignments[this.alignment[b]];
      names.push(a ? set.strings.get(a.key) : key(this.unit[b]));
      const m = this.members[b];
      const piece = `${bv.start[b0 + b]}|${bv.stop[b0 + b]}|${this.dir[b]}`;
      keys.push(`${names[b]}|${m ? Array.from(m, key) : ''}|${piece}|${this.span[b]}`);
    }
    this.breaks = Int32Array.from(breaks);
    this.cells = cells;
    this.names = names;
    const L = view.layers;
    this.facing = new Int8Array(L.count);
    this.lettered = new Uint8Array(L.count);
    // The section's tracks of residues (`letters`, `rows`), once: a template's thousand instances are not squared.
    const residues: number[] = [];
    for (let j = 0; j < L.count; j++) {
      const type = L.type[j] >= 0 ? set.strings.get(L.type[j]) : '';
      if (view.tracks.section[L.track[j]] === s && (type === 'letters' || type === 'rows')) residues.push(L.track[j]);
      // Layers of a track are consecutive: after a letters layer, its track's later layers are lettered.
      if (j > 0 && L.track[j - 1] === L.track[j])
        this.lettered[j] =
          this.lettered[j - 1] === 1 || (L.type[j - 1] >= 0 && set.strings.get(L.type[j - 1]) === 'letters') ? 1 : 0;
    }
    for (let k = 0; k < L.count; k++) {
      const t = L.track[k];
      if (view.tracks.section[t] !== s) continue;
      let best = Number.POSITIVE_INFINITY;
      for (const u of residues) {
        if (u === t) continue;
        // The nearer wins; at equal distance, the residues above.
        const d = Math.abs(u - t) * 2 + (u < t ? 0 : 1);
        if (d < best) {
          best = d;
          this.facing[k] = u < t ? 1 : -1;
        }
      }
    }
    this.keys = keys;
    this.gap = this.stacked ? 0 : n > 8 ? Math.min(gap, Math.max(4, Math.round((gap * 7) / (n - 1)))) : gap;
  }

  /**
   * Block b's segment slots from its layout map (block coordinate system → layout, segments by `toLo`), the filter's
   * ranks within each, its shown count, and its breaks: within a segment where consecutive shown positions are not
   * consecutive block indices, and at a cyclic origin unless both sides reach it. Edge rate.
   */
  private segments(
    b: number,
    m: { count: number; byTo: Int32Array; at: Int32Array; length: Int32Array; dir: Int8Array },
    breaks: number[],
  ) {
    if (m.count > SEGS) throw new Error(`Axis: a layout map of ${m.count} segments (at most ${SEGS})`);
    const p = this.positions[b];
    let k0 = 0;
    for (let q = 0; q < m.count; q++) {
      const g = m.byTo[q];
      const s = SEGS * b + q;
      const [i0, len] = [m.at[g], m.length[g]];
      const lb = p === null ? 0 : OrderedSet.lowerBound(p, i0);
      const ub = p === null ? len : OrderedSet.lowerBound(p, i0 + len);
      [this.sI0[s], this.sLen[s], this.sDir[s], this.sK0[s], this.sLb[s], this.sUb[s]] = [
        i0,
        len,
        m.dir[g],
        k0,
        lb,
        ub,
      ];
      this.sN[s] = ub - lb;
      if (p !== null && q > 0 && this.sN[s] > 0 && this.sN[s - 1] > 0) {
        const into = (t: number, end: boolean) =>
          this.at(b, t, end === this.sDir[t] > 0 ? this.sUb[t] - 1 : this.sLb[t]);
        const reach = (t: number, end: boolean) =>
          into(t, end) === (end === this.sDir[t] > 0 ? this.sI0[t] + this.sLen[t] - 1 : this.sI0[t]);
        if (!reach(s - 1, true) || !reach(s, false)) breaks.push(k0);
      }
      for (let r = lb + 1; p !== null && typeof p !== 'number' && r < ub; r++) {
        if (p[r] !== p[r - 1] + 1) breaks.push(k0 + (m.dir[g] > 0 ? r - lb : ub - r));
      }
      k0 += this.sN[s];
    }
    this.segs[b] = m.count;
    this.shown[b] = k0;
    if (this.dir[b] < 0) {
      // No spread: a filter's breaks can pass V8's argument limit (150,000 runs, VL-02).
      const tail = breaks.splice(this.breakStart[b]).sort((x, y) => x - y);
      for (const x of tail) breaks.push(x);
    }
  }

  /** The block index at rank r of segment s (unfiltered: its offset). */
  private at(b: number, s: number, r: number): number {
    const p = this.positions[b];
    return p === null ? this.sI0[s] + r : OrderedSet.getAt(p, r);
  }

  /**
   * The data columns that fit `width`: they yield from the last while the data area would be narrower than `floor`
   * (none yield at 0); the header never does. Sets `left` and `insetX` (keeping the px right of the data area), and
   * `columnsShown`. The viewport calls it with a new size (before it decides the data area is empty), the layout before
   * a build, so both read the same left.
   */
  fitWidth(width: number): void {
    const lefts = this.lefts;
    if (!lefts) return;
    let k = lefts.length - 1;
    while (k > 0 && this.floor > 0 && width - lefts[k] - this.reserve < this.floor) k--;
    const right = this.insetX - this.left;
    this.columnsShown = k;
    this.left = lefts[k];
    this.insetX = this.left + right;
  }

  /** The same geometry besides the blocks: a changed gap or data area refreshes a kept viewport. */
  sameFrame(o: Axis): boolean {
    return o.gap === this.gap && o.left === this.left && o.insetX === this.insetX && o.fixedCell === this.fixedCell;
  }

  /** Shown index → block index; -1 outside the block. */
  indexOf(b: number, k: ShownIndex): number {
    let s = SEGS * b;
    if (this.segs[b] > 1 && k >= this.sK0[s + 1]) s++;
    const o = k - this.sK0[s];
    if (o < 0 || o >= this.sN[s]) return -1;
    return this.at(b, s, this.sDir[s] > 0 ? this.sLb[s] + o : this.sUb[s] - 1 - o);
  }
  /** Block index → shown index; -1: outside the piece or filtered out. */
  shownOf(b: number, i: number): ShownIndex {
    const p = this.positions[b];
    for (let s = SEGS * b, e = s + this.segs[b]; s < e; s++) {
      if (i < this.sI0[s] || i >= this.sI0[s] + this.sLen[s]) continue;
      const r = p === null ? i - this.sI0[s] : OrderedSet.indexOf(p, i);
      if (r < 0) return -1 as ShownIndex;
      return (this.sK0[s] + (this.sDir[s] > 0 ? r - this.sLb[s] : this.sUb[s] - 1 - r)) as ShownIndex;
    }
    return -1 as ShownIndex;
  }
  /** Block index → the shown index of the nearest shown position (the index itself when shown); -1: none shown. */
  nearestShown(b: number, i: number): ShownIndex {
    const p = this.positions[b];
    let best = -1;
    let gap = 0x7fffffff;
    for (let s = SEGS * b, e = s + this.segs[b]; s < e; s++) {
      const lb = this.sLb[s];
      const ub = this.sUb[s] - 1;
      if (ub < lb) continue;
      // The segment's kept rank nearest to i, by block index.
      let r = clamp(p === null ? i - this.sI0[s] : OrderedSet.lowerBound(p, i), lb, ub);
      if (p !== null && r > lb && OrderedSet.getAt(p, r) - i > i - OrderedSet.getAt(p, r - 1)) r--;
      const d = Math.abs((p === null ? this.sI0[s] + r : OrderedSet.getAt(p, r)) - i);
      if (d >= gap) continue;
      gap = d;
      best = this.sK0[s] + (this.sDir[s] > 0 ? r - lb : ub - r);
    }
    return best as ShownIndex;
  }
  /** A break before shown k: k follows a position the block coordinate system does not. */
  breakAt(b: number, k: ShownIndex): boolean {
    let lo = this.breakStart[b];
    let hi = this.breakStart[b + 1];
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.breaks[mid] < k) lo = mid + 1;
      else hi = mid;
    }
    return lo < this.breakStart[b + 1] && this.breaks[lo] === k;
  }
  /** The shown index where block b's second segment starts (a cyclic origin); -1: one segment. */
  wrapAt(b: number): number {
    return this.segs[b] > 1 ? this.sK0[SEGS * b + 1] : -1;
  }

  /** The filter's rank of the first position at or after block index i in segment s, within its ranks. */
  private rank(b: number, s: number, i: number): number {
    const p = this.positions[b];
    if (p === null) return clamp(i - this.sI0[s], 0, this.sLen[s]);
    return clamp(OrderedSet.lowerBound(p, i), this.sLb[s], this.sUb[s]);
  }
  /**
   * A block-system run [start, stop) → its shown runs into `out` (pairs), ascending and merged where they touch: at
   * most two (across a cyclic origin). Their count. Allocation-free.
   */
  runToShown(b: number, start: number, stop: number, out: Int32Array): number {
    let n = 0;
    for (let s = SEGS * b, e = s + this.segs[b]; s < e; s++) {
      const ra = this.rank(b, s, Math.max(start, this.sI0[s]));
      const rc = this.rank(b, s, Math.min(stop, this.sI0[s] + this.sLen[s]));
      if (rc <= ra) continue;
      const k0 = this.sK0[s] + (this.sDir[s] > 0 ? ra - this.sLb[s] : this.sUb[s] - rc);
      const k1 = k0 + rc - ra;
      if (n > 0 && out[2 * n - 1] === k0) out[2 * n - 1] = k1;
      else {
        out[2 * n] = k0;
        out[2 * n + 1] = k1;
        n++;
      }
    }
    return n;
  }
  /**
   * A shown run [k0, k1) → its block-system runs into `out` (pairs), in shown order (a reversed piece's descend): one
   * per filter run, split at breaks. Their count; past `out`'s room the rest are counted, not written. Allocation-free.
   */
  shownToRuns(b: number, k0: number, k1: number, out: Int32Array): number {
    let n = 0;
    const room = out.length >> 1;
    for (let s = SEGS * b, e = s + this.segs[b]; s < e; s++) {
      const lo = Math.max(k0, this.sK0[s]);
      const hi = Math.min(k1, this.sK0[s] + this.sN[s]);
      // Runs between breaks: a break at x starts a run there.
      let x = lo;
      let j = this.breakStart[b];
      let jh = this.breakStart[b + 1];
      while (j < jh) {
        const mid = (j + jh) >> 1;
        if (this.breaks[mid] <= lo) j = mid + 1;
        else jh = mid;
      }
      while (x < hi) {
        const y = j < this.breakStart[b + 1] && this.breaks[j] < hi ? this.breaks[j++] : hi;
        if (n < room) {
          const a = this.indexOf(b, x as ShownIndex);
          const c = this.indexOf(b, (y - 1) as ShownIndex);
          out[2 * n] = Math.min(a, c);
          out[2 * n + 1] = Math.max(a, c) + 1;
        }
        n++;
        x = y;
      }
    }
    return n;
  }

  /**
   * Equal axes (seqspec:interaction#viewports): the same blocks by key (pieces included) and the same shown positions.
   */
  equals(o: Axis): boolean {
    if (o.count !== this.count || o.stacked !== this.stacked || !sameList(o.keys, this.keys)) return false;
    return this.positions.every((p, b) => {
      const q = o.positions[b];
      return p === q || (p !== null && q !== null && OrderedSet.areEqual(p, q));
    });
  }
}

/** `AxisView.geometry`. */
export const AxisGeometry = { linear: 0, wrapped: 1, circular: 2 } as const;

/** The gap a cell leaves on its right (wide cells: one device px), so cells and the marks over them share edges. */
export const cellGap = (width: number, dpr: number) => (width > 6 ? 1 / dpr : 0);
/**
 * The x extent of shown cells [k0, k1) of a block as drawn (linear): from the first cell's edge to the last cell's
 * right, less its gap. Marks (hover, selection) cover exactly the cells' pixels.
 */
export function cellSpan(v: AxisView, b: number, k0: ShownIndex, k1: ShownIndex, out: Float64Array) {
  out[0] = v.x(b, k0);
  out[1] = v.x(b, k1) - cellGap(v.cell * v.span[b], v.dpr);
}

/**
 * A point on the navigator's map (`io[0]`, px) as the axis's units, the inverse of the window's map: into `io[0]` the
 * cells before it (spans counted, as `Axis.offset`) and `io[1]` the block gaps, so its content x at a cell c is `io[0]
 * * c + io[1] * axis.gap` (doubles in and out through `io`).
 */
export function navUnits(v: AxisView, io: Float64Array): void {
  const ax = v.axis;
  const x = io[0];
  io[0] = io[1] = 0;
  if (!ax || ax.count === 0 || v.navCell <= 0) return;
  let b = 0;
  while (!ax.stacked && b + 1 < ax.count && v.navX0s[b + 1] <= x) b++;
  const len = ax.shown[b] * ax.span[b];
  const t = Math.max(0, (x - v.navX0s[b]) / v.navCell);
  io[0] = (ax.stacked ? 0 : ax.offset[b]) + Math.min(t, ax.stacked ? ax.cells : len);
  if (ax.stacked || t <= len) io[1] = ax.stacked ? 0 : b;
  else io[1] = b + (b + 1 < ax.count && v.navGap > 0 ? Math.min(1, (x - v.navX1s[b]) / v.navGap) : 0);
}

/** A section's headless geometry under its viewport: hits use it, the engine's views extend it. */
export interface AxisView {
  /** `AxisGeometry`. */
  geometry: number;
  /** Per block: the index at shown position 0; 1, or -1 for a reversed piece. */
  first: Int32Array;
  dir: Int8Array;
  /** Per block: x of shown index 0 (linear); cells per position. */
  blockX: Float64Array;
  span: Float64Array;
  /** Per block: the filter's block indices; null: all. */
  positions: readonly (OrderedSet | null)[];
  /** Per block, shown indices [lo, hi). */
  visLo: Int32Array;
  visHi: Int32Array;
  cell: number;
  lod: number;
  reductions: boolean;
  dpr: number;
  top: number;
  height: number;
  left: number;
  right: number;
  gutterLeft: number;
  gutterRight: number;
  /** Stacked: the lane's block (wrapped: the line's); a multi-row layer's row. */
  onlyBlock: number;
  row: number;
  lane: number;
  /** Wrapped: the line, its shown range, its x from the data area's left (a chip's on its row) and the letters' size. */
  line: number;
  lineStart: number;
  lineStop: number;
  lineLeft: number;
  /**
   * Wrapped: the entry's track is drawn in full (a height of its own; seqspec:viewer#the-text-view): its layers draw
   * their track forms in the band, the line standing in for a linear view (`visLo`/`visHi` the line's positions, `cell`
   * its advance), not their lanes.
   */
  trackForm: boolean;
  fontPx: number;
  /**
   * Wrapped: an annotation lane's height under the text, a strip lane's, and the gap above a track's first
   * (`wrappedLaneHeight`, `wrappedStripHeight`, `wrappedTrackGap`).
   */
  laneHeight: number;
  stripHeight: number;
  trackGap: number;
  /** The visible rows region. */
  regionTop: number;
  regionBottom: number;
  scrollY: number;
  /** The navigator's window across the data area, and the vertical overview's down its column (gestures hit-test them). */
  navX0: number;
  navX1: number;
  /** The navigator's map of the whole axis from `left`: the data area's width, or the axis's when narrower. */
  navWidth: number;
  /**
   * The map as at the fit: per block its extent on it, and its px per cell and per block gap; the window is mapped
   * through the blocks, so neither moves with zoom.
   */
  navX0s: Float64Array;
  navX1s: Float64Array;
  navCell: number;
  navGap: number;
  vnavY0: number;
  vnavY1: number;
  /** The vertical overview's column across (0: none); its map spans the scrolled region (`regionTop`, `regionBottom`). */
  vnavWidth: number;
  /** The section's axis (its set, per-block tables). */
  axis: Axis | null;
  /**
   * The layer drawn or hit: its regions held by view (seqspec:sources#the-viewers-side, `ctx.streams.get(layer)`);
   * null: none.
   */
  stream: StreamView | null;
  /** The layer drawn or hit: where its residues are (`Axis.facing`). */
  facing: number;
  /** Device-snapped. */
  x(block: number, shown: ShownIndex): number;
  shownOf(block: number, index: number): ShownIndex;
  indexOf(block: number, shown: ShownIndex): number;
  breakAt(block: number, k: ShownIndex): boolean;
  runToShown(block: number, start: number, stop: number, out: Int32Array): number;
  shownToRuns(block: number, k0: number, k1: number, out: Int32Array): number;
}

/** What `SectionAxisView.sync` reads of the section's layout (LY1's `SectionLayout` has these). */
export interface AxisFrame {
  readonly geometry: number;
  readonly dataLeft: number;
  readonly dataRight: number;
  readonly scrolledTop: number;
  readonly height: number;
  readonly vnavWidth?: number;
  /** The circular geometry: the axis unrolled, CSS px (`SectionLayout.ringW`). */
  readonly ringW?: number;
}

const NONE = new Int32Array(0);
const NO_DIR = new Int8Array(0);
/** One per section, owned by the layout manager; `sync` brings it to the current layout and viewport. */
export class SectionAxisView implements AxisView {
  geometry: number = AxisGeometry.linear;
  first: Int32Array = NONE;
  dir: Int8Array = NO_DIR;
  blockX: Float64Array = new Float64Array(0);
  /**
   * Per block: x of shown index 0 before the scroll, and the scroll in whole device px. `x` rounds the unscrolled edge
   * and then subtracts the scroll, so a pan moves every edge by exactly its device px (rounding `blockX + k·cell`
   * flipped edges near .5 by floating-point error, and a copied column then differed from a full redraw's).
   */
  blockX0: Float64Array = new Float64Array(0);
  scrollDev = 0;
  span: Float64Array = new Float64Array(0);
  positions: readonly (OrderedSet | null)[] = [];
  visLo: Int32Array = NONE;
  visHi: Int32Array = NONE;
  cell = 1;
  lod = 0;
  reductions = true;
  dpr = 1;
  top = 0;
  height = 0;
  left = 0;
  right = 0;
  gutterLeft = 0;
  gutterRight = 0;
  onlyBlock = -1;
  row = -1;
  lane = -1;
  line = -1;
  lineStart = 0;
  lineStop = 0;
  lineLeft = 0;
  fontPx = 0;
  laneHeight = 5;
  stripHeight = 3;
  trackGap = 4;
  regionTop = 0;
  regionBottom = 0;
  scrollY = 0;
  navX0 = 0;
  navX1 = 0;
  navWidth = 0;
  navX0s = new Float64Array(0);
  navX1s = new Float64Array(0);
  navCell = 0;
  navGap = 0;
  vnavY0 = 0;
  vnavY1 = 0;
  vnavWidth = 0;
  axis: Axis | null = null;
  stream: StreamView | null = null;
  facing = 0;
  /** Wrapped: each block's prepared run (the layout's), and px between groups. */
  runs: readonly PreparedRun[] = [];
  trackForm = false;
  private lineLo = new Int32Array(0);
  private lineHi = new Int32Array(0);
  /** The view's own, while a line stands in for it (`setTrackForm`). */
  private keptLo: Int32Array = NONE;
  private keptHi: Int32Array = NONE;
  private keptCell = 1;
  /**
   * Wrapped: draws (and hits) the entry's line as a linear view where its track is drawn in full (`trackForm`), the
   * view's own range and cell kept and put back after. Set after the line (`onlyBlock`, `lineStart`, `lineStop`).
   * Allocation-free once grown to the section's block count.
   */
  setTrackForm(full: boolean): void {
    if (!full && !this.trackForm) return;
    if (!this.trackForm) {
      this.keptLo = this.visLo;
      this.keptHi = this.visHi;
      this.keptCell = this.cell;
    }
    if (!full) {
      this.visLo = this.keptLo;
      this.visHi = this.keptHi;
      this.cell = this.keptCell;
      this.trackForm = false;
      return;
    }
    const n = this.span.length;
    if (this.lineLo.length < n) {
      this.lineLo = new Int32Array(n);
      this.lineHi = new Int32Array(n);
    }
    this.lineLo.fill(0);
    this.lineHi.fill(0);
    const b = this.onlyBlock;
    if (b >= 0 && b < n) {
      this.lineLo[b] = this.lineStart;
      this.lineHi[b] = this.lineStop;
      const run = this.runs[b];
      if (run && this.span[b] > 0) this.cell = run.advance / this.span[b];
    }
    this.visLo = this.lineLo;
    this.visHi = this.lineHi;
    this.trackForm = true;
  }
  gap = 0;
  private readonly xs = new Float64Array(1);
  /** The layout and viewport versions it was synced at. */
  readonly at = new Float64Array(2);

  /** Allocation-free unless the block count changed. Wrapped x is the layout's runs' (`runs`, `lineStart`). */
  sync(ax: Axis, l: AxisFrame, s: ViewportState) {
    this.trackForm = false;
    this.axis = ax;
    this.geometry = l.geometry;
    if (this.blockX.length !== ax.count) {
      this.blockX = new Float64Array(ax.count);
      this.blockX0 = new Float64Array(ax.count);
      this.navX0s = new Float64Array(ax.count);
      this.navX1s = new Float64Array(ax.count);
    }
    this.first = ax.first;
    this.dir = ax.dir;
    this.span = ax.span;
    this.positions = ax.positions;
    this.visLo = s.visLo;
    this.visHi = s.visHi;
    this.cell = s.cell;
    this.lod = s.lod;
    this.dpr = s.dpr;
    this.left = l.dataLeft;
    this.right = l.dataRight;
    this.regionTop = this.top = l.scrolledTop;
    this.regionBottom = l.height;
    this.height = l.height - l.scrolledTop;
    this.vnavWidth = l.vnavWidth ?? 0;
    this.scrollY = s.scrollY;
    this.scrollDev = Math.round(s.scrollX * s.dpr);
    for (let b = 0; b < ax.count; b++) {
      this.blockX0[b] = l.dataLeft + ax.offset[b] * s.cell + b * ax.gap;
      this.blockX[b] = this.blockX0[b] - s.scrollX;
    }
    if (l.geometry === AxisGeometry.circular) this.unrolled(ax, l.dataLeft, l.ringW ?? 0);
    const w = l.dataRight - l.dataLeft;
    const gaps = ax.count > 1 ? (ax.count - 1) * ax.gap : 0;
    // The navigator maps the whole axis as at the fit cell, never the empty data area past an axis narrower than it.
    const fit = s.minCell > 0 ? s.minCell : s.cell;
    const fw = ax.cells * fit + gaps;
    const scale = fw > w && fw > 0 ? w / fw : 1;
    this.navCell = fit * scale;
    this.navGap = ax.gap * scale;
    for (let b = 0; b < ax.count; b++) {
      this.navX0s[b] = l.dataLeft + (ax.stacked ? 0 : ax.offset[b] * this.navCell + b * this.navGap);
      this.navX1s[b] = this.navX0s[b] + ax.shown[b] * ax.span[b] * this.navCell;
    }
    this.navWidth = fw > 0 ? fw * scale : w;
    this.xs[0] = s.scrollX;
    this.mapToNav(ax);
    this.navX0 = this.xs[0];
    this.xs[0] = Math.min(s.scrollX + w, ax.cells * s.cell + gaps);
    this.mapToNav(ax);
    this.navX1 = this.xs[0];
    const rh = l.height - l.scrolledTop;
    const ch = s.contentHeight;
    this.vnavY0 = l.scrolledTop + (ch > 0 ? (s.scrollY / ch) * rh : 0);
    this.vnavY1 = l.scrolledTop + (ch > 0 ? Math.min(1, (s.scrollY + rh) / ch) * rh : rh);
  }

  /**
   * The circular geometry (seqspec:viewer#the-render-backend): the whole axis laid along the outer ring's circumference
   * `w` from `left`, block gaps included (a gap closes the ring too: as many gaps as blocks, but one block's), every
   * position in view, no scroll. Draws and hits then read x as along the ring; `Ring` maps it to an angle. Its own
   * `visLo`/`visHi`.
   */
  private unrolled(ax: Axis, left: number, w: number) {
    if (this.ringLo.length !== ax.count) {
      this.ringLo = new Int32Array(ax.count);
      this.ringHi = new Int32Array(ax.count);
    }
    const gaps = ax.count > 1 ? ax.count * ax.gap : 0;
    this.cell = ax.cells > 0 ? Math.max(0, w - gaps) / ax.cells : 1;
    this.scrollDev = 0;
    // viewport.ts's `Lod` rule at the section's `minLetterCell` (letters, cells, thin, aggregate; VL-17).
    const c = this.cell;
    this.lod = c >= this.minLetterCell ? 0 : c >= 3 ? 1 : c >= 1 ? 2 : 3;
    for (let b = 0; b < ax.count; b++) {
      this.blockX[b] = this.blockX0[b] = left + ax.offset[b] * this.cell + b * (ax.count > 1 ? ax.gap : 0);
      this.ringHi[b] = ax.shown[b];
    }
    this.visLo = this.ringLo;
    this.visHi = this.ringHi;
    this.right = left + w;
  }
  private ringLo = new Int32Array(0);
  private ringHi = new Int32Array(0);
  /** The section's `minLetterCell` (render props), for the rings' LOD; set by the layout before `sync`. */
  minLetterCell = 8;

  /**
   * Content x (from the axis's start, at the view's cell) in `xs[0]` to its x on the navigator's map, in place: through
   * the block or the gap it falls in (a double crosses no call).
   */
  private mapToNav(ax: Axis): void {
    const cx = this.xs[0];
    const c = this.cell;
    if (ax.count === 0) {
      this.xs[0] = this.left;
      return;
    }
    let b = 0;
    while (!ax.stacked && b + 1 < ax.count && ax.offset[b + 1] * c + (b + 1) * ax.gap <= cx) b++;
    const t = cx - (ax.stacked ? 0 : ax.offset[b] * c + b * ax.gap);
    const len = ax.shown[b] * ax.span[b] * c;
    const inGap = t > len && ax.gap > 0 && !ax.stacked;
    this.xs[0] = inGap
      ? this.navX1s[b] + (Math.min(t - len, ax.gap) / ax.gap) * this.navGap
      : this.navX0s[b] + (c > 0 ? (t * this.navCell) / c : 0);
  }

  /**
   * Linear: shown k of block b (the floor of an unsnapped x) moved to the cell drawn on device px `pd`, cell k covering
   * [xd(k), xd(k + 1)) as `DrawView.xd` places it, so a hit or a drag names the cell under the pointer's device px
   * (GG-01; seqspec:interaction#hover). Ints in and out (no double crosses the call).
   */
  drawnShown(b: number, k: number, pd: number): number {
    const w = this.span[b] * this.cell;
    const x0 = this.blockX0[b];
    if (!(w > 0)) return k;
    // In closed form (LC-29: stepping one cell at a time cost 60,000 steps a hit at chr2's overview): the last cell
    // whose left edge rounds to `pd` or before, (x0 + j·w)·dpr < pd + scrollDev + ½, then exactly.
    const j = Math.ceil(((pd + this.scrollDev + 0.5) / this.dpr - x0) / w) - 1;
    if (j < 0x7fffffff) k = j < -1 ? -1 : j;
    while (k > 0 && Math.round((x0 + k * w) * this.dpr) - this.scrollDev > pd) k--;
    while (Math.round((x0 + (k + 1) * w) * this.dpr) - this.scrollDev <= pd) k++;
    return k;
  }
  /**
   * Linear: the first cell of block b counted on device px `pd` where `k` is the cell drawn there (`drawnShown`): below
   * a device px a cell, the cells whose left edge rounds onto `pd` (the ones a density column holds, `DrawView.xd`),
   * else `k`. Ints in and out.
   */
  drawnFirst(b: number, k: number, pd: number): number {
    const w = this.span[b] * this.cell;
    const x0 = this.blockX0[b];
    if (!(w * this.dpr < 1) || !(w > 0)) return k;
    // The first cell whose left edge rounds to `pd` or past it: (x0 + j·w)·dpr ≥ pd + scrollDev − ½, then exactly.
    let j = Math.max(0, Math.min(k, Math.ceil(((pd + this.scrollDev - 0.5) / this.dpr - x0) / w))) | 0;
    while (j > 0 && Math.round((x0 + (j - 1) * w) * this.dpr) - this.scrollDev >= pd) j--;
    while (j < k && Math.round((x0 + j * w) * this.dpr) - this.scrollDev < pd) j++;
    return j;
  }
  x(block: number, shown: ShownIndex): number {
    if (this.geometry === AxisGeometry.wrapped) {
      runXInto(this.runs[block], this.lineStart, shown, this.gap, this.xs, 0);
      return this.left + this.lineLeft + this.xs[0];
    }
    return (
      (Math.round((this.blockX0[block] + shown * this.span[block] * this.cell) * this.dpr) - this.scrollDev) / this.dpr
    );
  }
  shownOf(block: number, index: number): ShownIndex {
    return this.axis ? this.axis.shownOf(block, index) : (-1 as ShownIndex);
  }
  indexOf(block: number, shown: ShownIndex): number {
    return this.axis ? this.axis.indexOf(block, shown) : -1;
  }
  breakAt(block: number, k: ShownIndex): boolean {
    return this.axis ? this.axis.breakAt(block, k) : false;
  }
  runToShown(block: number, start: number, stop: number, out: Int32Array): number {
    return this.axis ? this.axis.runToShown(block, start, stop, out) : 0;
  }
  shownToRuns(block: number, k0: number, k1: number, out: Int32Array): number {
    return this.axis ? this.axis.shownToRuns(block, k0, k1, out) : 0;
  }
}
