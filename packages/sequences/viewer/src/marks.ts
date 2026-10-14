/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the channels light (seqspec:interaction#channels, the marks; CP1a): per channel a `MarkSet`, per block of a
// section its sorted, disjoint block coordinate system runs, merged as they are added, so duplicates collapse (a
// nucleotide and its own codon light 77 once) and nothing caps them. A section's channels are its `SectionMarks`, which
// the overlay reads as `OverlayMarks`. Every buffer grows by doubling and is kept: rebuilding a set allocates nothing
// once it has room, and no query allocates.

import type { RunSink } from '@molstar/seq-core/format';
import type { ShownIndex } from './axis.js';

/**
 * A channel's bit; the pending selection draws from its own numbers. `row` flags a `RowMarks` rectangle of a row's
 * residues (their columns' hull, gaps included): its row is lit whole, header and data columns too.
 */
export const MarkBit = { selection: 1, highlight: 2, accent: 4, hover: 8, band: 16, row: 32 } as const;

/**
 * Block-system runs per block (section-local block indices): block b's runs are `[at[b], at[b + 1])` of `start` and
 * `stop`, sorted and disjoint once the block is closed. Filled block by block in ascending order: `reset`, then per
 * block `begin`, `add` (a `RunSink`: runs in any order), `close`. A block never begun lights nothing.
 */
export class MarkSet implements RunSink {
  at = new Int32Array(2);
  start = new Int32Array(16);
  stop = new Int32Array(16);
  /** Runs in all blocks. */
  count = 0;
  /** The last block begun; -1: none. */
  private open = -1;
  /** The open block's runs are in order and disjoint (adds kept them so). */
  private sorted = true;
  /** `has`'s last run: consecutive queries in order read it, not a search. */
  private cursor = 0;
  /**
   * Carets: boundaries lit (seqspec:interaction#channels), each a block and the index it precedes; block b's are
   * `[caretFrom(b), caretTo(b))`, sorted and distinct once the block is closed.
   */
  caretBlock = new Int32Array(4);
  caretAt = new Int32Array(4);
  carets = 0;
  private caretStart = new Int32Array(2);
  /** The open block's carets are in ascending order (boundaries kept them so). */
  private caretsSorted = true;

  /** Empty, for `blocks` blocks. */
  reset(blocks: number): void {
    if (this.at.length < blocks + 1) this.at = new Int32Array(Math.max(2 * this.at.length, blocks + 1));
    if (this.caretStart.length < blocks + 1)
      this.caretStart = new Int32Array(Math.max(2 * this.caretStart.length, blocks + 1));
    this.at[0] = 0;
    this.caretStart[0] = 0;
    this.count = 0;
    this.open = -1;
    this.cursor = 0;
    this.carets = 0;
  }

  /** Starts block b (after every block below it); the blocks skipped are empty. */
  begin(b: number): void {
    for (let k = this.open + 1; k <= b; k++) {
      this.at[k] = this.count;
      this.caretStart[k] = this.carets;
    }
    this.at[b + 1] = this.count;
    this.caretStart[b + 1] = this.carets;
    this.open = b;
    this.sorted = true;
    this.caretsSorted = true;
  }

  /** A run of the open block, merged into its last run when they overlap or touch, in either direction. */
  add(s: number, t: number): void {
    if (s >= t || this.open < 0) return;
    const first = this.at[this.open];
    const n = this.count;
    if (n > first) {
      const last = n - 1;
      if (s <= this.stop[last] && t >= this.start[last]) {
        if (t > this.stop[last]) this.stop[last] = t;
        if (s < this.start[last]) {
          this.start[last] = s;
          if (last > first && s <= this.stop[last - 1]) this.sorted = false;
        }
        return;
      }
      if (s < this.start[last]) this.sorted = false;
    }
    if (n === this.start.length) this.grow();
    this.start[n] = s;
    this.stop[n] = t;
    this.count = n + 1;
  }

  /**
   * A caret of the open block, before index `at` (a boundary at `hi`: after its last position), appended: a repeat of
   * the last skipped, the rest sorted and made distinct by `close`.
   */
  boundary(at: number): void {
    const b = this.open;
    if (b < 0) return;
    const n = this.carets;
    if (n > this.caretStart[b]) {
      if (this.caretAt[n - 1] === at) return;
      if (at < this.caretAt[n - 1]) this.caretsSorted = false;
    }
    if (this.carets === this.caretAt.length) {
      const [x, y] = [new Int32Array(2 * this.carets), new Int32Array(2 * this.carets)];
      x.set(this.caretBlock);
      y.set(this.caretAt);
      this.caretBlock = x;
      this.caretAt = y;
    }
    this.caretBlock[this.carets] = b;
    this.caretAt[this.carets++] = at;
  }

  /**
   * Ends the open block: its runs sorted only when out of order, merged, and clipped to `[lo, hi)`; its carets
   * likewise, distinct and clipped to `[lo, hi]`.
   */
  close(lo: number, hi: number): void {
    const b = this.open;
    if (b < 0) return;
    const c0 = this.caretStart[b];
    // The block field rides along (all b): one sort for both arrays.
    if (!this.caretsSorted) sortRuns(this.caretAt, this.caretBlock, c0, this.carets);
    let c = c0;
    for (let k = c0; k < this.carets; k++) {
      const at = this.caretAt[k];
      if (at >= lo && at <= hi && (c === c0 || this.caretAt[c - 1] !== at)) this.caretAt[c++] = at;
    }
    this.carets = c;
    this.caretStart[b + 1] = c;
    this.caretsSorted = true;
    const first = this.at[b];
    if (!this.sorted) sortRuns(this.start, this.stop, first, this.count);
    let w = first;
    for (let r = first; r < this.count; r++) {
      const s = this.start[r] < lo ? lo : this.start[r];
      const t = this.stop[r] > hi ? hi : this.stop[r];
      if (s >= t) continue;
      if (w > first && s <= this.stop[w - 1]) {
        if (t > this.stop[w - 1]) this.stop[w - 1] = t;
      } else {
        this.start[w] = s;
        this.stop[w] = t;
        w++;
      }
    }
    this.count = w;
    this.at[b + 1] = w;
    this.sorted = true;
  }

  /** Block b's runs as `[from, to)` indices into `start`/`stop`: `to - from`, its run count (0 for a block not begun). */
  from(b: number): number {
    return b < 0 || b > this.open ? 0 : this.at[b];
  }
  to(b: number): number {
    return b < 0 || b > this.open ? 0 : this.at[b + 1];
  }

  /** Block b's carets as `[caretFrom, caretTo)` indices into `caretAt` (none for a block not begun). */
  caretFrom(b: number): number {
    return b < 0 || b > this.open ? 0 : this.caretStart[b];
  }
  caretTo(b: number): number {
    return b < 0 || b > this.open ? 0 : this.caretStart[b + 1];
  }

  /** Whether block-system index i of block b is lit: the cursor's run or the next, else a binary search. */
  has(b: number, i: number): boolean {
    const lo = this.from(b);
    const hi = this.to(b);
    let c = this.cursor;
    if (c >= lo && c < hi && this.start[c] <= i) {
      if (i < this.stop[c]) return true;
      if (c + 1 < hi && i < this.start[c + 1]) return false;
      if (c + 1 < hi && i < this.stop[c + 1]) {
        this.cursor = c + 1;
        return true;
      }
    }
    c = this.after(lo, hi, i);
    if (c >= hi) return false;
    this.cursor = c;
    return this.start[c] <= i;
  }

  /** Whether any index of `[s, t)` on block b is lit. */
  any(b: number, s: number, t: number): boolean {
    const hi = this.to(b);
    const r = this.after(this.from(b), hi, s);
    return s < t && r < hi && this.start[r] < t;
  }

  /** Block b's first lit index; -1: none. */
  first(b: number): number {
    return this.to(b) > this.from(b) ? this.start[this.from(b)] : -1;
  }

  /** The first run in `[lo, hi)` that ends after i (`hi`: none). */
  private after(lo: number, hi: number, i: number): number {
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.stop[mid] <= i) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }

  private grow() {
    const n = 2 * this.start.length;
    const [s, t] = [new Int32Array(n), new Int32Array(n)];
    s.set(this.start);
    t.set(this.stop);
    this.start = s;
    this.stop = t;
  }
}

/** A heapsort of the runs `[a, b)` by start: in place, no comparator, no buffer. */
function sortRuns(start: Int32Array, stop: Int32Array, a: number, b: number) {
  const n = b - a;
  for (let i = (n >> 1) - 1; i >= 0; i--) sift(start, stop, a, i, n);
  for (let e = n - 1; e > 0; e--) {
    swap(start, stop, a, a + e);
    sift(start, stop, a, 0, e);
  }
}
function sift(start: Int32Array, stop: Int32Array, a: number, i: number, n: number) {
  for (;;) {
    let m = i;
    const l = 2 * i + 1;
    if (l < n && start[a + l] > start[a + m]) m = l;
    if (l + 1 < n && start[a + l + 1] > start[a + m]) m = l + 1;
    if (m === i) return;
    swap(start, stop, a + i, a + m);
    i = m;
  }
}
function swap(start: Int32Array, stop: Int32Array, i: number, j: number) {
  const s = start[i];
  const t = stop[i];
  start[i] = start[j];
  stop[i] = stop[j];
  start[j] = s;
  stop[j] = t;
}

/**
 * Rectangles in a section of an alignment's rows: rows `[rowStart, rowStop)` by block-system columns `[start, stop)` on
 * `block`, per channel `bit`. A rectangle extending the last one by the next row with the same columns joins it, so a
 * band over many rows stays one.
 */
export class RowMarks {
  count = 0;
  bit = new Uint8Array(8);
  block = new Int32Array(8);
  rowStart = new Int32Array(8);
  rowStop = new Int32Array(8);
  start = new Int32Array(8);
  stop = new Int32Array(8);

  /** Drops the rectangles of the bits in `bits` (a channel rebuilt), keeping the others in order. */
  clear(bits: number): void {
    let w = 0;
    for (let r = 0; r < this.count; r++) {
      if (this.bit[r] & bits) continue;
      if (w !== r) this.copy(r, w);
      w++;
    }
    this.count = w;
  }

  add(bit: number, block: number, row: number, s: number, t: number): void {
    const n = this.count;
    const l = n - 1;
    if (s >= t) return;
    const same = l >= 0 && this.bit[l] === bit && this.block[l] === block;
    if (same && this.rowStop[l] === row && this.start[l] === s && this.stop[l] === t) {
      this.rowStop[l] = row + 1;
      return;
    }
    if (n === this.bit.length) this.grow();
    this.set(n, bit, block, row, row + 1, s, t);
    this.count = n + 1;
  }

  private copy(r: number, w: number) {
    this.set(w, this.bit[r], this.block[r], this.rowStart[r], this.rowStop[r], this.start[r], this.stop[r]);
  }
  private set(n: number, bit: number, block: number, r0: number, r1: number, s: number, t: number) {
    this.bit[n] = bit;
    this.block[n] = block;
    this.rowStart[n] = r0;
    this.rowStop[n] = r1;
    this.start[n] = s;
    this.stop[n] = t;
  }

  private grow() {
    const n = 2 * this.bit.length;
    const bit = new Uint8Array(n);
    bit.set(this.bit);
    this.bit = bit;
    for (const f of ['block', 'rowStart', 'rowStop', 'start', 'stop'] as const) {
      const a = new Int32Array(n);
      a.set(this[f]);
      this[f] = a;
    }
  }
}

/** What the overlay and navigator passes read: reused, allocation-free. */
export interface OverlayMarks {
  /** `MarkBit` flags of a block index. */
  bitsAt(block: number, index: number): number;
  /** Whether any of `[start, stop)` is lit with one of `bits`. */
  any(block: number, start: number, stop: number, bits: number): boolean;
  /** One channel's runs (`MarkBit`: one bit); the navigator draws them. */
  runs(bit: number): MarkSet;
  /** The hover's first lit index the layout shows; -1: none. */
  primary(block: number): number;
  /** The hovered item's index when it is layer `layer`'s (the pairs overlay lights that pair alone); -1: none. */
  item(layer: number): number;
  /** That item's annotation (a layer drawn from two annotations has an item of one id in each); -1: none. */
  itemAnnotation(layer: number): number;
  /**
   * The hover's own hit on layer `layer`'s track (the pointer's or the keyboard's; never an external hover) into `out`:
   * its block, model index, state (a matrix band; -1) and part (-1), and where `out` has room a dot plot's partner's
   * unit, index and block (-1); false: none (seqspec:interaction#hover).
   */
  element(layer: number, out: Int32Array): boolean;
}

/** The axis reads the primary needs (`Axis`, `AxisView`). */
export interface ShownRuns {
  runToShown(block: number, start: number, stop: number, out: Int32Array): number;
  indexOf(block: number, k: ShownIndex): number;
}

/**
 * A section's marks: one `MarkSet` per channel bit, the rows' rectangles and the hover's primary per block. A set is
 * rebuilt only when its own channel or the view changes (the `Projector` fills it), never on a layout change.
 */
export class SectionMarks implements OverlayMarks {
  readonly selection = new MarkSet();
  readonly highlight = new MarkSet();
  readonly accent = new MarkSet();
  readonly hover = new MarkSet();
  readonly band = new MarkSet();
  readonly rows = new RowMarks();
  private primaries = new Int32Array(1);
  /** The hovered item: its layer, annotation and index (-1: none), set with the hover's rebuild. */
  hoverLayer = -1;
  hoverAnnotation = -1;
  hoverItem = -1;
  /** The hover's own hit (`element`): its track, block, index, state and part, and the view's track per layer. */
  hoverTrack = -1;
  hoverBlock = -1;
  hoverIndex = -1;
  hoverState = -1;
  hoverPart = -1;
  hoverPairUnit = -1;
  hoverPairIndex = -1;
  hoverPairBlock = -1;
  layerTrack: Int32Array | null = null;
  /**
   * Split rows on units of their own: those units, none a block of the section shows, and per unit the channels less
   * the others' runs, which a row on it draws (a run on one lights its rows alone).
   */
  ownUnits = new Int32Array(0);
  readonly own: SectionMarks[] = [];
  private readonly shown = new Int32Array(4);

  /** The marks a split row on unit u draws; null: this section's (no row of its own there). */
  ownOf(u: number): SectionMarks | null {
    for (let k = 0; u >= 0 && k < this.ownUnits.length; k++) if (this.ownUnits[k] === u) return this.own[k];
    return null;
  }

  bitsAt(block: number, index: number): number {
    let bits = 0;
    if (this.selection.has(block, index)) bits |= MarkBit.selection;
    if (this.highlight.has(block, index)) bits |= MarkBit.highlight;
    if (this.accent.has(block, index)) bits |= MarkBit.accent;
    if (this.hover.has(block, index)) bits |= MarkBit.hover;
    if (this.band.has(block, index)) bits |= MarkBit.band;
    return bits;
  }

  any(block: number, start: number, stop: number, bits: number): boolean {
    for (let bit = 1; bit <= MarkBit.band; bit <<= 1)
      if (bits & bit && this.runs(bit).any(block, start, stop)) return true;
    return false;
  }

  runs(bit: number): MarkSet {
    if (bit === MarkBit.selection) return this.selection;
    if (bit === MarkBit.highlight) return this.highlight;
    if (bit === MarkBit.accent) return this.accent;
    if (bit === MarkBit.band) return this.band;
    return this.hover;
  }

  item(layer: number): number {
    return layer >= 0 && layer === this.hoverLayer ? this.hoverItem : -1;
  }
  itemAnnotation(layer: number): number {
    return layer >= 0 && layer === this.hoverLayer ? this.hoverAnnotation : -1;
  }
  element(layer: number, out: Int32Array): boolean {
    const t = this.layerTrack;
    if (this.hoverTrack < 0 || this.hoverBlock < 0 || t === null || layer < 0 || layer >= t.length) return false;
    if (t[layer] !== this.hoverTrack) return false;
    out[0] = this.hoverBlock;
    out[1] = this.hoverIndex;
    out[2] = this.hoverState;
    out[3] = this.hoverPart;
    if (out.length < 6) return true;
    out[4] = this.hoverPairUnit;
    out[5] = this.hoverPairIndex;
    if (out.length > 6) out[6] = this.hoverPairBlock;
    return true;
  }

  primary(block: number): number {
    return block >= 0 && block < this.primaries.length ? this.primaries[block] : -1;
  }

  /**
   * The hover's primary per block: its first lit index that the layout shows, the lowest shown
   * index of its runs, so a filtered-out position gives no tick. After the hover set or the axis changes.
   */
  updatePrimary(axis: ShownRuns, blocks: number): void {
    if (this.primaries.length < blocks) this.primaries = new Int32Array(blocks);
    this.primaries.fill(-1);
    const h = this.hover;
    for (let b = 0; b < blocks; b++) {
      let k = 0x7fffffff;
      for (let r = h.from(b); r < h.to(b); r++)
        if (axis.runToShown(b, h.start[r], h.stop[r], this.shown) > 0 && this.shown[0] < k) k = this.shown[0];
      if (k < 0x7fffffff) this.primaries[b] = axis.indexOf(b, k as ShownIndex);
    }
  }
}
