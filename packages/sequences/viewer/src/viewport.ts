/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Viewports, one per section (seqspec:interaction#viewports), and scroll-group membership
// (seqspec:interaction#linked-sections). A viewport's `live` record is mutated per change and its version bumps only
// when something changed; `committed` emits when a change settles.
// The geometry is the section's axis (axis.ts): blocks side by side (or stacked) at `cell` px per shown position, in
// the data area (`width` less the axis's `insetX`); the layout gives the scrolled rows' extent through `setRows`.
// x arguments (`zoomAt`) are px from the data area's left. Block ranges (`zoomTo`, `committed.blocks`) are in the
// block coordinate system; the axis's run converters take them to shown positions.

import type { ResolvedView, SeqViewSpecDoc } from '@molstar/seq-core/format';
import { Extents, type SeqSet } from '@molstar/seq-core/model';
import { Axis, clamp, type ShownIndex, sameList } from './axis.js';
import { type ConfigManager, SeqViewerConfig } from './config.js';
import type { SeqLociAny } from './lineage.js';
import { BlockHulls, type Coordinates, Projector } from './projection.js';
import type { Scheduler } from './scheduler.js';
import { live, type Origin, pulled, type SeqState, state, type WritableLive } from './state.js';

/**
 * Every section of every view of a document, the layouts' and the viewports' lifecycle (seqspec:interaction#viewports).
 */
export function sectionsOf(doc: SeqViewSpecDoc | null): string[] {
  const out: string[] = [];
  for (const v of doc?.root.children ?? [])
    if (v.kind === 'view') for (const c of v.children ?? []) if (c.kind === 'section' && c.ref) out.push(c.ref);
  return out;
}

export const Lod = { letters: 0, cells: 1, thin: 2, aggregate: 3 } as const;

/** The render props a viewport reads (LY1's `RenderProps` has them; seqspec:viewer#render-props). */
export interface ViewportProps {
  readonly maxCell: number;
  readonly minLetterCell: number;
  readonly maxRowScale: number;
  readonly blockGap: number;
  readonly ligandCells: number;
}
const DefaultViewportProps: ViewportProps = {
  maxCell: 28,
  minLetterCell: 8,
  maxRowScale: 3,
  blockGap: 28,
  ligandCells: 3,
};

/** Reused; CSS px unless noted. */
export interface ViewportState {
  width: number;
  height: number;
  dpr: number;
  /** Px per shown position. */
  cell: number;
  minCell: number;
  maxCell: number;
  /** Snapped to device pixels. */
  scrollX: number;
  scrollY: number;
  rowScale: number;
  minRowScale: number;
  maxRowScale: number;
  contentWidth: number;
  contentHeight: number;
  /** `Lod`. */
  lod: number;
  /** Per block, shown indices [lo, hi); reallocated only when the block count changes. */
  visLo: Int32Array;
  visHi: Int32Array;
  firstEntry: number;
  lastEntry: number;
  animating: boolean;
  interacting: boolean;
  /** Of the last change. */
  origin: Origin;
}
export interface CommittedViewport {
  readonly origin: Origin;
  readonly cell: number;
  readonly scrollX: number;
  readonly scrollY: number;
  readonly rowScale: number;
  /** Visible blocks by unit or alignment key, block coordinate system indices; a piece across a cyclic origin twice. */
  readonly blocks: readonly { readonly ref: string; readonly start: number; readonly stop: number }[];
  readonly rows: readonly [number, number];
}
export interface MoveOpts {
  animate?: boolean;
  durationMs?: number;
  origin?: Origin;
  /** A fraction of the range added on each side. */
  pad?: number;
}
export interface SectionViewport {
  readonly ref: string;
  readonly live: ViewportLive;
  readonly committed: SeqState<CommittedViewport>;
  scrollTo(x: number, y: number, origin?: Origin): void;
  panBy(dx: number, dy: number, origin?: Origin): void;
  /** Anchored at x. */
  zoomAt(factor: number, x: number, origin?: Origin): void;
  /** Block coordinate system indices. */
  zoomTo(block: number, start: number, stop: number, opts?: MoveOpts): void;
  /**
   * Several blocks' ranges at once, (block, start, stop) triples (a focus's `blocks`): the view from the first one's
   * start to the last one's stop, at the largest cell that shows them all (the fit when none does).
   */
  zoomToBlocks(ranges: Int32Array, opts?: MoveOpts): void;
  /** Loci or `Coordinates` through the projector (no closure), as `focusExtents`; false: nothing shown here. */
  focusLoci(target: SeqLociAny | Coordinates, opts?: MoveOpts): boolean;
  /**
   * The hull of the target's projection (no closure) on each block of this section it lights: (block, start, stop)
   * triples in block indices, `stop` exclusive; a view of scratch, valid until the next call. Edge-rate.
   */
  hullOf(target: SeqLociAny | Coordinates): Int32Array;
  /** Whether the view already is what `focusLoci(target, { pad })` would move to; false: nothing shown here. */
  showsFocus(target: SeqLociAny | Coordinates, pad?: number): boolean;
  /** The hull of their projection on every block it lights (seqspec:interaction#focus); false: nothing shown here. */
  focusExtents(ext: Extents, opts?: MoveOpts): boolean;
  fit(opts?: MoveOpts): void;
  fitRows(): void;
  setRowScale(s: number, anchorY: number): void;
  /** Scroll groups. */
  setView(cell: number, scrollX: number, origin: Origin): void;
  /** The row scale and the vertical scroll at once (a window's top or bottom dragged). */
  setRowView(rowScale: number, scrollY: number, origin: Origin): void;
  /** Another viewport's zoom and horizontal scroll (scroll groups), read from its live record. */
  follow(from: SectionViewport, origin: Origin): void;
  /** The center's block (-1: none), its index with the fraction of its cell at `out[0]`, px a position at `out[1]`. */
  centerElement(out: Float64Array): number;
  /** Block `block`'s index `a[0]` (as `centerElement` gives) centered at `a[1]` px a shown position (`follow`). */
  centerOn(block: number, a: Float64Array, origin: Origin): void;
  /** From the engine or a test. A data width of 0 keeps the view; moves made meanwhile wait for a width. */
  setSize(width: number, height: number, dpr: number): void;
  /**
   * From the layout: `rows` holds the scrolled rows' height, the row scale that fits them, and the px above them (an
   * array: a double passed to a call V8 does not inline is boxed, each row zoom's event).
   */
  setRows(rows: Float64Array): void;
  /** A pointer gesture ended: commits now. */
  release(): void;
}
export type ViewportLive = { readonly current: ViewportState; readonly version: SeqState<number> };
export interface ViewportManager {
  get(section: string): SectionViewport | null;
  /** Every section of every view of the shown document. */
  readonly sections: SeqState<readonly string[]>;
}

export interface ViewportHost {
  readonly scheduler: Scheduler;
  readonly config: ConfigManager;
  active(): boolean;
}
/** What a viewport asks its manager: a section's hull of extents (`Viewports.hull`), loci as extents. */
interface ViewportOwner extends ViewportHost {
  hull(section: number, ext: Extents, out: BlockHulls): number;
  extents(target: SeqLociAny | Coordinates, out: Extents): number;
}

const EMPTY = new Int32Array(0);
const INITIAL: CommittedViewport = {
  origin: 'api',
  cell: 1,
  scrollX: 0,
  scrollY: 0,
  rowScale: 1,
  blocks: [],
  rows: [0, -1],
};
/** A pending intent: none, a block range, a cell centered at a view x. */
const Intent = { none: 0, range: 1, centered: 2 } as const;
/** Shown runs of a range (two across a cyclic origin), and a hull. */
const runs = new Int32Array(4);
const hullOut = BlockHulls.create();
/** One block's range, as the many-block verbs take it. */
const one = new Int32Array(3);

export class Viewport implements SectionViewport {
  readonly ref: string;
  readonly live: WritableLive<ViewportState>;
  /** Built from the last commit when read or delivered: a programmatic move nobody reads allocates nothing. */
  readonly committed = pulled<CommittedViewport>(INITIAL, () => this.publish());
  axis: Axis | null = null;
  /** The section's render props, set by the layout (a field: the clamp reads it per move). */
  props: ViewportProps = DefaultViewportProps;
  /** The section's index in the current view; -1: not in it. */
  index = -1;
  private readonly host: ViewportOwner;
  private fitted = true;
  /** A change not committed yet. */
  private dirty = false;
  private timer = -1;
  private lastInput = 0;
  private settleOrigin: Origin = 'api';
  private frame = -1;
  /** An animation's cell and scroll from and to; `y1` NaN: the vertical scroll is not animated. */
  private readonly anim = {
    c0: 1,
    x0: 0,
    c1: 1,
    x1: 0,
    y0: Number.NaN,
    y1: Number.NaN,
    t0: 0,
    d: 0,
    origin: 'api' as Origin,
  };
  /** `locate`'s: the view x to find (in), then the shown position (fractional) under it (out). */
  private readonly at = new Float64Array(2);
  /** The plain verbs' numbers, passed to their `…From` variants. */
  private readonly args = new Float64Array(2);
  /** `focusLoci`'s extents, made on first use. */
  private ext: Extents | null = null;
  private readonly prev = new Float64Array(4);
  /** Px above the scrolled rows (navigator, pinned entries), from the layout. */
  private insetY = 0;
  /** The visible content's [y0, y1), then its entries [first, last]: no double crosses a call. */
  private readonly range = new Float64Array(2);
  /** The pending intent (`Intent`) and its numbers: block, start, stop, pad; or cell, x. */
  private pending: number = Intent.none;
  private readonly want = new Float64Array(4);
  /** A pending range move's (block, start, stop) triples; `want[0]` of them. */
  private wantRanges = new Int32Array(3);
  private wantOrigin: Origin = 'api';
  /** The last commit, not yet in `committed`: cell, scrollX, scrollY, rowScale, first and last entry; runs; axis. */
  private readonly kept = new Float64Array(6);
  private keptLo = EMPTY;
  private keptHi = EMPTY;
  private keptAxis: Axis | null = null;
  private keptOrigin: Origin = 'api';
  private unpublished = false;
  /** A row zoom's anchored scrollY and its row scale, applied on the next `setRows`; NaN: none. */
  private readonly rowAnchor = new Float64Array([Number.NaN, 1]);

  constructor(ref: string, host: ViewportOwner) {
    this.ref = ref;
    this.host = host;
    this.live = live<ViewportState>({
      width: 0,
      height: 0,
      dpr: 1,
      cell: 1,
      minCell: 1,
      maxCell: 1,
      scrollX: 0,
      scrollY: 0,
      rowScale: 1,
      minRowScale: 1,
      maxRowScale: 1,
      contentWidth: 0,
      contentHeight: 0,
      lod: Lod.aggregate,
      visLo: EMPTY,
      visHi: EMPTY,
      firstEntry: 0,
      lastEntry: -1,
      animating: false,
      interacting: false,
      origin: 'api',
    });
  }

  scrollTo(x: number, y: number, origin: Origin = 'api') {
    this.axis?.ensure?.();
    const s = this.start();
    s.scrollX = x;
    s.scrollY = y;
    this.moved(origin);
  }
  /**
   * The vertical scroll to `y`, animated as a zoom to a focus is (`opts.animate`, `durationMs`; reduced motion and an
   * inactive context jump): wrapped text's scroll to a focus (`Layouts.reveal`). Edge-rate.
   */
  scrollYTo(y: number, opts: MoveOpts = {}) {
    const s = this.live.current;
    this.moveTo(s.cell, s.scrollX, opts, y);
  }
  panBy(dx: number, dy: number, origin: Origin = 'pointer') {
    this.args[0] = dx;
    this.args[1] = dy;
    this.panByFrom(this.args, origin);
  }
  /**
   * The verbs that carry a computed double have a `…From` variant reading its numbers from a Float64Array its caller
   * writes once and reuses: a double passed to a call V8 does not inline is boxed, and the wheel, a window's drag and
   * the keyboard make these calls at pointer rate. The plain verbs call them (edge rate).
   * `panByFrom`: dx at `a[0]`, dy at `a[1]`.
   */
  panByFrom(a: Float64Array, origin: Origin = 'pointer') {
    this.axis?.ensure?.();
    const s = this.start();
    s.scrollX += a[0];
    s.scrollY += a[1];
    this.moved(origin);
  }
  /** A drag's move (`dragPan`), read from its record here: a double passed to a call V8 does not inline is boxed. */
  panDrag(d: { readonly dx: number; readonly dy: number }, x: boolean, y: boolean) {
    this.axis?.ensure?.();
    const s = this.start();
    if (x) s.scrollX -= d.dx;
    if (y) s.scrollY -= d.dy;
    this.moved('pointer');
  }
  zoomAt(factor: number, x: number, origin: Origin = 'pointer') {
    this.args[0] = factor;
    this.args[1] = x;
    this.zoomAtFrom(this.args, origin);
  }
  /** `zoomAt`: the factor at `a[0]`, the view x it is anchored at `a[1]`. */
  zoomAtFrom(a: Float64Array, origin: Origin = 'pointer') {
    this.at[1] = a[1];
    const b = this.locate();
    const ax = this.axis;
    if (!ax || b < 0) return;
    const s = this.start();
    s.cell = Math.max(s.minCell, Math.min(s.maxCell, s.cell * a[0]));
    const x0 = ax.offset[b] * s.cell + b * ax.gap;
    s.scrollX = x0 + this.at[0] * ax.span[b] * s.cell - a[1];
    this.moved(origin);
  }
  zoomTo(block: number, start: number, stop: number, opts: MoveOpts = {}) {
    one[0] = block;
    one[1] = start;
    one[2] = stop;
    this.zoomToBlocks(one, opts);
  }
  zoomToBlocks(ranges: Int32Array, opts: MoveOpts = {}) {
    if (this.deferRanges(ranges, opts.pad ?? 0, opts.origin)) return;
    const cell = this.target(ranges, ranges.length / 3, opts.pad ?? 0);
    if (cell > 0) this.moveTo(cell, this.anim.x1, opts);
  }
  /** To `cell` (clamped), the position under view x centered (the navigator's double click). */
  zoomCentered(cell: number, x: number, opts: MoveOpts = {}) {
    if (this.defer(Intent.centered, cell, x, 0, 0, opts.origin)) return;
    this.at[1] = x;
    const b = this.locate();
    const ax = this.axis;
    if (!ax || b < 0) return;
    const s = this.live.current;
    const c = clamp(cell, s.minCell, s.maxCell);
    this.moveTo(c, this.xOf(b, this.at[0], c) - (s.width - ax.insetX) / 2, opts);
  }
  focusLoci(target: SeqLociAny | Coordinates, opts?: MoveOpts): boolean {
    this.ext ??= Extents.create(4);
    const ext = this.ext;
    ext.count = 0;
    this.host.extents(target, ext);
    return this.focusExtents(ext, opts);
  }
  hullOf(target: SeqLociAny | Coordinates): Int32Array {
    this.ext ??= Extents.create(4);
    const ext = this.ext;
    ext.count = 0;
    this.host.extents(target, ext);
    const n = this.index < 0 ? 0 : this.host.hull(this.index, ext, hullOut);
    return hullOut.data.subarray(0, 3 * n);
  }
  /**
   * Whether the view already is what `focusLoci(target, { pad })` would move to (the cell and the scroll within a
   * device px); false where the target is shown nowhere here. Edge-rate (a double click).
   */
  showsFocus(target: SeqLociAny | Coordinates, pad = 0): boolean {
    this.ext ??= Extents.create(4);
    const ext = this.ext;
    ext.count = 0;
    this.host.extents(target, ext);
    const n = this.index < 0 ? 0 : this.host.hull(this.index, ext, hullOut);
    if (n === 0) return false;
    const s = this.live.current;
    const keep = this.anim.x1;
    const cell = this.target(hullOut.data.subarray(0, 3 * n), n, pad);
    const x = this.anim.x1;
    this.anim.x1 = keep;
    return cell > 0 && Math.abs(cell - s.cell) <= cell * 1e-3 && Math.abs(x - s.scrollX) <= 1 / s.dpr + 0.5;
  }
  focusExtents(ext: Extents, opts?: MoveOpts): boolean {
    const n = this.index < 0 ? 0 : this.host.hull(this.index, ext, hullOut);
    if (n > 0) this.zoomToBlocks(hullOut.data.subarray(0, 3 * n), opts);
    return n > 0;
  }
  fit(opts: MoveOpts = {}) {
    this.pending = Intent.none;
    this.moveTo(this.live.current.minCell, 0, opts);
    this.fitted = true;
  }
  fitRows() {
    this.axis?.ensure?.();
    const s = this.start();
    s.scrollY = 0;
    s.rowScale = s.minRowScale;
    this.moved('api');
  }
  setRowScale(r: number, anchorY: number) {
    this.args[0] = r;
    this.args[1] = anchorY;
    this.setRowScaleFrom(this.args);
  }
  /** `setRowScale`: the scale at `a[0]`, the anchor's y at `a[1]`. */
  setRowScaleFrom(a: Float64Array) {
    this.axis?.ensure?.();
    const s = this.start();
    const next = clamp(a[0], s.minRowScale, s.maxRowScale);
    s.scrollY = ((s.scrollY + a[1]) * next) / s.rowScale - a[1];
    s.rowScale = next;
    // Clamped below against the rows' height at the old scale; the layout's next `setRows` restores the anchor.
    this.rowAnchor[0] = s.scrollY;
    this.rowAnchor[1] = next;
    this.moved('pointer');
  }
  setView(cell: number, scrollX: number, origin: Origin) {
    this.args[0] = cell;
    this.args[1] = scrollX;
    this.setViewFrom(this.args, origin);
  }
  /** `setView`: the cell at `a[0]`, scrollX at `a[1]`. */
  setViewFrom(a: Float64Array, origin: Origin) {
    const s = this.start();
    s.cell = a[0];
    s.scrollX = a[1];
    this.moved(origin);
  }
  setRowView(rowScale: number, scrollY: number, origin: Origin) {
    this.args[0] = rowScale;
    this.args[1] = scrollY;
    this.setRowViewFrom(this.args, origin);
  }
  /** `setRowView`: the row scale at `a[0]`, scrollY at `a[1]`. */
  setRowViewFrom(a: Float64Array, origin: Origin) {
    this.axis?.ensure?.();
    const s = this.start();
    s.rowScale = a[0];
    s.scrollY = a[1];
    this.moved(origin);
  }
  follow(from: SectionViewport, origin: Origin) {
    // Read here, not passed: a double passed to a call V8 does not inline is boxed.
    const f = from.live.current;
    const s = this.start();
    s.cell = f.cell;
    s.scrollX = f.scrollX;
    this.moved(origin);
  }
  /**
   * The element at the data area's center (`follow`, LK1): its block (-1: no axis, width or position there); `out[0]`
   * its block coordinate system index with the fraction of its cell, read in index order (`i + f`, `i` the element),
   * `out[1]` the px a shown position takes. Allocation-free; no double crosses a call.
   */
  centerElement(out: Float64Array): number {
    const ax = this.axis;
    const s = this.live.current;
    if (!ax || s.width - ax.insetX <= 0) return -1;
    this.at[1] = (s.width - ax.insetX) / 2;
    const b = this.locate();
    if (b < 0 || ax.shown[b] === 0) return -1;
    const k = Math.min(Math.floor(this.at[0]), ax.shown[b] - 1);
    const i = ax.indexOf(b, k as ShownIndex);
    if (i < 0) return -1;
    const f = this.at[0] - k;
    out[0] = ax.dir[b] > 0 ? i + f : i + 1 - f;
    out[1] = ax.span[b] * s.cell;
    return b;
  }
  /**
   * Centers block `b`'s index `a[0]` (as `centerElement` gives one: an element and its fraction, in index order) at
   * `a[1]` px per shown position (`follow`). An index not shown centers the nearest shown position.
   */
  centerOn(b: number, a: Float64Array, origin: Origin) {
    const ax = this.axis;
    const s = this.live.current;
    if (!ax || b < 0 || b >= ax.count || s.width - ax.insetX <= 0) return;
    const i = Math.floor(a[0]);
    let k: number = ax.shownOf(b, i);
    let g = a[0] - i;
    if (k < 0) {
      k = ax.nearestShown(b, i);
      g = 0.5;
    }
    if (k < 0) return;
    const pos = ax.dir[b] > 0 ? k + g : k + 1 - g;
    const cell = Math.max(s.minCell, Math.min(s.maxCell, a[1] / ax.span[b]));
    const x = ax.offset[b] * cell + b * ax.gap + pos * ax.span[b] * cell - (s.width - ax.insetX) / 2;
    // Where it already is (snapped and clamped as `clamp` would), nothing: an animation or a pending intent goes on.
    const max = s.contentWidth - (s.width - ax.insetX);
    if (cell === s.cell && Math.round(Math.max(0, Math.min(max, x)) * s.dpr) / s.dpr === s.scrollX) return;
    this.start();
    s.cell = cell;
    s.scrollX = x;
    this.moved(origin);
  }
  setSize(width: number, height: number, dpr: number) {
    const s = this.live.current;
    // The data columns that fit the width first: a narrower width can still leave a data area.
    this.axis?.fitWidth(width);
    const hidden = width - (this.axis?.insetX ?? 0) <= 0;
    const taller = width === s.width && height !== s.height;
    s.width = width;
    s.height = height;
    s.dpr = dpr;
    // Hidden: the view stays as it is, for the next size with room, and the layout waits at its width. A height alone
    // is the layout's to take (else a box that grew after the first layout keeps that layout's rows: one). What the
    // clamp moves is no one's move (`graph`, seqspec:viewer#state): a follower whose inset or DPR changed
    // at the same width leads nothing (VL-03).
    if (!hidden) this.refresh('graph');
    else if (taller) this.live.bump();
  }
  setRows(rows: Float64Array) {
    const contentHeight = rows[0];
    const minRowScale = rows[1];
    const insetY = rows[2];
    const s = this.live.current;
    // A row zoom's anchor (at this row scale) waits for the rows' new height, then is clamped against it.
    const anchored = !Number.isNaN(this.rowAnchor[0]) && this.rowAnchor[1] === s.rowScale;
    if (!anchored && s.contentHeight === contentHeight && s.minRowScale === minRowScale && this.insetY === insetY)
      return;
    this.begin();
    s.contentHeight = contentHeight;
    s.minRowScale = minRowScale;
    this.insetY = insetY;
    if (anchored) s.scrollY = this.rowAnchor[0];
    this.rowAnchor[0] = Number.NaN;
    // A clamp that moved the view commits (unless a gesture will); new extents alone only bump.
    if (this.end('graph')) {
      if (!s.interacting && !s.animating) this.commit('graph');
    } else this.live.bump();
  }
  /**
   * An export's run (seqspec:viewer#the-render-backend): `committed` emits nothing until the returned call. Then a view
   * put back as it
   * was commits nothing (a change pending before is pending again), and a view left elsewhere commits once.
   */
  quiet(): () => void {
    const s = this.live.current;
    const was = Float64Array.of(s.cell, s.scrollX, s.scrollY, s.rowScale);
    const pending = this.dirty;
    this.muted++;
    return () => {
      if (--this.muted > 0) return;
      const t = this.live.current;
      if (t.cell === was[0] && t.scrollX === was[1] && t.scrollY === was[2] && t.rowScale === was[3])
        this.dirty = pending;
      else this.commit('api');
    };
  }
  private muted = 0;
  release() {
    if (this.timer >= 0) this.host.scheduler.clearTimer(this.timer);
    this.timer = -1;
    this.settled();
  }

  /**
   * A new axis from the graph: an equal axis keeps the viewport; a changed one keeps the center's
   * block location and the zoom (`keepCenter`, when that location is still shown and the view did not show the whole
   * axis) or fits; a new section opens on
   * its focus, else fits. A focus before the section has a width waits for one.
   */
  resolve(axis: Axis, focus: { block: number; start: number; stop: number } | null, mode: 'keepCenter' | 'fit') {
    const old = this.axis;
    const s = this.live.current;
    if (old?.equals(axis)) {
      this.axis = axis;
      if (!old.sameFrame(axis)) this.refresh('graph');
      return;
    }
    let key = '';
    let index = -1;
    let frac = 0;
    // A view of the whole axis (fitted, or zoomed out to it) shows the whole new one: a filter's back to everything.
    const whole = this.fitted || s.cell <= s.minCell;
    if (old && mode === 'keepCenter' && !whole && s.width - old.insetX > 0) {
      this.at[1] = (s.width - old.insetX) / 2;
      const b = this.locate();
      if (b >= 0 && old.shown[b] > 0) {
        const k = Math.min(Math.floor(this.at[0]), old.shown[b] - 1);
        [key, index, frac] = [old.keys[b], old.indexOf(b, k as ShownIndex), this.at[0] - k];
      }
    }
    this.stop();
    this.axis = axis;
    if (focus) {
      one[0] = focus.block;
      one[1] = focus.start;
      one[2] = focus.stop;
    }
    if (!old && focus && this.deferRanges(one, 0, 'graph')) return this.refresh('graph');
    this.clamp();
    const nb = index < 0 ? -1 : axis.keys.indexOf(key);
    const k = nb < 0 ? -1 : axis.shownOf(nb, index);
    const cell = !old && focus ? this.target(one, 1, 0) : 0;
    this.fitted = k < 0 && cell <= 0 && this.pending === Intent.none;
    if (cell > 0) [s.cell, s.scrollX] = [cell, this.anim.x1];
    else if (k >= 0) s.scrollX = this.xOf(nb, k + frac, s.cell) - (s.width - axis.insetX) / 2;
    this.refresh('graph');
  }

  /** A move made while the data width is 0 waits for a width, the last one winning; true: deferred. */
  private defer(kind: number, a: number, b: number, c: number, pad: number, origin: Origin = 'api'): boolean {
    if (this.live.current.width - (this.axis?.insetX ?? 0) > 0) return false;
    this.stop();
    this.pending = kind;
    [this.want[0], this.want[1], this.want[2], this.want[3]] = [a, b, c, pad];
    this.wantOrigin = origin;
    this.fitted = false;
    return true;
  }

  /** A range move made while the data width is 0: its ranges kept (copied), as `defer`'s numbers are. */
  private deferRanges(ranges: Int32Array, pad: number, origin: Origin = 'api'): boolean {
    if (!this.defer(Intent.range, ranges.length / 3, 0, 0, pad, origin)) return false;
    if (this.wantRanges.length < ranges.length) this.wantRanges = new Int32Array(ranges.length);
    this.wantRanges.set(ranges);
    return true;
  }

  /**
   * The cell and (in `anim.x1`) the scroll that center the hull of `n` block ranges, (block, start, stop) triples:
   * from the leftmost shown start to the rightmost shown stop, gaps between blocks included, at the largest cell that
   * shows it all (the fit when even that is too wide); 0: nothing shown. x = cell·(offset + shown·span) + block·gap.
   */
  private target(ranges: Int32Array, n: number, pad: number): number {
    const ax = this.axis;
    if (!ax) return 0;
    // Plain locals, no destructuring: a pan's step may land here (a zoom per frame), and must allocate nothing.
    let lo = Number.POSITIVE_INFINITY;
    let hi = Number.NEGATIVE_INFINITY;
    let gLo = 0;
    let gHi = 0;
    for (let i = 0; i < n; i++) {
      const b = ranges[3 * i];
      if (b < 0 || b >= ax.count) continue;
      const s = ranges[3 * i + 1];
      const t = ranges[3 * i + 2];
      // A range across a cyclic origin (start > stop) is its two runs, [start, L) then [0, stop).
      const wraps = s > t;
      for (let part = 0; part < (wraps ? 2 : 1); part++) {
        const m = ax.runToShown(b, wraps && part ? 0 : s, wraps && !part ? ax.base[b] + ax.length[b] : t, runs);
        if (m === 0) continue;
        const a0 = ax.offset[b] + runs[0] * ax.span[b];
        const a1 = ax.offset[b] + runs[2 * m - 1] * ax.span[b];
        if (a0 < lo) {
          lo = a0;
          gLo = b * ax.gap;
        }
        if (a1 > hi) {
          hi = a1;
          gHi = b * ax.gap;
        }
      }
    }
    if (!(hi > lo)) return 0;
    const s = this.live.current;
    const w = s.width - ax.insetX;
    const cell = clamp((w / (1 + 2 * pad) - (gHi - gLo)) / (hi - lo), s.minCell, s.maxCell);
    this.anim.x1 = clamp((cell * (lo + hi) + gLo + gHi) / 2 - w / 2, 0, this.maxScroll(cell));
    return cell;
  }

  /**
   * To `cell` and scroll x, and (unless NaN) scroll y, over the animation `opts` ask for (seqspec:interaction#focus: a
   * zoom to a focus, wrapped text's scroll to one).
   */
  private moveTo(cell: number, x: number, opts: MoveOpts, y = Number.NaN) {
    const s = this.live.current;
    const origin = opts.origin ?? 'api';
    const config = this.host.config;
    const d = opts.durationMs ?? config.get(SeqViewerConfig.Viewport.AnimationMs);
    const reduced = config.get(SeqViewerConfig.Viewport.ReducedMotion) === true;
    this.start();
    if (!opts.animate || d <= 0 || reduced || !this.host.active()) {
      s.cell = cell;
      s.scrollX = x;
      if (!Number.isNaN(y)) s.scrollY = y;
      this.moved(origin);
      return;
    }
    const c1 = clamp(cell, s.minCell, s.maxCell);
    const x1 = clamp(x, 0, this.maxScroll(c1));
    Object.assign(this.anim, { c0: s.cell, x0: s.scrollX, c1, x1, y0: s.scrollY, y1: y, d, origin });
    this.anim.t0 = this.host.scheduler.now();
    // Interacting through the animation; its last frame clears it.
    s.animating = true;
    s.interacting = true;
    this.fitted = false;
    this.frame = this.host.scheduler.requestFrame(this.onFrame);
  }

  private readonly onFrame = () => {
    this.frame = -1;
    const a = this.anim;
    const s = this.live.current;
    const k = Math.min(1, (this.host.scheduler.now() - a.t0) / a.d);
    const e = k < 0.5 ? 2 * k * k : 1 - ((2 - 2 * k) * (2 - 2 * k)) / 2;
    this.begin();
    s.cell = k >= 1 ? a.c1 : Math.exp(Math.log(a.c0) + (Math.log(a.c1) - Math.log(a.c0)) * e);
    // A zoom about the one view x where the start and the end show the same position, so every position moves one
    // way: the scroll linear in `e` while the cell grows exponentially swept past the target and back once the cell
    // grew more than e-fold (a double click on a long sequence's full view). A pan (one cell) is linear.
    const dc = a.c1 - a.c0;
    s.scrollX =
      k >= 1
        ? a.x1
        : dc > 1e-3 * a.c0 || dc < -1e-3 * a.c0
          ? a.x0 + ((a.x1 - a.x0) * (s.cell - a.c0)) / dc
          : a.x0 + (a.x1 - a.x0) * e;
    if (!Number.isNaN(a.y1)) s.scrollY = k >= 1 ? a.y1 : a.y0 + (a.y1 - a.y0) * e;
    if (k < 1) {
      this.end(a.origin);
      this.frame = this.host.scheduler.requestFrame(this.onFrame);
      return;
    }
    s.animating = false;
    if (!this.end(a.origin, true)) {
      s.interacting = false;
      this.live.bump();
    }
    this.commit(a.origin);
  };

  /** Ends an animation where it is (a later move replaces it). */
  private stop() {
    if (this.frame < 0) return;
    this.host.scheduler.cancelFrame(this.frame);
    this.frame = -1;
    this.live.current.animating = false;
    this.live.current.interacting = false;
  }
  /** Jumps a running animation to its end (`setActive(false)`). */
  finish() {
    if (this.frame < 0) return;
    this.host.scheduler.cancelFrame(this.frame);
    this.anim.t0 = Number.NEGATIVE_INFINITY;
    this.onFrame();
  }

  /**
   * A change is `begin` (the view as it was, kept in a typed array), writes to the live record, then `end`. Doubles
   * never cross a call as arguments or results here: a call V8 does not inline would box them.
   */
  private begin(): Float64Array {
    const s = this.live.current;
    const p = this.prev;
    p[0] = s.cell;
    p[1] = s.scrollX;
    p[2] = s.scrollY;
    p[3] = s.rowScale;
    return p;
  }
  /** A user or host move begins: an animation ends where it is, and a pending intent or row anchor is dropped. */
  private start(): ViewportState {
    this.stop();
    this.pending = Intent.none;
    this.rowAnchor[0] = Number.NaN;
    this.begin();
    return this.live.current;
  }
  private moved(origin: Origin) {
    this.fitted = false;
    if (this.end(origin)) this.settle(origin);
  }
  /** Clamped and snapped; bumps only when it changed. An animation's last frame (`last`) clears `interacting`. */
  private end(origin: Origin, last = false): boolean {
    const s = this.live.current;
    const p = this.prev;
    this.clamp();
    if (p[0] === s.cell && p[1] === s.scrollX && p[2] === s.scrollY && p[3] === s.rowScale) return false;
    s.origin = origin;
    if (last) s.interacting = false;
    else if (!s.animating && !immediate(origin)) s.interacting = true;
    this.dirty = true;
    this.live.bump();
    return true;
  }

  /** `committed` emits at once for `api`, `graph` and `host:…`; else after `Viewport.SettleMs` of no change. */
  private settle(origin: Origin) {
    if (immediate(origin)) return this.commit(origin);
    this.lastInput = this.host.scheduler.now();
    this.settleOrigin = origin;
    if (this.timer < 0) this.timer = this.host.scheduler.setTimer(this.onSettle, this.settleMs());
  }
  private readonly onSettle = () => {
    this.timer = -1;
    const left = this.lastInput + this.settleMs() - this.host.scheduler.now();
    if (left > 0) this.timer = this.host.scheduler.setTimer(this.onSettle, left);
    else this.settled();
  };
  private settled() {
    const s = this.live.current;
    if (s.interacting) {
      s.interacting = false;
      this.live.bump();
    }
    this.commit(this.settleOrigin);
  }
  private settleMs() {
    return this.host.config.get(SeqViewerConfig.Viewport.SettleMs);
  }

  /**
   * After a size or axis change: re-clamped (a fitted view refits; a pending intent applies once there is a width),
   * bumped and committed.
   */
  private refresh(origin: Origin) {
    const s = this.live.current;
    const ax = this.axis;
    this.clamp();
    let o = origin;
    if (this.pending !== Intent.none && ax && s.width - ax.insetX > 0) {
      const w = this.want;
      const kind = this.pending;
      this.pending = Intent.none;
      o = this.wantOrigin;
      const cell = kind === Intent.range ? this.target(this.wantRanges, w[0], w[3]) : 0;
      if (cell > 0) [s.cell, s.scrollX] = [cell, this.anim.x1];
      this.at[1] = w[1];
      const b = kind === Intent.centered ? this.locate() : -1;
      if (b >= 0) {
        s.cell = clamp(w[0], s.minCell, s.maxCell);
        s.scrollX = this.xOf(b, this.at[0], s.cell) - (s.width - ax.insetX) / 2;
      }
      this.fitted = cell <= 0 && b < 0;
      this.clamp();
    }
    if (this.fitted) {
      s.cell = s.minCell;
      s.scrollX = 0;
      this.clamp();
    }
    s.origin = o;
    this.dirty = true;
    this.live.bump();
    if (!s.animating) this.commit(o);
  }

  /** The view as committed, kept; built into `committed` at once only when it has a subscriber (else when read). */
  private commit(origin: Origin) {
    if (!this.dirty || this.muted > 0) return;
    this.dirty = false;
    const s = this.live.current;
    const k = this.kept;
    k[0] = s.cell;
    k[1] = s.scrollX;
    k[2] = s.scrollY;
    k[3] = s.rowScale;
    k[4] = s.firstEntry;
    k[5] = s.lastEntry;
    if (this.keptLo.length !== s.visLo.length) {
      this.keptLo = new Int32Array(s.visLo.length);
      this.keptHi = new Int32Array(s.visLo.length);
    }
    this.keptLo.set(s.visLo);
    this.keptHi.set(s.visHi);
    this.keptAxis = this.axis;
    this.keptOrigin = origin;
    this.unpublished = true;
    if (this.committed.subscriberCount > 0) this.publish();
  }
  private publish() {
    if (!this.unpublished) return;
    this.unpublished = false;
    const k = this.kept;
    const ax = this.keptAxis;
    const blocks: { ref: string; start: number; stop: number }[] = [];
    const put = (b: number, lo: number, hi: number) => {
      if (hi <= lo || !ax) return;
      const a = ax.indexOf(b, lo as ShownIndex);
      const c = ax.indexOf(b, (hi - 1) as ShownIndex);
      blocks.push({ ref: ax.names[b], start: Math.min(a, c), stop: Math.max(a, c) + 1 });
    };
    for (let b = 0; ax && b < ax.count; b++) {
      const w = ax.wrapAt(b);
      const [lo, hi] = [this.keptLo[b], this.keptHi[b]];
      if (w > lo && w < hi) {
        put(b, lo, w);
        put(b, w, hi);
      } else put(b, lo, hi);
    }
    const [cell, scrollX, scrollY, rowScale] = [k[0], k[1], k[2], k[3]];
    const origin = this.keptOrigin;
    this.committed.set({ origin, cell, scrollX, scrollY, rowScale, blocks, rows: [k[4], k[5]] }, origin);
  }

  private gaps() {
    const ax = this.axis;
    return ax && ax.count > 1 ? (ax.count - 1) * ax.gap : 0;
  }
  private maxScroll(cell: number) {
    const ax = this.axis;
    return (ax?.cells ?? 0) * cell + this.gaps() - this.live.current.width + (ax?.insetX ?? 0);
  }
  /** Content x of a (fractional) shown position of a block at a cell width. */
  private xOf(b: number, shown: number, cell: number): number {
    const ax = this.axis as Axis;
    return ax.offset[b] * cell + b * ax.gap + shown * ax.span[b] * cell;
  }
  /**
   * The block under view x (read from `at[1]`: a double passed to a call is boxed; sets `at[0]`, the shown position
   * there, fractional); -1 without an axis. Block x is computed inline, here and in `clamp`, for the reason given there.
   */
  private locate(): number {
    const ax = this.axis;
    if (!ax || ax.count === 0) return -1;
    const s = this.live.current;
    const px = s.scrollX + this.at[1];
    let b = 0;
    let hi = ax.count - 1;
    if (ax.stacked) {
      // Stacked, every block starts at x 0: the first block reaching past px (on each, x maps to the same shown
      // index), else the widest (VL-01).
      let wide = 0;
      for (b = 0; b < ax.count; b++) {
        const w = ax.shown[b] * ax.span[b];
        if (w * s.cell > px) break;
        if (w > ax.shown[wide] * ax.span[wide]) wide = b;
      }
      if (b === ax.count) b = wide;
      hi = b;
    }
    // A binary search over offsets, inline (see above).
    while (b < hi) {
      const mid = (b + hi + 1) >> 1;
      if (ax.offset[mid] * s.cell + mid * ax.gap <= px) b = mid;
      else hi = mid - 1;
    }
    const x0 = ax.offset[b] * s.cell + b * ax.gap;
    this.at[0] = Math.max(0, Math.min(ax.shown[b], (px - x0) / (ax.span[b] * s.cell)));
    return b;
  }

  /** Limits, snapping, the LOD, the visible ranges and entries. Allocation-free unless the block count changed. */
  private clamp() {
    const s = this.live.current;
    const p = this.props;
    const ax = this.axis;
    const n = ax?.count ?? 0;
    const cells = ax?.cells ?? 0;
    const gaps = this.gaps();
    const dw = s.width - (ax ? ax.insetX : 0);
    const fixed = ax ? ax.fixedCell : 0;
    s.maxCell = fixed > 0 ? fixed : p.maxCell;
    // `minCell` (seqspec:model#pieces-and-layout-maps): when the gaps take more than three quarters, a quarter stays
    // for the positions.
    const room = Math.max(dw - gaps, dw / 4);
    s.minCell = fixed > 0 ? fixed : cells > 0 && dw > 0 ? Math.min(p.maxCell, room / cells) : p.maxCell;
    // Math's own calls only: a double returned from a call V8 does not inline is boxed (a zoom would allocate).
    s.cell = Math.max(s.minCell, Math.min(s.maxCell, s.cell));
    s.contentWidth = cells * s.cell + gaps;
    s.scrollX = fixed > 0 ? 0 : Math.round(Math.max(0, Math.min(s.contentWidth - dw, s.scrollX)) * s.dpr) / s.dpr;
    s.maxRowScale = p.maxRowScale;
    s.rowScale = Math.max(s.minRowScale, Math.min(s.maxRowScale, s.rowScale));
    const vh = s.height - this.insetY;
    s.scrollY = Math.round(Math.max(0, Math.min(s.contentHeight - vh, s.scrollY)) * s.dpr) / s.dpr;
    s.lod = s.cell >= p.minLetterCell ? Lod.letters : s.cell >= 3 ? Lod.cells : s.cell >= 1 ? Lod.thin : Lod.aggregate;
    if (s.visLo.length !== n) {
      s.visLo = new Int32Array(n);
      s.visHi = new Int32Array(n);
    }
    for (let b = 0; ax && b < n; b++) {
      const w = ax.span[b] * s.cell;
      const x0 = ax.offset[b] * s.cell + b * ax.gap;
      s.visLo[b] = fixed > 0 ? 0 : Math.max(0, Math.min(ax.shown[b], Math.floor((s.scrollX - x0) / w)));
      s.visHi[b] = fixed > 0 ? ax.shown[b] : Math.max(0, Math.min(ax.shown[b], Math.ceil((s.scrollX + dw - x0) / w)));
    }
    const pinned = ax ? ax.pinnedCount : 0;
    if (ax?.entriesIn) {
      this.range[0] = s.scrollY;
      this.range[1] = s.scrollY + vh;
      ax.entriesIn(this.range);
      s.firstEntry = this.range[0];
      s.lastEntry = this.range[1];
    } else {
      s.firstEntry = pinned;
      s.lastEntry = pinned - 1;
    }
  }

  /** The section left every view, or the context is disposed: timers cleared, states completed. */
  close() {
    this.stop();
    if (this.timer >= 0) this.host.scheduler.clearTimer(this.timer);
    this.live.complete();
    this.committed.complete();
  }
}

const immediate = (o: Origin) => o === 'api' || o === 'graph' || o.startsWith('host:');

export class Viewports implements ViewportManager {
  readonly sections = state<readonly string[]>([], sameList);
  private readonly map = new Map<string, Viewport>();
  /** Sections that left the shown document, closed at the end of the update unless they return in it. */
  private readonly leaving = new Map<string, Viewport>();
  private refs: string[] = [];
  private all: string[] = [];
  private readonly host: ViewportOwner;
  /** The edge-rate projector (seqspec:interaction#channels): focus and scroll-to, over the current view. */
  readonly projector = new Projector();

  constructor(host: ViewportHost) {
    this.host = {
      ...host,
      hull: (s, e, out) => this.projector.extents(e, s, out),
      extents: (t, out) => this.projector.fromLoci(t, out),
    };
  }

  get(section: string): Viewport | null {
    return this.map.get(section) ?? null;
  }
  sectionIndex(ref: string) {
    return this.refs.indexOf(ref);
  }
  sectionRef(i: number) {
    return this.refs[i] ?? null;
  }
  /** The current view's section refs, by number. */
  get current(): readonly string[] {
    return this.refs;
  }

  /**
   * Lifecycle: per-key states for every section of every view of the shown document. A ref that leaves and returns
   * within one update (`setDocument` after `setData`) keeps its objects; `publish` emits `sections`.
   */
  sync(doc: SeqViewSpecDoc | null) {
    const all = sectionsOf(doc);
    for (const ref of all) {
      if (this.map.has(ref)) continue;
      const back = this.leaving.get(ref);
      this.leaving.delete(ref);
      this.map.set(ref, back ?? new Viewport(ref, this.host));
    }
    for (const [ref, vp] of this.map) {
      if (all.includes(ref)) continue;
      this.map.delete(ref);
      if (this.leaving.size === 0) void Promise.resolve().then(this.drop);
      this.leaving.set(ref, vp);
    }
    this.all = all;
  }
  /** `sections`, after the layout and ui have synced. */
  publish() {
    this.sections.set(this.all, 'graph');
  }
  private readonly drop = () => {
    for (const vp of this.leaving.values()) vp.close();
    this.leaving.clear();
  };

  /** The current view's axes (from the layout, else built here with each section's props), on every graph emission. */
  update(view: ResolvedView | null, set: SeqSet | null, axisOf?: (section: number) => Axis | null) {
    this.refs = [];
    this.projector.update(view, set);
    for (const vp of this.map.values()) vp.index = -1;
    if (!view || !set) return;
    const s = view.sections;
    const mode = this.host.config.get(SeqViewerConfig.Viewport.OnAxisChange);
    for (let i = 0; i < s.count; i++) {
      const ref = set.strings.get(s.ref[i]);
      this.refs.push(ref);
      const vp = this.map.get(ref);
      if (!vp) continue;
      vp.index = i;
      const f = s.focus[i];
      const focus = f && { block: f.block - s.blockStart[i], start: f.start, stop: f.stop };
      const p = vp.props;
      const axis = axisOf ? axisOf(i) : new Axis(view, set, i, p.blockGap, p.ligandCells);
      if (axis) vp.resolve(axis, focus, mode);
    }
  }

  /**
   * The hull of the extents' projection (without closure, links and rows included) on every block of section `s` it
   * lights, in axis order, into `out` (block coordinate systems); the count of blocks lit.
   */
  hull(section: number, ext: Extents, out: BlockHulls): number {
    return this.projector.extents(ext, section, out);
  }

  finish() {
    for (const vp of this.map.values()) vp.finish();
  }
  dispose() {
    for (const vp of this.map.values()) vp.close();
    this.drop();
    this.map.clear();
    this.sections.complete();
  }
}

export interface ScrollGroupManager {
  /** Members by group: the current view's sections with one `scroll_group`, and joins from other contexts. */
  readonly groups: SeqState<ReadonlyMap<string, readonly SectionViewport[]>>;
  groupOf(section: string): string | null;
  /** Adds another context's viewport under the same rule (equal axes); returns the leave function. */
  join(name: string, viewport: SectionViewport): () => void;
}

interface Join {
  readonly name: string;
  readonly viewport: SectionViewport;
  leave(): void;
}

export class ScrollGroups implements ScrollGroupManager {
  readonly groups = state<ReadonlyMap<string, readonly SectionViewport[]>>(new Map());
  private local = new Map<string, string>();
  private readonly joins = new Set<Join>();
  private readonly viewports: Viewports;
  private readonly report: (code: string, message: string) => void;

  constructor(viewports: Viewports, report: (code: string, message: string) => void) {
    this.viewports = viewports;
    this.report = report;
  }

  update(view: ResolvedView | null, set: SeqSet | null) {
    this.local = new Map();
    for (let i = 0; view && set && i < view.sections.count; i++) {
      const g = view.sections.scrollGroup[i];
      if (g >= 0) this.local.set(set.strings.get(view.sections.ref[i]), set.strings.get(g));
    }
    this.publish();
  }
  groupOf(section: string) {
    return this.local.get(section) ?? null;
  }

  join(name: string, viewport: SectionViewport): () => void {
    const first = this.members(name)[0];
    const a = first instanceof Viewport ? first.axis : null;
    const b = viewport instanceof Viewport ? viewport.axis : null;
    if (a && b && !a.equals(b)) {
      this.report('svs.section.scroll_group', `"${viewport.ref}" cannot join "${name}": its axis differs`);
      return () => {};
    }
    let sub: { unsubscribe(): void } | null = null;
    const j: Join = {
      name,
      viewport,
      leave: () => {
        if (!this.joins.delete(j)) return;
        sub?.unsubscribe();
        this.publish();
      },
    };
    this.joins.add(j);
    // The other context's section leaving, or its dispose, completes the viewport: the join is left there too.
    sub = viewport.live.version.subscribe({ complete: j.leave });
    this.publish();
    return j.leave;
  }

  private members(name: string): SectionViewport[] {
    const out: SectionViewport[] = [];
    for (const [ref, g] of this.local) {
      const vp = g === name && this.viewports.get(ref);
      if (vp) out.push(vp);
    }
    for (const j of this.joins) if (j.name === name) out.push(j.viewport);
    return out;
  }
  private publish() {
    const names = new Set([...this.local.values(), ...[...this.joins].map((j) => j.name)]);
    this.groups.set(new Map([...names].map((n) => [n, this.members(n)])));
  }

  dispose() {
    for (const j of [...this.joins]) j.leave();
    this.groups.complete();
  }
}
