/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keyboard behavior and the cursor it keeps (seqspec:viewer#behaviors, seqspec:viewer#state): one tab stop per
// section, a cursor over the positions, rows and headers inside it (`ui.keyCursor`, and the hover with owner
// `keyboard`), zoom, row zoom, select as a click, activate, the track menu's key. The cursor is kept in the block
// coordinate system and repaired against the axis and layout it is read from, so a filter, a view switch or a relayout
// never leaves a stale index.

import { type Closure, Lanes, type ResolvedView } from "@molstar/seq-core/format";
import type { SeqSet, Tree } from "@molstar/seq-core/model";
import { type Axis, AxisGeometry, cellGap, type SectionAxisView, type ShownIndex } from "./axis.js";
import type { SeqViewerContext } from "./context.js";
import { defineBehavior, internal, viewportOf } from "./default-behaviors.js";
import { Hit, HitRegion, HoverKind } from "./hit.js";
import { type ClickEvt, Inputs, KeyBinding, type KeyEvt, Modifier, PointerType } from "./input.js";
import { EntryKind, type SectionLayout, SectionLayouts } from "./layout.js";
import { lanePad } from "./renderer-helpers.js";
import { Ring } from "./ring.js";
import { ringWindowZoom } from "./ring-window.js";
import { runXInto } from "./runs.js";
import type { KeyCursor } from "./ui.js";
import type { Viewport } from "./viewport.js";

const BINDINGS = {
  left: [{ key: "ArrowLeft" }],
  right: [{ key: "ArrowRight" }],
  up: [{ key: "ArrowUp" }],
  down: [{ key: "ArrowDown" }],
  home: [{ key: "Home" }],
  end: [{ key: "End" }],
  zoomIn: [{ key: "+" }, { key: "=" }],
  zoomOut: [{ key: "-" }, { key: "_" }],
  rowZoomIn: [{ code: "Equal", alt: true }],
  rowZoomOut: [{ code: "Minus", alt: true }],
  select: [{ key: "Enter" }],
  extend: [{ key: "Enter", shift: true }],
  toggle: [{ key: " " }],
  menu: [{ key: "F10", shift: true }, { key: "ContextMenu" }],
  clear: [{ key: "Escape" }],
} satisfies Record<string, KeyBinding>;
type Action = keyof typeof BINDINGS;
const ACTIONS = Object.keys(BINDINGS) as Action[];
/** Zoom steps, in and out; a key press is a step, not a wheel's continuous delta. */
const ZOOM_IN = 1.25;
const ZOOM_OUT = 0.8;
const NEAR = new Int32Array(4);
const XS = new Float64Array(1);
const NONE: KeyBinding = [];

/** The shown index nearest to block index `i` of block `b`; -1: nothing shown. Allocation-free. */
function nearest(ax: Axis, b: number, i: number): number {
  let best = -1;
  let gap = Number.POSITIVE_INFINITY;
  for (let side = 0; side < 2; side++) {
    // The runs of positions at or after i, then before it: each run's ends are its nearest and farthest.
    const n = side === 0 ? ax.runToShown(b, i, 1 << 30, NEAR) : ax.runToShown(b, 0, i, NEAR);
    for (let j = 0; j < 2 * n; j++) {
      const k = j & 1 ? NEAR[j] - 1 : NEAR[j];
      const d = Math.abs(ax.indexOf(b, k as ShownIndex) - i);
      if (d < gap) {
        best = k;
        gap = d;
      }
    }
  }
  return best;
}

/** The line of a wrapped section holding shown position k of block b: the last starting at or before it. */
function lineOf(l: SectionLayout, b: number, k: number): number {
  const t = l.lines;
  let lo = 0;
  for (let hi = t.count - 1; lo < hi; ) {
    const mid = (lo + hi + 1) >> 1;
    if (t.block[mid] < b || (t.block[mid] === b && t.start[mid] <= k)) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/** The clade a leaf's row folds into (its parent node); -1: none. Edge-rate (Space on a header). */
function cladeOf(tree: Tree | null, row: number): number {
  for (let n = 0; tree && row >= 0 && n < tree.count; n++) if (tree.leafRow[n] === row) return tree.parent[n];
  return -1;
}

/**
 * What a cursor stands on in the section's current layout, as numbers: the track `t`, the entry `e`, the layout's row
 * `d`, the block `b` and the shown index `k` (-1 on a header). One instance per user, reused; a stale cursor is
 * repaired to the nearest block, position, entry and row still there (`resolve(…, true)`), or refused.
 */
class Spot {
  ref = "";
  sec = -1;
  t = -1;
  e = -1;
  d = -1;
  b = 0;
  k = -1;
  header = false;
  /** The current `resolve` may move the cursor to the nearest place. */
  repair = false;
  l = null as unknown as SectionLayout;
  v = null as unknown as SectionAxisView;
  ax = null as unknown as Axis;
  vp = null as unknown as Viewport;
  rv = null as unknown as ResolvedView;
  set = null as unknown as SeqSet;
  /** The section's `rowBorders`: an entry's content stands above its border (`SectionLayouts.entryBorder`). */
  rowBorders = false;

  /** Binds the section's layout, axis view and viewport; false: nothing to stand on. */
  bind(ctx: SeqViewerContext, ref: string): boolean {
    const vp = viewportOf(ctx, ctx.layout.sectionIndex(ref));
    const a = ctx.graph.assembled.value?.view;
    const l = ctx.layout.get(ref)?.current;
    const v = ctx.layout.axis(ref);
    if (!vp?.axis || !a || !l || !v || vp.axis.count === 0) return false;
    this.ref = ref;
    this.sec = ctx.layout.sectionIndex(ref);
    this.vp = vp;
    this.ax = vp.axis;
    this.rv = a.resolved;
    this.set = a.set;
    this.l = l;
    this.v = v;
    this.rowBorders = ctx.layout.props(ref).rowBorders;
    return true;
  }
  drop() {
    this.l = this.v = this.ax = this.vp = this.rv = this.set = null as never;
  }

  /** The rows layer's display rows are this entry's (a lazy run is one entry for all of them). */
  inRows(e: number): boolean {
    const r = this.l.rows;
    return (
      r.layer >= 0 &&
      this.l.geometry === AxisGeometry.linear &&
      !this.ax.stacked &&
      (r.lazy ? e === r.entry : e >= r.entry && e < r.entry + r.count)
    );
  }
  /** A display row the cursor can stand on: a collapsed clade's wedge holds no element, except for a header. */
  live(d: number): boolean {
    const r = this.l.rows;
    return d >= 0 && d < r.count && (this.header || SectionLayouts.laneOf(r, d) >= 0);
  }
  private kind(e: number): boolean {
    const k = this.l.entryKind[e];
    return k === EntryKind.track || (this.header && k === EntryKind.group);
  }

  /** The first entry the cursor can stand on: the scrolled rows', then the pinned ones'. */
  first(): number {
    const l = this.l;
    for (let n = 0; n < l.entryCount; n++) {
      const i = (l.pinnedCount + n) % l.entryCount;
      if (this.kind(i)) return i;
    }
    return -1;
  }
  /**
   * The display row of data row `row` in the rows layer (a row folded into a clade: the clade's wedge); `d` itself while
   * it still shows that row, or for a lazy run. A collapse or a filter reorders display rows; the data row is the cursor's.
   */
  private displayOf(row: number, d: number): number {
    const r = this.l.rows;
    const lanes = this.rv.lanes;
    if (row < 0 || r.lazy || (d >= 0 && d < r.count && Lanes.rowOf(lanes, SectionLayouts.laneOf(r, d)) === row))
      return d;
    for (let k = 0; k < r.laneCount; k++) {
      if (Lanes.rowOf(lanes, r.laneStart + k) !== row) continue;
      for (let j = k; j >= 0; j--) {
        const x = SectionLayouts.displayOfLane(r, r.laneStart + j);
        if (x >= 0) return x;
      }
    }
    return d;
  }
  /**
   * Track t's entry at data row `row`, display row (or line) `d`; repairing, the nearest still there (by track, then
   * row): a collapsed split track's rows fall to its first. -1: none.
   */
  private entryFor(t: number, row: number, d: number): number {
    const l = this.l;
    const r = l.rows;
    if (this.inRows(r.entry) && t === r.track) {
      const want = this.displayOf(row, d);
      const x = Math.max(0, Math.min(r.count - 1, want));
      return x === want || this.repair ? SectionLayouts.entryOfDisplay(r, x) : -1;
    }
    const want = l.geometry === AxisGeometry.wrapped && !this.header ? lineOf(l, this.b, Math.max(0, this.k)) : d;
    let best = -1;
    let gap = Number.POSITIVE_INFINITY;
    for (let i = 0; i < l.entryCount && gap > 0; i++) {
      if (l.entryKind[i] !== EntryKind.track || !this.here(i)) continue;
      const g = Math.abs(l.entryTrack[i] - t) * 2 ** 24 + Math.abs(l.entryRow[i] - want);
      if (g < gap) {
        best = i;
        gap = g;
      }
    }
    return gap === 0 || this.repair ? best : -1;
  }
  /** A stacked section's entry is its block's (pinned ones, -1, are every block's). */
  private here(e: number): boolean {
    return !this.ax.stacked || this.l.entryBlock[e] < 0 || this.l.entryBlock[e] === this.b;
  }

  /** Sets the entry and what follows from it: the track, the row, and a stacked section's block. False: no live row. */
  at(e: number, dir: number): boolean {
    const l = this.l;
    const r = l.rows;
    this.e = e;
    if (l.entryKind[e] === EntryKind.group) {
      this.t = -1;
      this.d = l.entryGroup[e];
      return true;
    }
    this.t = l.entryTrack[e];
    let d = l.entryRow[e];
    if (this.inRows(e)) {
      d = r.lazy ? (dir > 0 ? 0 : r.count - 1) : d;
      while (r.lazy && d >= 0 && d < r.count && !this.live(d)) d += dir;
    }
    this.d = d;
    if (this.inRows(e) && !this.live(d)) return false;
    if (this.ax.stacked && l.entryBlock[e] >= 0) this.setBlock(l.entryBlock[e]);
    return true;
  }
  /** A new block keeps the shown index, clamped. */
  private setBlock(b: number) {
    this.b = b;
    if (this.k >= 0) this.k = Math.min(this.k, this.ax.shown[b] - 1);
  }

  /**
   * The cursor as numbers. Exact (`repair` false: any stale part refuses), or moved to the nearest block, position,
   * entry and row still shown (never reading a stale index).
   */
  resolve(ctx: SeqViewerContext, c: KeyCursor, repair: boolean): boolean {
    if (!this.bind(ctx, c.section)) return false;
    const { ax, l } = this;
    this.repair = repair;
    this.header = c.header;
    this.b = c.block < 0 ? 0 : c.block;
    if (this.b >= ax.count) {
      if (!repair) return false;
      this.b = ax.count - 1;
    }
    this.k = -1;
    if (!c.header) {
      this.k = c.index >= 0 ? ax.shownOf(this.b, c.index) : -1;
      if (this.k < 0) {
        if (!repair) return false;
        for (let n = 0; n < ax.count && ax.shown[this.b] === 0; n++) this.b = (this.b + 1) % ax.count;
        this.k = nearest(ax, this.b, c.index);
        this.header = this.k < 0;
      }
    }
    let e = -1;
    if (c.track === null) {
      for (let i = 0; i < l.entryCount && this.header && e < 0; i++)
        if (l.entryKind[i] === EntryKind.group && l.entryGroup[i] === c.displayRow) e = i;
    } else {
      let t = -1;
      for (let i = 0; i < this.rv.tracks.count && t < 0; i++)
        if (this.set.strings.get(this.rv.tracks.ref[i]) === c.track) t = i;
      if (t >= 0) e = this.entryFor(t, c.row, c.displayRow);
    }
    // A removed entry: the first one the layout has.
    if (e < 0 && repair) e = this.first();
    if (e < 0) return false;
    // A row that holds no element (a clade collapsed over it): the nearest entry that does.
    if (!this.at(e, 1) && !(repair && (this.step(1) || this.step(-1)))) return false;
    // A lazy run's row the layout no longer has, or that holds no element, is the nearest live one.
    if (this.e === e && this.inRows(e) && l.rows.lazy) {
      const want = Math.max(0, Math.min(l.rows.count - 1, c.displayRow));
      if (want !== c.displayRow && !repair) return false;
      this.d = want;
      for (let n = 1; !this.live(this.d) && n <= l.rows.count; n++)
        this.d = Math.max(0, Math.min(l.rows.count - 1, want + (n & 1 ? (n + 1) >> 1 : -(n >> 1))));
      if (!this.live(this.d) || (this.d !== c.displayRow && !repair)) return false;
    }
    return true;
  }

  /** The numbers the published cursor holds. */
  get index(): number {
    return this.header ? -1 : this.ax.indexOf(this.b, this.k as ShownIndex);
  }
  get block(): number {
    return this.header && !this.ax.stacked ? -1 : this.b;
  }
  get row(): number {
    const r = this.l.rows;
    return this.t >= 0 && this.inRows(this.e) ? Lanes.rowOf(this.rv.lanes, SectionLayouts.laneOf(r, this.d)) : -1;
  }
  get track(): string | null {
    return this.t < 0 ? null : this.set.strings.get(this.rv.tracks.ref[this.t]);
  }
  same(c: KeyCursor): boolean {
    return (
      c.section === this.ref &&
      c.header === this.header &&
      c.displayRow === this.d &&
      c.row === this.row &&
      c.index === this.index &&
      c.block === this.block &&
      c.track === this.track
    );
  }
  cursor(): KeyCursor {
    return {
      section: this.ref,
      track: this.track,
      header: this.header,
      block: this.block,
      index: this.index,
      row: this.row,
      displayRow: this.d,
    };
  }

  /**
   * The cursor's cell (or header) in surface px, not clipped to the view: x, y, w, h, then the point a hit test reads.
   * On rings (seqspec:viewer#the-render-backend) the cell's box on the surface and its middle, or with `frame` its rect
   * in the unrolled frame the rings are drawn from (the overlay's cursor ring). False: no such place.
   */
  rect(out: Float64Array, frame = false): boolean {
    const { l, v, ax, e } = this;
    let h = l.entryHeight[e];
    let y = e < l.pinnedCount ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - v.scrollY;
    if (this.inRows(e) && l.rows.lazy) {
      h = l.rows.height;
      y += this.d * h;
    } else if (SectionLayouts.entryBorder(l, e, this.rowBorders)) h -= (Math.max(1, Math.round(v.dpr)) | 0) / v.dpr;
    let w = l.headerWidth;
    let x = 0;
    out[4] = Math.min(12, w / 2);
    if (!this.header) {
      w = ax.span[this.b] * v.cell;
      if (l.geometry === AxisGeometry.wrapped) {
        const run = v.runs[this.b];
        if (!run) return false;
        runXInto(run, l.lines.start[this.d], this.k, v.gap, XS, 0);
        x = v.left + XS[0];
        w = run.advance;
        out[4] = x + w / 2;
      } else {
        x = v.x(this.b, this.k as ShownIndex);
        out[4] = v.blockX[this.b] + (this.k + 0.5) * w;
        // The cell as drawn, so the ring sits evenly around it: less its gap on the right (`cellGap`) and, off the
        // rows, the lane's padding above and below (`lanePad`).
        w -= cellGap(w, v.dpr);
        const pad = this.inRows(e) ? 0 : lanePad(h);
        y += pad;
        h -= 2 * pad;
      }
    }
    out[0] = x;
    out[1] = y;
    out[2] = w;
    out[3] = h;
    out[5] = y + h / 2;
    if (l.geometry === AxisGeometry.circular && !frame) onRing(l, v, out);
    return Number.isFinite(out[0] + out[1] + out[2] + out[3] + out[4]);
  }
  /** Whether the rect (from `rect`) meets the part of the surface its entry is drawn in. */
  inView(out: Float64Array): boolean {
    const l = this.l;
    // Rings show every position.
    if (l.geometry === AxisGeometry.circular) return true;
    const across = this.header || (out[0] + out[2] > l.dataLeft && out[0] < l.dataRight);
    return across && (this.e < l.pinnedCount || (out[1] + out[3] > l.scrolledTop && out[1] < l.height));
  }

  /** Pans the viewport so the rect is in view, centered when it was more than a view away. True: it moved. */
  reveal(out: Float64Array, by: Float64Array): boolean {
    const l = this.l;
    if (l.geometry === AxisGeometry.circular) return false;
    let dx = 0;
    let dy = 0;
    const area = l.dataRight - l.dataLeft;
    const rows = l.height - l.scrolledTop;
    if (!this.header && l.geometry === AxisGeometry.linear) {
      dx =
        out[0] < l.dataLeft ? out[0] - l.dataLeft : out[0] + out[2] > l.dataRight ? out[0] + out[2] - l.dataRight : 0;
      if (Math.abs(dx) > area) dx = out[0] + out[2] / 2 - (l.dataLeft + l.dataRight) / 2;
    }
    if (this.e >= l.pinnedCount) {
      dy =
        out[1] < l.scrolledTop ? out[1] - l.scrolledTop : out[1] + out[3] > l.height ? out[1] + out[3] - l.height : 0;
      if (Math.abs(dy) > rows) dy = out[1] + out[3] / 2 - (l.scrolledTop + l.height) / 2;
    }
    if (dx === 0 && dy === 0) return false;
    by[0] = dx;
    by[1] = dy;
    this.vp.panByFrom(by, "keyboard");
    return true;
  }

  /** On data at shown position k of block b (a wrapped section's entry is its line's); never a group's or a wedge's. */
  moveTo(b: number, k: number) {
    if (this.t < 0 || (this.inRows(this.e) && SectionLayouts.laneOf(this.l.rows, this.d) < 0)) return;
    this.b = b;
    this.k = k;
    this.header = false;
    const e = this.l.geometry === AxisGeometry.wrapped ? this.entryFor(this.t, -1, this.d) : -1;
    if (e >= 0) this.at(e, 1);
  }
  /** The first shown position in view (the first block with one; a stacked section's own block). */
  firstInView() {
    const { ax, vp } = this;
    const s = vp.live.current;
    for (let b = ax.stacked ? this.b : 0; b < (ax.stacked ? this.b + 1 : ax.count); b++)
      if (s.visHi[b] > s.visLo[b]) return this.moveTo(b, s.visLo[b]);
    for (let b = 0; b < ax.count; b++) if (ax.shown[b] > 0) return this.moveTo(b, 0);
  }
  /** The entry's last shown position. */
  last() {
    const { ax } = this;
    let b = ax.stacked ? this.b : ax.count - 1;
    while (b > 0 && !ax.stacked && ax.shown[b] === 0) b--;
    if (ax.shown[b] > 0) this.moveTo(b, ax.shown[b] - 1);
  }
  /**
   * One position along, across block gaps; off the first it is the header, off the last it stays. On rings (no header)
   * it goes round: off the last onto the first, and back.
   */
  along(dir: number) {
    const { ax } = this;
    const k = this.k + dir;
    const round = this.l.geometry === AxisGeometry.circular;
    if (k >= 0 && k < ax.shown[this.b]) return this.moveTo(this.b, k);
    let n = this.b + dir;
    for (let i = 0; i <= ax.count && !ax.stacked; i++, n += dir) {
      if (round) n = (n + ax.count) % ax.count;
      if (n < 0 || n >= ax.count || ax.shown[n] > 0) break;
    }
    if (!ax.stacked && n >= 0 && n < ax.count && ax.shown[n] > 0) this.moveTo(n, dir > 0 ? 0 : ax.shown[n] - 1);
    else if (k < 0 && !round) this.header = true;
  }
  /** The entry above (-1) or below (1) the cursor's. False: at the edge, and nothing moved. */
  step(dir: number): boolean {
    const l = this.l;
    const r = l.rows;
    if (this.inRows(this.e) && r.lazy)
      for (let d = this.d + dir; d >= 0 && d < r.count; d += dir)
        if (this.live(d)) {
          this.d = d;
          return true;
        }
    const e0 = this.e;
    const t0 = this.t;
    const d0 = this.d;
    const b0 = this.b;
    const k0 = this.k;
    for (let e = e0 + dir; e >= 0 && e < l.entryCount; e += dir) {
      if (l.geometry === AxisGeometry.wrapped && l.entryRow[e] !== l.entryRow[e0]) break;
      if (this.kind(e) && this.at(e, dir)) return true;
    }
    this.e = e0;
    this.t = t0;
    this.d = d0;
    this.b = b0;
    this.k = k0;
    return false;
  }
}

/** The one the draw side reads. */
const probe = new Spot();
const ring = new Ring();
const at = new Float64Array(2);

/** A rect of the unrolled frame (`Spot.rect`) onto the rings' surface: its corners' box, and its middle. */
function onRing(l: SectionLayout, v: SectionAxisView, out: Float64Array) {
  ring.from(l, v, 1);
  const x = out[0];
  const y = out[1];
  const w = out[2];
  const h = out[3];
  let x0 = Number.POSITIVE_INFINITY;
  let y0 = Number.POSITIVE_INFINITY;
  let x1 = Number.NEGATIVE_INFINITY;
  let y1 = Number.NEGATIVE_INFINITY;
  for (let c = 0; c < 4; c++) {
    ring.toSurface(x + (c & 1 ? w : 0), y + (c & 2 ? h : 0), at);
    x0 = Math.min(x0, at[0]);
    x1 = Math.max(x1, at[0]);
    y0 = Math.min(y0, at[1]);
    y1 = Math.max(y1, at[1]);
  }
  ring.toSurface(x + w / 2, y + h / 2, at);
  out[0] = x0;
  out[1] = y0;
  out[2] = x1 - x0;
  out[3] = y1 - y0;
  out[4] = at[0];
  out[5] = at[1];
}

/**
 * The cursor's place in the section's surface px (x, y, w, h into `out`, which may hold six numbers); false when it is
 * not in view or the layout has no such place. The overlay's ring and the keyboard's tooltip anchor read it; on rings
 * the overlay's with `frame` (the unrolled frame it draws in).
 */
export function cursorRect(ctx: SeqViewerContext, cursor: KeyCursor, out: Float64Array, frame = false): boolean {
  try {
    return probe.resolve(ctx, cursor, false) && probe.rect(out, frame) && probe.inView(out);
  } finally {
    probe.drop();
  }
}

export const KeyboardBehavior = defineBehavior(
  "keyboard",
  "interaction",
  0,
  { bindings: BINDINGS as Readonly<Record<Action, KeyBinding>>, closure: "containing" as Closure },
  "A cursor over positions, rows and headers (the hover, owner `keyboard`): arrows, Home and End; +/- zoom, Alt+=/Alt+- " +
    "row zoom; Enter and Space select or activate; Shift+F10 or ContextMenu asks for the track's menu (`events.menu`); " +
    "Escape clears the selection and puts the cursor away; a press puts it away too. Tab leaves the surface.",
  ({ ctx, h, on, add }) => {
    /** The key's cursor; `rs` is the settle's own, since it runs inside the other's `land`. */
    const r = new Spot();
    const rs = new Spot();
    const hit = Hit.create();
    /** The cell as a pointer record: a keyboard Enter is a click on it. */
    const click: ClickEvt = { input: Inputs.pointer(), hit, count: 1 };
    const rect = new Float64Array(6);
    const by = new Float64Array(2);
    const memo = new Map<string, KeyCursor>();
    let busy = false;
    /** A press came before the focus it gives (no cursor on it). */
    let pressed = false;
    let watched = "";
    let layoutSub: { unsubscribe(): void } | null = null;

    const drop = () => {
      if (ctx.ui.keyCursor.value) ctx.ui.setKeyCursor(null);
      ctx.hover.clear("keyboard");
    };
    /** The hover follows the cursor's cell; none where it is not in view. */
    const show = (s: Spot) => {
      const p = h.forSection(s.sec);
      const ok = s.rect(rect);
      hit.x = rect[4];
      hit.y = rect[5];
      if (p && ok && s.inView(rect) && ctx.layout.hitTestAt(s.sec, hit) && hit.kind !== HoverKind.none)
        ctx.hover.setFromHit(hit, "keyboard", p.closure);
      else ctx.hover.clear("keyboard");
    };
    /**
     * Publishes the spot (when it differs from the cursor) and follows it with the hover; `settled` (a relayout, not a
     * key) leaves a hover another owner holds alone: a later owner replaces an earlier (seqspec:interaction#hover). The
     * cursor's hover comes back at the next relayout or key once nobody holds one, not when the pointer leaves.
     */
    const publish = (s: Spot, settled = false) => {
      const c = ctx.ui.keyCursor.value;
      if (!c || !s.same(c)) {
        const next = s.cursor();
        memo.set(s.ref, next);
        ctx.ui.setKeyCursor(next);
      }
      const owner = ctx.hover.live.current.owner;
      if (!settled || owner === "keyboard" || owner === "") show(s);
    };
    /** The spot becomes the cursor, brought into view first. */
    const land = () => {
      ctx.layout.axis(r.ref);
      if (r.rect(rect) && r.reveal(rect, by)) ctx.layout.axis(r.ref);
      publish(r);
    };
    /** A change of axis or layout: the cursor moves to the nearest place still there, or goes. */
    const settle = () => {
      const c = ctx.ui.keyCursor.value;
      if (!c || busy) return;
      busy = true;
      try {
        if (rs.resolve(ctx, c, true)) publish(rs, true);
        else drop();
      } finally {
        busy = false;
        rs.drop();
      }
    };
    /** From the section's remembered cursor, else the first shown position in view of its first scrolled entry. */
    const start = (ref: string) => {
      const c = memo.get(ref);
      if (c && r.resolve(ctx, c, true)) return true;
      if (!r.bind(ctx, ref)) return false;
      r.header = false;
      r.repair = true;
      r.b = 0;
      r.k = -1;
      const e = r.first();
      // The first row may be a wedge: the first below it.
      if (e < 0 || (!r.at(e, 1) && !r.step(1))) return false;
      r.firstInView();
      r.header = r.k < 0;
      return true;
    };

    on(ctx.interaction.focus, (f) => {
      const ref = ctx.layout.sectionRef(f.section);
      const c = ctx.ui.keyCursor.value;
      if (ref === null || !h.forSection(f.section)) return;
      if (!f.focused) {
        // `publish` keeps the section's last cursor for the next focus.
        pressed = false;
        if (c?.section === ref) drop();
      } else if (pressed) pressed = false;
      else if (start(ref)) land();
      r.drop();
    });
    // A press puts the cursor away (the pointer is in use), and the focus it gives shows none until a key.
    on(ctx.interaction.pointer, (e) => {
      if (e.input.type !== PointerType.down) return;
      pressed = true;
      if (ctx.ui.keyCursor.value) drop();
    });
    // `r` lets go of the set and the layout after each key.
    on(ctx.interaction.key, (e) => {
      try {
        key(e);
      } finally {
        r.drop();
      }
    });
    const key = (e: KeyEvt) => {
      const k = e.input;
      const p = h.forSection(k.section);
      let act: Action | null = null;
      for (let i = 0; i < ACTIONS.length && act === null && p; i++)
        if (KeyBinding.match(p.bindings[ACTIONS[i]] ?? NONE, k)) act = ACTIONS[i];
      const ref = ctx.layout.sectionRef(k.section);
      if (e.handled || !p || act === null || ref === null) return;
      if (act === "clear") {
        // The cursor goes too; the next key brings it back where it was (`memo`).
        ctx.selection.clear("keyboard");
        drop();
        e.handled = true;
        return;
      }
      const c = ctx.ui.keyCursor.value;
      if (!(c?.section === ref ? r.resolve(ctx, c, true) : start(ref))) return;
      e.handled = true;
      if (act === "menu") {
        // The track menu is the UI's (`events.menu`): the cursor's track opens it; a group's header has none.
        if (r.t >= 0) internal(ctx).openMenu(r.ref, r.track as string, "keyboard");
        return;
      }
      const l = r.l;
      switch (act) {
        case "left":
        case "right":
          if (!r.header) r.along(act === "left" ? -1 : 1);
          else if (act === "left") break;
          // A wrapped header is its line's: onto the line's first position.
          else if (l.geometry === AxisGeometry.wrapped) r.moveTo(l.lines.block[r.d], l.lines.start[r.d]);
          else r.firstInView();
          break;
        case "up":
        case "down":
          r.step(act === "up" ? -1 : 1);
          break;
        case "home":
          // Rings have no header: their first position.
          if (l.geometry === AxisGeometry.circular) r.firstInView();
          else r.header = true;
          break;
        case "end":
          r.last();
          break;
        case "zoomIn":
        case "zoomOut":
          by[0] = act === "zoomIn" ? ZOOM_IN : ZOOM_OUT;
          // Rings show their whole axis: the keys zoom the section their window shows, about its center (CG3).
          if (l.geometry === AxisGeometry.circular) {
            ringWindowZoom(ctx, r.sec, by);
            break;
          }
          // Anchored at the cursor's cell, on a header at the view's center.
          r.rect(rect);
          by[1] = r.header ? (r.vp.live.current.width - r.ax.insetX) / 2 : rect[4] - r.v.left;
          r.vp.zoomAtFrom(by, "keyboard");
          break;
        case "rowZoomIn":
        case "rowZoomOut":
          r.rect(rect);
          by[0] = r.vp.live.current.rowScale * (act === "rowZoomIn" ? ZOOM_IN : ZOOM_OUT);
          by[1] = rect[5] - l.scrolledTop;
          r.vp.setRowScaleFrom(by);
          break;
      }
      land();
      if (act !== "select" && act !== "extend" && act !== "toggle") return;
      if (!r.header) {
        // On data, a click on the cell: the `select` behavior's policy, Shift for `extend`.
        const input = click.input;
        input.type = PointerType.up;
        input.section = r.sec;
        input.button = 0;
        input.x = hit.x;
        input.y = hit.y;
        input.modifiers = act === "extend" ? Modifier.shift : 0;
        if (hit.kind !== HoverKind.none && hit.region !== HitRegion.none && r.inView(rect))
          internal(ctx).input.synthClick(click);
      } else if (act !== "toggle") {
        if (r.t >= 0) internal(ctx).activate(r.ref, r.track as string, "keyboard");
      } else {
        // Space on a header expands or collapses a group, a split track, or a clade: a wedge's, a leaf's parent's.
        const str = (id: number) => r.set.strings.get(id);
        const aln = r.inRows(r.e) ? r.set.alignments[l.rows.alignment] : null;
        const node = aln ? (l.rows.node?.[r.d] ?? -1) : -1;
        const clade = aln && node < 0 ? cladeOf(aln.tree, r.row) : node;
        if (l.entryKind[r.e] === EntryKind.group) ctx.ui.toggleGroup(str(r.rv.groups.ref[l.entryGroup[r.e]]));
        else if (aln && clade >= 0) ctx.ui.toggleClade(str(aln.key), clade);
        else if (r.rv.tracks.split[r.t] >= 0) ctx.ui.toggleGroup(str(r.rv.tracks.ref[r.t]));
      }
    };
    add(ctx.graph.assembled.subscribe(settle));
    add(
      ctx.ui.keyCursor.subscribe((c) => {
        // The layout of the cursor's section; moves within it do not resubscribe.
        if ((c?.section ?? "") === watched) return;
        watched = c?.section ?? "";
        layoutSub?.unsubscribe();
        layoutSub = c ? (ctx.layout.get(c.section)?.version.subscribe(settle) ?? null) : null;
      }),
    );
    add(() => {
      layoutSub?.unsubscribe();
      r.drop();
      drop();
    });
  },
);
