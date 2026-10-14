/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The default gesture behaviors (seqspec:viewer#behaviors, seqspec:interaction#selection): the scrollbar, navigator and
// vertical overview (claimed on the press, so they skip the slop), headers, selection by click and drag, drag to pan,
// the wheel and row zoom (the keyboard is keyboard.ts). Who claims a drag is params (`dragRegions`, `dragOn`), never a
// rule in the dispatcher. Hot callers write their numbers into a reused Float64Array and call the viewport's `…From`
// verbs.

import { type Closure, Height, Lanes, SelectionGranularity } from "@molstar/seq-core/format";
import { Loci, OrderedSet } from "@molstar/seq-core/model";
import { AxisGeometry, navUnits, type SectionAxisView } from "./axis.js";
import { ChannelMarks } from "./channels.js";
import { SeqViewerConfig } from "./config.js";
import type { SeqViewerContext } from "./context.js";
import { defineBehavior, internal, RegionName, viewportOf } from "./default-behaviors.js";
import type { ViewModel } from "./graph.js";
import { CHEVRON_GAP, CHEVRON_HALF, CHEVRON_ROOM, CHEVRON_SLOP, splitChevron } from "./header-toggle.js";
import { type Hit, HitRegion, HoverKind, wrappedPositionNear } from "./hit.js";
import { Binding, type DragEvt, DragPhase, Modifier, PointerType, type Trigger } from "./input.js";
import { EntryKind, type SectionLayout, SectionLayouts } from "./layout.js";
import { carry, type SeqLociAny } from "./lineage.js";
import { MarkBit } from "./marks.js";
import { ofRuns } from "./projection.js";
import { RegionState } from "./streams.js";
import { WINDOW_FRAME, WindowPart, windowEnds, windowFrame, windowIdle, windowPart } from "./windows.js";

const W = WindowPart;

/**
 * The navigator's params: a double click toggles between the full view and a zoom centered there, showing about
 * `zoomPositions` positions across the data area (so the cell follows the view's width; at least the letters
 * threshold), or `zoomCell` px per position when set; `none` turns it off.
 */
const NAVIGATOR = { doubleClick: "toggle" as "toggle" | "none", zoomPositions: 60, zoomCell: 0 };
type NavigatorParams = typeof NAVIGATOR;
/** The double click's moves: animated as other viewport changes are. */
const TOGGLE_MOVE = { animate: true, origin: "pointer" } as const;

/**
 * The navigator (`kind` 0: x over the data area), the vertical overview (1: x and y over its column) and the
 * scrollbar (2: y). A press is claimed: on an edge it zooms, inside the window it pans, elsewhere it jumps there
 * (the scrollbar pages). Window fractions come from the layout's axis view (`navX0`…, `vnavY0`…).
 */
function windowBehavior<P extends object>(
  name: string,
  region: number,
  kind: 0 | 1 | 2,
  description: string,
  params: P,
) {
  return defineBehavior(name, "interaction", 100, params, description, ({ ctx, h, on }) => {
    let mode: number = W.none;
    /** The window's frame (windows.ts `windowFrame`), and at the press: its frame, and the cell, scrollX, rowScale,
     * scrollY and content height (an end's drag is from the press, by the pointer's travel since). */
    const g = new Float64Array(WINDOW_FRAME);
    /** The section's axis view the frame is from: the navigator's map goes through its blocks (`navUnits`). */
    let av: SectionAxisView | null = null;
    const frame = (section: number) => {
      const ref = ctx.layout.sectionRef(section);
      const v = ref === null ? null : ctx.layout.axis(ref);
      const l = ref === null ? null : ctx.layout.get(ref)?.current;
      if (!v || !l) return false;
      windowFrame(v, l, kind, g);
      av = v;
      return true;
    };
    const g0 = new Float64Array(WINDOW_FRAME);
    /** At the press: cell, scrollX, rowScale, scrollY, content height. */
    const at = new Float64Array(5);
    /** A point on the navigator's map in, its units out (`navUnits`). */
    const io = new Float64Array(2);
    /** `windowEnds`' travel and least in, the window's fractions out. */
    const ends = new Float64Array(2);
    /** The numbers a viewport verb reads. */
    const buf = new Float64Array(2);
    on(ctx.interaction.pointer, (e) => {
      const x = e.input.x;
      const y = e.input.y;
      if (e.input.type !== PointerType.down || e.hit.region !== region || e.owner !== "") return;
      const vp = viewportOf(ctx, e.hit.section);
      if (!vp || !h.forSection(e.hit.section) || !frame(e.hit.section)) return;
      // A navigator with nowhere to go (the whole axis in view at the zoom limit) is inert: no press, jump or drag.
      if (kind === 0 && windowIdle(vp.live.current, vp.axis)) return;
      e.claim(name);
      // The same classification the cursor shows (windows.ts).
      mode = windowPart(kind, e.input, g, ctx.layout.props(vp.ref).verticalZoom);
      const s = vp.live.current;
      g0.set(g);
      at[0] = s.cell;
      at[1] = s.scrollX;
      at[2] = s.rowScale;
      at[3] = s.scrollY;
      at[4] = s.contentHeight;
      if (mode !== W.page) return;
      const wy0 = g[2] + g[6] * g[3];
      const vh = g[3];
      if (kind === 2) vp.scrollTo(s.scrollX, s.scrollY + (y < wy0 ? -vh : vh), "pointer");
      else {
        io[0] = x;
        if (kind === 0) navUnits(av as SectionAxisView, io);
        const gap = vp.axis ? vp.axis.gap : 0;
        const cx = kind === 0 ? io[0] * s.cell + io[1] * gap : ((x - g[0]) / g[1]) * s.contentWidth;
        // Centered there: half the data area (the navigator's map is not linear in content).
        const half =
          kind === 0 ? Math.min(s.width - (vp.axis?.insetX ?? 0), s.contentWidth) : (g[5] - g[4]) * s.contentWidth;
        const sx = cx - half * 0.5;
        const sy = kind === 1 ? ((y - g[2]) / g[3]) * s.contentHeight - vh / 2 : s.scrollY;
        vp.scrollTo(sx, sy, "pointer");
        mode = W.pan;
      }
      vp.release();
    });
    on(ctx.interaction.drag, (e) => {
      const vp = viewportOf(ctx, e.start.section);
      if (e.owner !== name || !vp || !frame(e.start.section)) return;
      if (e.phase >= DragPhase.end) return vp.release();
      const s = vp.live.current;
      const ax = vp.axis;
      if (mode === W.pan && kind === 0 && ax && av) {
        // The navigator's window follows the pointer on the map, through its blocks: its real start moved by the drawn
        // one's travel (a floored window travels less than its share).
        const wd = g[9] - g[8];
        io[0] = av.navX0 + (g[1] > wd ? (e.dx * (g[1] - (av.navX1 - av.navX0))) / (g[1] - wd) : 0);
        navUnits(av, io);
        buf[0] = s.cell;
        buf[1] = io[0] * s.cell + io[1] * ax.gap;
        vp.setViewFrom(buf, "pointer");
      } else if (mode === W.pan) {
        // Per px of the window's travel (a window drawn longer than its share travels less than its span).
        buf[0] = kind === 2 ? 0 : e.dx * g[12] * s.contentWidth;
        buf[1] = kind === 0 ? 0 : e.dy * g[13] * s.contentHeight;
        vp.panByFrom(buf, "pointer");
      } else if ((mode === W.left || mode === W.right) && ax && ax.cells > 0) {
        // The moved end from the press by the pointer's travel, the other held (in content), both clamped: the
        // zoom limit (`maxCell`) as the window's least share of the axis.
        const area = s.width - ax.insetX;
        const gaps = s.contentWidth - ax.cells * s.cell;
        ends[0] = e.totalX;
        ends[1] = Math.min(1, area / (ax.cells * s.maxCell + gaps));
        if (!windowEnds(g0, false, mode === W.left, ends)) {
          buf[0] = at[0];
          buf[1] = at[1];
          return vp.setViewFrom(buf, "pointer");
        }
        // The ends as the axis's units: on the navigator through its map's blocks, on the overview's column linearly;
        // cells u and gaps q, content x = u·cell + q·gap.
        io[0] = g0[0] + ends[0] * g0[1];
        if (kind === 0) navUnits(av as SectionAxisView, io);
        else {
          io[0] = ends[0] * ax.cells;
          io[1] = ends[0] * (ax.count - 1);
        }
        const u0 = io[0];
        const q0 = io[1];
        io[0] = g0[0] + ends[1] * g0[1];
        if (kind === 0) navUnits(av as SectionAxisView, io);
        else {
          io[0] = ends[1] * ax.cells;
          io[1] = ends[1] * (ax.count - 1);
        }
        // At a limit when within rounding of it (an end stopped at the zoom limit is at it exactly); inline, as below.
        const c = (area - (io[1] - q0) * ax.gap) / (io[0] - u0);
        const cell = c >= s.maxCell * (1 - 1e-9) ? s.maxCell : c <= s.minCell * (1 + 1e-9) ? s.minCell : c;
        buf[0] = cell;
        buf[1] = mode === W.left ? io[0] * cell + io[1] * ax.gap - area : u0 * cell + q0 * ax.gap;
        vp.setViewFrom(buf, "pointer");
      } else if ((mode === W.top || mode === W.bottom) && at[4] > 0) {
        // The same for rows: the content's height taken as proportional to the row scale.
        const rh = g0[3];
        ends[0] = e.totalY;
        ends[1] = Math.min(1, (rh * at[2]) / (at[4] * s.maxRowScale));
        if (!windowEnds(g0, true, mode === W.top, ends)) {
          buf[0] = at[2];
          buf[1] = at[3];
          return vp.setRowViewFrom(buf, "pointer");
        }
        const v = (rh * at[2]) / (at[4] * (ends[1] - ends[0]));
        const r = v >= s.maxRowScale * (1 - 1e-9) ? s.maxRowScale : v <= s.minRowScale * (1 + 1e-9) ? s.minRowScale : v;
        const ch = (at[4] * r) / at[2];
        buf[0] = r;
        buf[1] = mode === W.top ? ends[1] * ch - rh : ends[0] * ch;
        vp.setRowViewFrom(buf, "pointer");
      }
    });
    if (kind !== 0) return;
    // The navigator's double click: zoomed in to the focus window or closer, back to the full view; wider than it (the
    // full view, or part way in), in to the focus window at the position.
    on(ctx.interaction.click, (e) => {
      const sec = e.hit.section;
      const p = h.forSection(sec) as NavigatorParams | null;
      if (e.count !== 2 || e.hit.region !== region || e.input.button !== 0 || p?.doubleClick !== "toggle") return;
      const vp = viewportOf(ctx, sec);
      if (!vp || !frame(sec) || windowIdle(vp.live.current, vp.axis)) return;
      const s = vp.live.current;
      // The navigator maps the whole axis: the position's view x, px from the data area's left.
      io[0] = e.input.x;
      navUnits(av as SectionAxisView, io);
      toggleZoom(vp, p, io[0] * s.cell + io[1] * (vp.axis?.gap ?? 0) - s.scrollX);
    });
    on(ctx.interaction.doubleClick, (e) => {
      const p = h.forSection(e.hit.section) as NavigatorParams | null;
      if (e.hit.region !== region || e.input.button !== 0 || p?.doubleClick !== "toggle") return;
      const vp = viewportOf(ctx, e.hit.section);
      if (vp && !windowIdle(vp.live.current, vp.axis)) e.claim(name);
    });
  });
}

/**
 * A double click's zoom (the navigator's and the data's): zoomed in to the focus window (`zoomCell`, else about
 * `zoomPositions` across the data area, at least the letters threshold) or closer, back to the full view; wider (the
 * full view, or part way in), in to the focus window centered at `x`, px from the data area's left. Animated.
 */
function toggleZoom(
  vp: NonNullable<ReturnType<typeof viewportOf>>,
  p: { zoomPositions: number; zoomCell: number },
  x: number,
) {
  const s = vp.live.current;
  const area = s.width - (vp.axis?.insetX ?? 0);
  const cell = Math.min(
    s.maxCell,
    p.zoomCell > 0 ? p.zoomCell : Math.max(vp.props.minLetterCell, area / Math.max(1, p.zoomPositions)),
  );
  if (s.cell >= cell * (1 - 1e-9) && s.cell > s.minCell * (1 + 1e-9)) vp.fit(TOGGLE_MOVE);
  else vp.zoomCentered(cell, x, TOGGLE_MOVE);
}
/** The focus window's cell (`toggleZoom`'s). */
function focusCell(vp: NonNullable<ReturnType<typeof viewportOf>>, p: { zoomPositions: number; zoomCell: number }) {
  const s = vp.live.current;
  const area = s.width - (vp.axis?.insetX ?? 0);
  return Math.min(
    s.maxCell,
    p.zoomCell > 0 ? p.zoomCell : Math.max(vp.props.minLetterCell, area / Math.max(1, p.zoomPositions)),
  );
}
/** The padding around an element a double click focuses (a fraction of the data area's width on each side). */
const ITEM_PAD = 0.05;

export const ScrollbarBehavior = windowBehavior(
  "scrollbar",
  HitRegion.scrollbar,
  2,
  "Dragging the thumb scrolls the rows; a press beside it pages.",
  {},
);
export const NavigatorBehavior = windowBehavior(
  "navigator",
  HitRegion.navigator,
  0,
  "A press jumps; dragging the window pans; dragging an edge zooms; a double click toggles the full view and a zoom there.",
  NAVIGATOR,
);
export const VerticalOverviewBehavior = windowBehavior(
  "verticalOverview",
  HitRegion.vnav,
  1,
  "The window on both axes: top and bottom handles zoom rows, left and right columns; a press jumps.",
  {},
);

const refOf = (ctx: SeqViewerContext, section: number) => ctx.layout.sectionRef(section) ?? "";

/**
 * A split track's chevron, after its label's text (header-toggle.ts: the numbers the composer draws it by): the hit is
 * within its box and a few px around it, only where it is drawn (a scrolled entry whose split has values).
 */
function splitToggle(ctx: SeqViewerContext, view: ViewModel, l: SectionLayout, hit: Hit): boolean {
  const e = hit.entry;
  if (!splitChevron(view, l, e)) return false;
  const rv = view.resolved;
  let depth = 0;
  for (let g = l.entryGroup[e]; g >= 0 && depth < 16; g = rv.groups.parent[g]) depth++;
  const t = rv.tracks.label[hit.track] >= 0 ? rv.tracks.label[hit.track] : rv.tracks.ref[hit.track];
  const px = Math.max(
    8,
    Math.min(ctx.layout.props(refOf(ctx, hit.section)).labelFontPx, Math.floor(l.entryHeight[e] * 0.8)),
  );
  const x = 10 + depth * 14;
  const w = ctx.text.label(view.set.strings.get(t), ctx.text.font("ui", px, 400)).width;
  const cx = x + Math.min(w, l.headerWidth - x - 6 - CHEVRON_ROOM) + CHEVRON_GAP + CHEVRON_HALF;
  return Math.abs(hit.x - cx) <= CHEVRON_HALF + CHEVRON_SLOP;
}

export const HeadersBehavior = defineBehavior(
  "headers",
  "interaction",
  80,
  { activate: "doubleClick" as "click" | "doubleClick" | "none" },
  "A group's header, a split track's chevron or a clade toggles on a single click; activating a track header emits `events.activate`.",
  ({ ctx, h, on }) => {
    on(ctx.interaction.click, (e) => {
      const hit = e.hit;
      const p = h.forSection(hit.section);
      const a = ctx.graph.assembled.value;
      const rv = a?.view?.resolved;
      const l = ctx.layout.get(refOf(ctx, hit.section))?.current;
      if (hit.region !== HitRegion.header || !p || !a || !rv || !l || e.input.button !== 0) return;
      const str = (id: number) => a.set.strings.get(id);
      // A toggle is a single click's: a double click toggles once, not twice.
      const once = e.count === 1;
      if (hit.treeNode >= 0 && l.rows.alignment >= 0) {
        if (once) ctx.ui.toggleClade(str(a.set.alignments[l.rows.alignment]?.key as number), hit.treeNode);
      } else if (hit.entry >= 0 && l.entryKind[hit.entry] === EntryKind.group) {
        if (once) ctx.ui.toggleGroup(str(rv.groups.ref[l.entryGroup[hit.entry]]));
      } else if (
        hit.track >= 0 &&
        rv.tracks.split[hit.track] >= 0 &&
        l.entryRow[hit.entry] < 0 &&
        splitToggle(ctx, a.view as ViewModel, l, hit)
      ) {
        if (once) ctx.ui.toggleGroup(str(rv.tracks.ref[hit.track]));
      } else if (hit.track >= 0 && e.count === (p.activate === "click" ? 1 : p.activate === "doubleClick" ? 2 : 0)) {
        internal(ctx).activate(refOf(ctx, hit.section), str(rv.tracks.ref[hit.track]), "pointer");
      }
    });
    // A track header's double click activates it (`activate: "doubleClick"`).
    on(ctx.interaction.doubleClick, (e) => {
      const p = h.forSection(e.hit.section);
      if (
        e.hit.region === HitRegion.header &&
        e.hit.track >= 0 &&
        e.input.button === 0 &&
        p?.activate === "doubleClick"
      )
        e.claim("headers");
    });
  },
);

const GRANULARITY: Record<string, number> = {
  position: SelectionGranularity.position,
  range: SelectionGranularity.range,
  unit: SelectionGranularity.unit,
  copyGroup: SelectionGranularity.copyGroup,
  column: SelectionGranularity.column,
  row: SelectionGranularity.row,
};
/**
 * A failed streamed region of the hit's track under it: its track's failed regions asked again; a gated one (the zoom
 * gate's "zoom in, or click to load"): asked again whole; true. Edge-rate.
 */
function retried(ctx: SeqViewerContext, hit: Hit): boolean {
  const m = ctx.graph.assembled.value?.view;
  if (!m || hit.track < 0 || hit.block < 0) return false;
  const rv = m.resolved;
  for (let k = 0; k < rv.layers.count; k++) {
    if (rv.layers.track[k] !== hit.track) continue;
    const ref = m.set.strings.get(rv.layers.ref[k]);
    if (ctx.streams.force(ref, hit.block, hit.index)) return true;
    // A rows layer's tile row is the hit's lane's (its lane within the layer's rows).
    const rows = ctx.layout.get(ctx.layout.sectionRef(hit.section) ?? "")?.current.rows;
    const lane = rows && rows.layer === k && hit.lane >= 0 ? hit.lane - rows.laneStart : -1;
    if (ctx.streams.state(ref, hit.block, hit.index, lane) !== RegionState.failed) continue;
    return ctx.actions.run("stream.retry", m.set.strings.get(rv.tracks.ref[hit.track]), "pointer");
  }
  return false;
}

/** The hit's element as loci; null: none. Edge-rate. */
function lociOfHit(hit: Hit): SeqLociAny | null {
  const K = HoverKind;
  // A streamed item (ST4b): its range on the unit (it has no row in the set to select).
  if (hit.kind === K.seq && hit.span[1] > hit.span[0] && hit.seq.set)
    return {
      kind: "seq-loci",
      set: hit.seq.set,
      elements: [{ unit: hit.seq.unit, indices: OrderedSet.ofRange(hit.span[0], hit.span[1]) }],
    } as SeqLociAny;
  const l =
    hit.kind === K.seq
      ? hit.seq
      : hit.kind === K.aligned
        ? hit.aligned
        : hit.kind === K.column
          ? hit.column
          : hit.kind === K.item
            ? hit.itemLoc
            : hit.kind === K.row
              ? hit.row
              : null;
  const loci = l ? (Loci.fromLocation(l) as SeqLociAny) : null;
  // A dot plot's cell: its partner under y too (seqspec:interaction#hover).
  if (!loci || hit.kind !== K.seq || hit.pairUnit < 0) return loci;
  const pair = {
    kind: "seq-loci",
    set: hit.seq.set,
    elements: [{ unit: hit.pairUnit, indices: OrderedSet.ofSingleton(hit.pairIndex) }],
  };
  return Loci.union(loci as never, pair as never) as SeqLociAny;
}

/** Plain: the primary button with no modifier (or only `also`'s). */
const plain = (button: number, mods: number, also?: Trigger) =>
  button === 0 && (mods === 0 || (also !== undefined && Binding.matches(also, button, mods)));
const claims = (
  p: { dragRegions: readonly string[]; dragOn: Binding },
  start: Hit,
  button: number,
  mods: number,
  also?: Trigger,
) =>
  (p.dragRegions.includes(RegionName[start.region]) && plain(button, mods, also)) ||
  (start.region === HitRegion.data && Binding.match(p.dragOn, button, mods));

export const SelectBehavior = defineBehavior(
  "select",
  "interaction",
  50,
  {
    click: [{ button: "primary" }] as Binding,
    extend: [{ button: "primary", shift: true }] as Binding,
    /** Adds the element, or removes it when selected. */
    toggle: [
      { button: "primary", ctrl: true },
      { button: "primary", meta: true },
    ] as Binding,
    /**
     * Shift with a click or a drag on positions or columns (seqspec:interaction#selection): `extend`, the range from the
     * anchor; `toggle`, adds or removes (a drag decided at its press: on a selected position it removes its range).
     */
    shiftPositions: "extend" as "extend" | "toggle",
    /** `toggle` on an item or a row: Shift-click adds or removes it; `set`: selects it alone. */
    shiftItems: "toggle" as "toggle" | "set",
    /**
     * Where a double click at the point means something else (a behavior claims `doubleClick`), a single click's change
     * waits for the double-click time, so a double click never touches the selection; false: it applies at once.
     */
    doubleClickWait: true,
    /** A double or triple click: `ignore` (its first click did all), or `click`, that click again (the text presets). */
    repeatClicks: "ignore" as "ignore" | "click",
    clearOnEmpty: true,
    toggleWhole: true,
    /**
     * What the element is closed by before it is widened (seqspec:interaction#channels); a click selects what the hover
     * lit.
     */
    closure: "containing" as Closure,
    granularity: "position" as "position" | "range" | "unit" | "copyGroup" | "column" | "row",
    dragRegions: ["column"] as readonly RegionName[],
    dragOn: [{ button: "primary", alt: true }] as Binding,
    /** A drag on data over whole columns (positions on a blocks axis), as `shiftPositions` says. */
    dragColumns: [{ button: "primary", shift: true }] as Binding,
    /**
     * Past an edge during a select drag the view scrolls (text and rows in a box section down, a linear section across),
     * at `autoscrollMin` px a second plus `autoscrollRate` for each px past, up to `autoscrollMax`.
     */
    autoscroll: true,
    autoscrollMin: 120,
    autoscrollRate: 15,
    autoscrollMax: 1440,
  },
  "A click selects (again: clears a whole selection), held for the double-click time where a double click means something else (a behavior claims `doubleClick`: a zoom), which then never selects; Shift (a click or a drag) extends a range of positions or columns, and adds or removes an item or a row; Ctrl/⌘ adds or removes; a drag it claims selects a range, its end following the point (clamped to the surface) and the view scrolling past an edge. In text, a drag follows the pointer line by line.",
  ({ ctx, h, on, add }) => {
    const sel = ctx.selection;
    /** The element closed by the section's projector, then widened (a click and a release alike). */
    const closed = (loci: SeqLociAny, p: { closure: Closure; granularity: string }) =>
      sel.widen(sel.close(loci, p.closure), GRANULARITY[p.granularity] ?? 0);
    const layoutOf = (section: number) => ctx.layout.get(refOf(ctx, section))?.current ?? null;
    /**
     * Where Shift extends from in text: the last click's position or drag's start, as a position of the model (loci of
     * one position), so another section showing it (another view) finds it through its own projection; null: none.
     */
    let anchor: SeqLociAny | null = null;
    /** Block `b`'s shown position of the anchor in section `vp`'s projection; -1: not on that block. */
    const anchorOn = (vp: NonNullable<ReturnType<typeof viewportOf>>, b: number) => {
      const ax = vp.axis;
      const hull = anchor && ax ? vp.hullOf(anchor) : null;
      for (let q = 0; hull && ax && q < hull.length; q += 3) if (hull[q] === b) return ax.shownOf(b, hull[q + 1]);
      return -1;
    };
    /** The block of the section the anchor is on, other than `b`; -1: none. */
    const anchorBlock = (vp: NonNullable<ReturnType<typeof viewportOf>>, b: number) => {
      const hull = anchor ? vp.hullOf(anchor) : null;
      for (let q = 0; hull && q < hull.length; q += 3) if (hull[q] !== b) return hull[q];
      return -1;
    };
    /** The anchor at shown k of block b of a section. */
    const anchorAt = (section: number, b: number, k: number) => {
      const ax = viewportOf(ctx, section)?.axis;
      anchor = ax
        ? shownLoci(ctx, section, b, k, k + 1, ax.alignment[b] >= 0 ? HoverKind.column : HoverKind.seq)
        : null;
    };
    /** The selection's runs per section as the overlay lights them, made at the first Shift-click. */
    let channels: ChannelMarks | null = null;
    add(() => channels?.dispose());
    /** A run of shown positions' scratch (`runToShown`), grown when a selection needs more. */
    let shownRuns = new Int32Array(16);
    /**
     * Shift-click on a position or a column (seqspec:interaction#selection): the selection on the clicked block
     * (positions, items or columns, as the overlay lights it from the model's loci: each run as shown, a cyclic origin's
     * two parts, a filter's pieces) extended to the clicked position, one end kept: with several runs there and the
     * anchor on the block, the anchor; else outside it its far end, inside it the anchor where it is one of its ends (a
     * drag's start, a click's position), else the farther end. With no selection on the block: from the anchor's block's
     * end nearest it, on another block of the section; else false (a click).
     */
    const extendRange = (hit: Hit, p: { closure: Closure; granularity: string }): boolean => {
      const vp = viewportOf(ctx, hit.section);
      const ax = vp?.axis;
      if (!vp || !ax) return false;
      const b = hit.block;
      const k = hit.shown;
      // The selection's runs on the block as lit (the overlay's channel: positions, items, columns, in block indices; a
      // cyclic origin's arc in two), as shown runs, merged where they touch. Edge-rate: a click.
      const found: [number, number][] = [];
      channels ??= new ChannelMarks(ctx);
      const lit = channels.marks(hit.section)?.runs(MarkBit.selection);
      for (let r = lit ? lit.from(b) : 0; lit && r < lit.to(b); r++) {
        let n = ax.runToShown(b, lit.start[r], lit.stop[r], shownRuns);
        if (2 * n > shownRuns.length) {
          shownRuns = new Int32Array(4 * n);
          n = ax.runToShown(b, lit.start[r], lit.stop[r], shownRuns);
        }
        for (let j = 0; j < n; j++)
          if (shownRuns[2 * j + 1] > shownRuns[2 * j]) found.push([shownRuns[2 * j], shownRuns[2 * j + 1]]);
      }
      found.sort((x, y) => x[0] - y[0]);
      let runs = 0;
      let end = -1;
      for (const [s0, s1] of found) {
        if (s0 > end) runs++;
        end = Math.max(end, s1);
      }
      const lo = found.length > 0 ? found[0][0] : -1;
      const hi = end - 1;
      let k0: number;
      let k1: number;
      const other = lo < 0 ? anchorBlock(vp, b) : -1;
      if (lo >= 0) {
        const a = anchorOn(vp, b);
        const far = Math.abs(k - lo) >= Math.abs(k - hi) ? lo : hi;
        const fixed = runs > 1 && a >= 0 ? a : k < lo ? hi : k > hi ? lo : a === lo || a === hi ? a : far;
        k0 = Math.min(fixed, k);
        k1 = Math.max(fixed, k) + 1;
        anchorAt(hit.section, b, fixed);
      } else if (other >= 0) {
        // Another block's position: from that block's end nearest the anchor's, as a drag stops there.
        k0 = other < b ? 0 : k;
        k1 = other < b ? k + 1 : ax.shown[b];
      } else return false;
      const raw = shownLoci(ctx, hit.section, b, k0, k1, ax.alignment[b] >= 0 ? HoverKind.column : HoverKind.seq);
      if (raw) sel.set(closed(raw, p), "pointer");
      return true;
    };
    /**
     * A single click's change held while a double click there would mean something else (`doubleClick`): applied when
     * the double-click time passes with no second click, or at once at the next click elsewhere; dropped by a double
     * click, which so never touches the selection.
     */
    const held = { op: 0, loci: null as SeqLociAny | null, timer: -1 };
    const OP_SET = 1;
    const OP_REMOVE = 2;
    const OP_TOGGLE = 3;
    const OP_CLEAR = 4;
    const apply = (op: number, w: SeqLociAny | null) => {
      if (op === OP_CLEAR) sel.clear("pointer");
      else if (!w) return;
      else if (op === OP_SET) sel.set(w, "pointer");
      else if (op === OP_REMOVE) sel.remove(w, "pointer");
      else if (op === OP_TOGGLE) sel.toggle(w, "pointer");
    };
    const drop = () => {
      if (held.timer >= 0) ctx.scheduler.clearTimer(held.timer);
      held.timer = -1;
      held.op = 0;
      held.loci = null;
    };
    let applying = false;
    const flush = () => {
      const { op, loci } = held;
      drop();
      applying = true;
      try {
        // Its loci on the set shown now (a patch or an arrival may have assembled another meanwhile).
        const set = ctx.graph.assembled.value?.set;
        if (op) apply(op, loci && set ? carry(loci, set) : loci);
      } finally {
        applying = false;
      }
    };
    // A key meanwhile: Escape drops the held click, another applies it first (it acts on what the click did).
    on(ctx.interaction.key, (e) => {
      if (!held.op) return;
      if (e.input.key === "Escape") drop();
      else flush();
    });
    // The selection changed meanwhile (a host, a command): the held click is moot. A carry to a new set (origin
    // `graph`) changed nothing anyone did (seqspec:interaction#replays-and-carries).
    add(
      sel.current.subscribe(() => {
        if (held.op && !applying && sel.current.origin !== "graph") drop();
      }),
    );
    const settled = () => {
      held.timer = -1;
      flush();
    };
    on(ctx.interaction.click, (e) => {
      // A second click: a double click drops what its first click held; any other click applies it first.
      if (held.op) {
        if (e.count > 1) return drop();
        flush();
      }
      const p = h.forSection(e.hit.section);
      if (!p || (e.hit.region !== HitRegion.data && e.hit.region !== HitRegion.column)) return;
      // A click on a failed streamed region (its error band) retries it, on a gated one loads it whole, and selects
      // nothing (ST4b).
      if (e.input.button === 0 && e.count === 1 && retried(ctx, e.hit)) return;
      const extend = Binding.match(p.extend, e.input.button, e.input.modifiers);
      const toggle = Binding.match(p.toggle, e.input.button, e.input.modifiers);
      const set = Binding.match(p.click, e.input.button, e.input.modifiers);
      const hit = e.hit;
      // A double or triple click selects nothing more (`repeatClicks: "ignore"`): its first click did, or held it
      // (above); `click`: that click again, which changes nothing (text, where nothing takes it).
      if (e.count > 1 && p.repeatClicks !== "click") return;
      // Shift on a position or a column extends the range (seqspec:interaction#selection), on an item or a row adds or
      // removes it, by the params.
      const onRange = hit.kind === HoverKind.seq || hit.kind === HoverKind.column;
      const extendRangeHere = extend && onRange && p.shiftPositions === "extend";
      if (extendRangeHere && hit.block >= 0 && hit.shown >= 0 && extendRange(hit, p)) return;
      const loci = extend || toggle || set ? lociOfHit(hit) : null;
      let op = 0;
      let w: SeqLociAny | null = null;
      if (!loci) op = set && p.clearOnEmpty ? OP_CLEAR : 0;
      else {
        w = closed(loci, p);
        if (hit.block >= 0 && hit.shown >= 0 && (set || extend)) anchorAt(hit.section, hit.block, hit.shown);
        const shiftToggles = extend && (onRange ? p.shiftPositions === "toggle" : p.shiftItems === "toggle");
        if (toggle || shiftToggles) op = OP_TOGGLE;
        else if (p.toggleWhole && e.count === 1 && !extend && sel.isWhole(w)) op = OP_REMOVE;
        else if (e.count > 1 && sel.isWhole(w)) op = 0;
        else op = OP_SET;
      }
      if (!op) return;
      // Where a double click here means something else (a zoom), the change waits for the double-click time.
      if (p.doubleClickWait && e.count === 1 && internal(ctx).input.askDoubleClick(e) !== "") {
        held.op = op;
        held.loci = w;
        held.timer = ctx.scheduler.setTimer(settled, ctx.config.get(SeqViewerConfig.Input.DoubleClickMs));
        return;
      }
      apply(op, w);
    });
    // Removed with a click held: nothing applies later.
    add(drop);
    const pd = sel.pending.current;
    const bump = (sel.pending as { bump?: () => void }).bump as () => void;
    let adding = false;
    /** The drag's change: set, add or remove its range. */
    let dragOp = 0;
    // A drag's end follows the point (seqspec:interaction#selection), clamped to its surface during the press: in text
    // its nearest position line by line, past a line's end or between lines included; on a linear section the column
    // under its x (and on rows the row under its y), the edge's while the point is past it. Past an edge the view
    // scrolls, faster the further past (`overX`, `overY`): text and rows vertically in a box section, a linear section
    // across. The point, its unclamped x and y, its nearest line, block and shown, the scroll's numbers and its frame.
    let wrapped = false;
    /** Where the drag is: wrapped, linear (by its point), or on rings (by its hit). */
    let geom = 0;
    const at = new Float64Array(4);
    const near = new Int32Array(4);
    const buf = new Float64Array(2);
    let frame = -1;
    let scrolled = false;
    /** The scroll's last frame time (ms; -1: its first frame comes next), and the px it owes down and across. */
    const time = new Float64Array([-1, 0, 0]);
    /** The scroll's speeds, px a second across and down (`speeds`). */
    const vel = new Float64Array(2);
    /**
     * The press's position on either side of the point: on a position both are it; between two (a gap, past a line's
     * end, an inline label), the one before and the one after, the drag's direction picking which starts the range.
     */
    const startAt = new Int32Array(2);
    /** The point at `at` as the pending end: in text its nearest position on the pending block (else that block's
     * first or last); on a linear section the column at its x on the pending block, and on rows the row at its y. */
    const follow = () => {
      const l = layoutOf(pd.section);
      const v = ctx.layout.axis(refOf(ctx, pd.section));
      const ax = v?.axis;
      if (!l || !v || !ax) return;
      const b = pd.block;
      let k: number;
      let s = pd.startShown;
      let row = pd.rowEnd;
      if (geom === G_WRAPPED) {
        if (!wrappedPositionNear(l, v, at, near)) return;
        // Between two positions the end takes the one toward the start (forward the one before, backward the one after).
        k = near[1] === b ? near[2] : near[1] > b ? ax.shown[b] - 1 : 0;
        if (near[1] === b && near[3] === 1 && k < startAt[0]) k = Math.min(k + 1, ax.shown[b] - 1);
        s = k >= startAt[1] ? startAt[1] : startAt[0];
      } else {
        const x = at[0] < l.dataLeft ? l.dataLeft : at[0] >= l.dataRight ? l.dataRight - 0.5 : at[0];
        // The cell drawn under the point, as the hit names it (GG-01).
        k = v.drawnShown(b, Math.floor((x - v.blockX[b]) / (v.span[b] * v.cell)), Math.floor(x * v.dpr) | 0);
        k = k < 0 ? 0 : k >= ax.shown[b] ? ax.shown[b] - 1 : k;
        const r = l.rows;
        if (pd.kind === HoverKind.aligned && r.layer >= 0 && r.count > 0) {
          const d = Math.floor((at[1] - (l.scrolledTop + l.entryTop[r.entry] - v.scrollY)) / r.height);
          row = d < 0 ? 0 : d >= r.count ? r.count - 1 : d;
        }
      }
      if (k === pd.endShown && s === pd.startShown && row === pd.rowEnd) return;
      pd.endShown = k;
      pd.startShown = s;
      pd.rowEnd = row;
      bump();
    };
    /**
     * Px a second into `vel`, by the section's params (`autoscroll…`; off: 0), faster the further past an edge: across,
     * a linear section's pan while the point is left or right of its data area; down, the rows' while it is above or
     * below them (text or rows in a box section). Written, not returned: a double a call V8 does not inline returns or
     * takes is boxed (about 40 bytes an autoscroll frame on Node 24.19, 2026-10-14).
     */
    const speeds = () => {
      vel[0] = 0;
      vel[1] = 0;
      const l = layoutOf(pd.section);
      const p = h.forSection(pd.section);
      if (!l || !p?.autoscroll) return;
      if (geom === G_LINEAR) {
        const x = at[2];
        const past = x < l.dataLeft ? l.dataLeft - x : x > l.dataRight ? x - l.dataRight : 0;
        const v = Math.min(p.autoscrollMax, p.autoscrollMin + past * p.autoscrollRate);
        vel[0] = past <= 0 ? 0 : x < l.dataLeft ? -v : v;
      }
      if (l.heightMode === BOX && l.vscroll && (geom === G_WRAPPED || pd.kind === HoverKind.aligned)) {
        const y = at[3];
        const past = y < l.scrolledTop ? l.scrolledTop - y : y > l.height ? y - l.height : 0;
        const v = Math.min(p.autoscrollMax, p.autoscrollMin + past * p.autoscrollRate);
        vel[1] = past <= 0 ? 0 : y < l.scrolledTop ? -v : v;
      }
    };
    /**
     * A frame of the scroll: by the time since the last (at most `AUTOSCROLL_DT`), so its speed is the display's rate's
     * own; it stops when the point is back over the data, or the view is at its end (a move asks again).
     */
    const tick = (t: number) => {
      frame = -1;
      const vp = pd.active ? viewportOf(ctx, pd.section) : null;
      if (vp) speeds();
      const vy = vp ? vel[1] : 0;
      const vx = vp ? vel[0] : 0;
      if (!vp || (vx === 0 && vy === 0)) {
        time[0] = -1;
        return;
      }
      const dt = time[0] < 0 ? 0 : Math.max(0, Math.min(AUTOSCROLL_DT, t - time[0]));
      time[0] = t;
      // With what the snapped scroll left over before: the distance is the time's, at any frame rate.
      buf[0] = time[2] + (vx * dt) / 1000;
      buf[1] = time[1] + (vy * dt) / 1000;
      const x0 = vp.live.current.scrollX;
      const y0 = vp.live.current.scrollY;
      vp.panByFrom(buf, "pointer");
      const movedX = vp.live.current.scrollX - x0;
      const movedY = vp.live.current.scrollY - y0;
      time[2] = vx === 0 ? 0 : buf[0] - movedX;
      time[1] = vy === 0 ? 0 : buf[1] - movedY;
      if (movedX === 0 && movedY === 0 && Math.abs(buf[0]) + Math.abs(buf[1]) >= 1) {
        time[0] = -1;
        return;
      }
      scrolled = true;
      follow();
      frame = ctx.scheduler.requestFrame(tick);
    };
    /** A drag past an edge starts the scroll, unless it runs. */
    const scroll = () => {
      if (!pd.active || geom === G_RINGS || frame >= 0) return;
      speeds();
      if (vel[0] === 0 && vel[1] === 0) return;
      // Its first frame counts from now, at any display rate.
      if (time[0] < 0) {
        time[0] = ctx.scheduler.now();
        time[1] = time[2] = 0;
      }
      frame = ctx.scheduler.requestFrame(tick);
    };
    /** The scroll paused: its frame cancelled (the drag goes on; a move starts it again). */
    const pause = () => {
      if (frame >= 0) ctx.scheduler.cancelFrame(frame);
      frame = -1;
      time[0] = -1;
    };
    /** The drag is over: its scroll's frame cancelled, the scroll committed. */
    const stop = () => {
      pause();
      if (scrolled) viewportOf(ctx, pd.section)?.release();
      scrolled = wrapped = false;
    };
    // The drag's surface lost focus (another window, a tab hidden): the view stops scrolling.
    on(ctx.interaction.focus, (e) => {
      if (!e.focused && pd.active && e.section === pd.section) pause();
    });
    /**
     * The pending end at the drag event's point: the dispatcher's first drag event comes at a move past the slop and its
     * last at the release, either one far from the move before it (a fast pointer), so each is followed, not only the
     * moves between. On rings, the hit's position on the pending block.
     */
    const extend = (e: DragEvt) => {
      if (geom !== G_RINGS) {
        at[0] = e.input.x;
        at[1] = e.input.y;
        at[2] = e.input.x + e.input.overX;
        at[3] = e.input.y + e.input.overY;
        return follow();
      }
      const hit = e.hit;
      if (hit.section !== pd.section || hit.block !== pd.block || hit.shown < 0) return;
      const row = pd.kind === HoverKind.aligned && hit.displayRow >= 0 ? hit.displayRow : pd.rowEnd;
      if (hit.shown === pd.endShown && row === pd.rowEnd) return;
      pd.endShown = hit.shown;
      pd.rowEnd = row;
      bump();
    };
    on(ctx.interaction.drag, (e) => {
      const s = e.start;
      if (e.phase === DragPhase.start) {
        const p = h.forSection(s.section);
        if (e.owner !== "" || !p) return;
        const columns = s.region === HitRegion.data && Binding.match(p.dragColumns, e.button, e.input.modifiers);
        if (!columns && !claims(p, s, e.button, e.input.modifiers)) return;
        // Wrapped, a press on the text's data area or its inset starts at its nearest position (a gap, a line's end, a
        // labels band, the inset around the text).
        const l = s.region === HitRegion.data || s.region === HitRegion.gutter ? layoutOf(s.section) : null;
        const v = l?.geometry === WRAPPED ? ctx.layout.axis(refOf(ctx, s.section)) : null;
        at[0] = at[2] = s.x;
        at[1] = at[3] = s.y;
        const w = l !== null && v !== null && wrappedPositionNear(l, v, at, near);
        if (!w && (s.block < 0 || s.shown < 0)) return;
        // A press that drags is no double click: a click held before it applies first.
        flush();
        e.claim("select");
        adding = columns;
        dragOp = DRAG_SET;
        wrapped = w;
        geom = w
          ? G_WRAPPED
          : l?.geometry === CIRCULAR || layoutOf(s.section)?.geometry === CIRCULAR
            ? G_RINGS
            : G_LINEAR;
        pd.active = true;
        pd.section = s.section;
        pd.block = w ? near[1] : s.block;
        const onColumns = (viewportOf(ctx, s.section)?.axis?.alignment[pd.block] ?? -1) >= 0;
        pd.kind =
          columns || w
            ? onColumns
              ? HoverKind.column
              : HoverKind.seq
            : s.displayRow >= 0
              ? HoverKind.aligned
              : s.kind === HoverKind.column
                ? HoverKind.column
                : HoverKind.seq;
        pd.startShown = pd.endShown = w ? near[2] : s.shown;
        pd.rowStart = pd.rowEnd = w ? -1 : s.displayRow;
        startAt[0] = pd.startShown;
        startAt[1] = w && near[3] === 1 ? Math.min(near[2] + 1, (v?.axis?.shown[pd.block] ?? 1) - 1) : pd.startShown;
        // Shift with a drag as with a click (`shiftPositions`): `extend` from the anchor (on this block) to the drag's
        // end; `toggle` decided at the press: on a selected position the drag removes its range, else adds it.
        pd.remove = false;
        if (columns) {
          const svp = viewportOf(ctx, s.section);
          const sax = svp?.axis;
          if (p.shiftPositions === "extend") {
            const a = svp ? anchorOn(svp, pd.block) : -1;
            dragOp = DRAG_SET;
            adding = false;
            if (a >= 0) pd.startShown = pd.endShown = startAt[0] = startAt[1] = a;
          } else if (sax) {
            channels ??= new ChannelMarks(ctx);
            const lit = channels.marks(s.section)?.runs(MarkBit.selection);
            pd.remove = lit?.has(pd.block, sax.indexOf(pd.block, pd.startShown as never)) === true;
            dragOp = pd.remove ? DRAG_REMOVE : DRAG_ADD;
          }
        }
        bump();
        extend(e);
        scroll();
      } else if (e.owner === "select" && e.phase === DragPhase.move) {
        extend(e);
        scroll();
      } else if (e.owner === "select") {
        if (e.phase === DragPhase.end) extend(e);
        stop();
        const p = h.forSection(pd.section);
        const raw = e.phase === DragPhase.end && p ? rangeLoci(ctx) : null;
        pd.active = false;
        bump();
        const loci = raw && p ? closed(raw, p) : null;
        pd.remove = false;
        if (loci && dragOp === DRAG_REMOVE) sel.remove(loci, "pointer");
        else if (loci && adding) sel.add(loci, "pointer");
        else if (loci) sel.set(loci, "pointer");
        // Shift-click extends from the range's start, as the drag left it.
        if (loci) anchorAt(pd.section, pd.block, pd.startShown);
      }
    });
    // Removed mid-drag: nothing stays pending, nothing scrolls.
    add(() => {
      stop();
      if (!pd.active) return;
      pd.active = false;
      bump();
    });
  },
);

/** The most ms a frame counts (a frame after a pause never jumps). */
const AUTOSCROLL_DT = 50;
// Read once: a member read per pointer event allocated under the test runner's module wrapper.
const BOX = Height.box;
const WRAPPED = AxisGeometry.wrapped;
const CIRCULAR = AxisGeometry.circular;
/** A select drag's change to the selection. */
const DRAG_SET = 0;
const DRAG_ADD = 1;
const DRAG_REMOVE = 2;
/** A select drag's geometry: text, a linear section, rings. */
const G_WRAPPED = 0;
const G_LINEAR = 1;
const G_RINGS = 2;

/** A drag's block runs (pairs), grown only when a drag over a filtered, reversed or broken piece needs more. */
let runBuffer = new Int32Array(16);
/**
 * The pending range as loci: columns, a band of cells over rows (one rectangle), or a unit's positions. Edge-rate, and
 * O(runs) however long the drag: the shown range goes to block runs (`shownToRuns`, several for a reversed or
 * filtered piece), a band's rows from its display rows' lanes.
 */
function rangeLoci(ctx: SeqViewerContext): SeqLociAny | null {
  const pd = ctx.selection.pending.current;
  const k0 = Math.min(pd.startShown, pd.endShown);
  const k1 = Math.max(pd.startShown, pd.endShown) + 1;
  return shownLoci(ctx, pd.section, pd.block, k0, k1, pd.kind);
}
/** Shown [k0, k1) of block b as loci of `kind` (`HoverKind`: column, seq, or aligned over the pending rows). Edge-rate. */
function shownLoci(ctx: SeqViewerContext, section: number, b: number, k0: number, k1: number, kind: number) {
  const pd = ctx.selection.pending.current;
  const ax = viewportOf(ctx, section)?.axis;
  const model = ctx.graph.assembled.value?.view;
  if (!ax || !model) return null;
  let n = ax.shownToRuns(b, k0, k1, runBuffer);
  if (2 * n > runBuffer.length) {
    runBuffer = new Int32Array(2 * n);
    n = ax.shownToRuns(b, k0, k1, runBuffer);
  }
  const columns = ofRuns(runBuffer, n);
  const set = ax.set;
  const alignment = ax.alignment[b];
  if (kind === HoverKind.column) return { kind: "column-loci", set, alignment, columns } as unknown as SeqLociAny;
  if (kind === HoverKind.seq)
    return ax.unit[b] < 0
      ? null
      : ({ kind: "seq-loci", set, elements: [{ unit: ax.unit[b], indices: columns }] } as unknown as SeqLociAny);
  const r = ctx.layout.get(refOf(ctx, section))?.current.rows;
  const [d0, d1] = [Math.min(pd.rowStart, pd.rowEnd), Math.max(pd.rowStart, pd.rowEnd)];
  const found = new Int32Array(r ? d1 - d0 + 1 : 0);
  let m = 0;
  for (let d = d0; r && d <= d1; d++) {
    const lane = SectionLayouts.laneOf(r, d);
    if (lane >= 0) found[m++] = Lanes.rowOf(model.resolved.lanes, lane);
  }
  const sorted = found.subarray(0, m).sort();
  let u = 0;
  for (let i = 0; i < m; i++) if (i === 0 || sorted[i] !== sorted[i - 1]) sorted[u++] = sorted[i];
  const elements = { rows: OrderedSet.ofSortedArray(sorted.slice(0, u)), columns };
  return u > 0 ? ({ kind: "aligned-loci", set, alignment, elements } as unknown as SeqLociAny) : null;
}

export const DragPanBehavior = defineBehavior(
  "dragPan",
  "interaction",
  10,
  {
    dragRegions: ["data", "header", "gutter"] as readonly RegionName[],
    dragOn: [] as Binding,
    x: true,
    y: true,
    /** None by default: Shift-drag on data selects columns (`select`'s `dragColumns`). */
    dominantAxis: null as Trigger | null,
  },
  "Dragging pans both axes, each with a lock; with `dominantAxis` held (none by default), the dominant axis only.",
  ({ ctx, h, on }) => {
    on(ctx.interaction.drag, (e) => {
      const p = h.forSection(e.start.section);
      if (!p) return;
      if (e.phase === DragPhase.start) {
        if (e.owner !== "" || !claims(p, e.start, e.button, e.input.modifiers, p.dominantAxis ?? undefined)) return;
        e.claim("dragPan");
      } else if (e.owner !== "dragPan") return;
      const vp = viewportOf(ctx, e.start.section);
      if (!vp) return;
      if (e.phase >= DragPhase.end) return vp.release();
      let x = p.x;
      let y = p.y;
      if (x && y && p.dominantAxis !== null && Binding.matches(p.dominantAxis, e.button, e.input.modifiers)) {
        if (Math.abs(e.totalX) >= Math.abs(e.totalY)) y = false;
        else x = false;
      }
      vp.panDrag(e, x, y);
    });
  },
);

export const WheelBehavior = defineBehavior(
  "wheel",
  "interaction",
  10,
  { zoom: [{ ctrl: true }, { meta: true }] as Binding, zoomSpeed: 1, panSpeed: 1 },
  "Zoom at the cursor with Ctrl/⌘ and pinch; a horizontal wheel pans, and so does a vertical one with Shift; a vertical wheel scrolls a box section's rows.",
  ({ ctx, h, on }) => {
    /** The numbers a viewport verb reads. */
    const buf = new Float64Array(2);
    on(ctx.interaction.wheel, (e) => {
      const w = e.input;
      const p = h.forSection(e.hit.section);
      const vp = viewportOf(ctx, e.hit.section);
      const ref = refOf(ctx, e.hit.section);
      const v = ctx.layout.axis(ref);
      const l = ctx.layout.get(ref)?.current;
      if (e.handled || !p || !vp || !v || !l || e.hit.region === HitRegion.none) return;
      const s = vp.live.current;
      if (Binding.match(p.zoom, -1, w.modifiers)) {
        buf[0] = Math.exp(-w.deltaY * 0.004 * p.zoomSpeed);
        buf[1] = w.x - v.left;
        vp.zoomAtFrom(buf, "pointer");
        e.handled = true;
        return;
      }
      // With Shift a vertical wheel pans, by whichever delta the platform set (macOS swaps them under Shift).
      const shift = w.modifiers === Modifier.shift;
      const across = shift && w.deltaX === 0 ? w.deltaY : w.deltaX;
      if (l.geometry === AxisGeometry.linear && Math.abs(across) > (shift ? 0 : Math.abs(w.deltaY))) {
        buf[0] = across * p.panSpeed;
        buf[1] = 0;
        vp.panByFrom(buf, "pointer");
        e.handled = true;
        return;
      }
      // A box section's rows scroll while they can; past the end, and in a content section, the page scrolls.
      const max = s.contentHeight - (l.height - l.scrolledTop);
      const dy = w.deltaY;
      if (l.heightMode === Height.box && ((dy < 0 && s.scrollY > 0) || (dy > 0 && s.scrollY < max))) {
        buf[0] = 0;
        buf[1] = dy * p.panSpeed;
        vp.panByFrom(buf, "pointer");
        e.handled = true;
      }
    });
  },
);

export const DoubleClickZoomBehavior = defineBehavior(
  "doubleClickZoom",
  "interaction",
  40,
  { doubleClick: "toggle" as "toggle" | "none", zoomPositions: 60, zoomCell: 0 },
  "A double click on an element shows its whole extent, then the focus window there, then the full view; on the data elsewhere, the focus window there, or zoomed in to it or closer, the full view.",
  ({ ctx, h, on }) => {
    /** The section the last "focus window at the position" step zoomed (-1: none); another double click there fits. */
    const atPosition = new Int32Array([-1]);
    on(ctx.interaction.click, (e) => {
      const p = h.forSection(e.hit.section);
      const i = e.input;
      if (p?.doubleClick !== "toggle" || e.count !== 2 || i.button !== 0 || i.modifiers !== 0) return;
      if (e.hit.region !== HitRegion.data) return;
      const ref = refOf(ctx, e.hit.section);
      const v = ctx.layout.axis(ref);
      const vp = viewportOf(ctx, e.hit.section);
      if (!v || !vp || ctx.layout.get(ref)?.current.geometry !== AxisGeometry.linear) return;
      const x = i.x - v.left;
      // On an element wider than a position (a range, an item): first its whole extent; shown exactly already, the
      // position's focus window there; at that, the full view.
      const loci = e.hit.kind === HoverKind.item ? lociOfHit(e.hit) : null;
      if (loci) {
        const cell = focusCell(vp, p);
        // Still at the focus window the last step zoomed to (not zoomed since): the full view.
        if (atPosition[0] === e.hit.section && Math.abs(vp.live.current.cell - cell) <= cell * 1e-3) {
          atPosition[0] = -1;
          return vp.fit(TOGGLE_MOVE);
        }
        atPosition[0] = -1;
        if (!vp.showsFocus(loci, ITEM_PAD)) return void vp.focusLoci(loci, { ...TOGGLE_MOVE, pad: ITEM_PAD });
        // Shown exactly already: the focus window at the position.
        atPosition[0] = e.hit.section;
        return vp.zoomCentered(cell, x, TOGGLE_MOVE);
      }
      atPosition[0] = -1;
      toggleZoom(vp, p, x);
    });
    // Its double click at a plain primary click on a linear section's data.
    on(ctx.interaction.doubleClick, (e) => {
      const p = h.forSection(e.hit.section);
      const i = e.input;
      if (p?.doubleClick !== "toggle" || i.button !== 0 || i.modifiers !== 0 || e.hit.region !== HitRegion.data) return;
      if (ctx.layout.get(refOf(ctx, e.hit.section))?.current.geometry === AxisGeometry.linear)
        e.claim("doubleClickZoom");
    });
  },
);

export const RowZoomBehavior = defineBehavior(
  "rowZoom",
  "interaction",
  20,
  { zoom: [{ alt: true }] as Binding, speed: 1 },
  "Zooms rows vertically, every row at once.",
  ({ ctx, h, on }) => {
    /** The numbers a viewport verb reads. */
    const buf = new Float64Array(2);
    on(ctx.interaction.wheel, (e) => {
      const w = e.input;
      const p = h.forSection(e.hit.section);
      const vp = viewportOf(ctx, e.hit.section);
      const l = ctx.layout.get(refOf(ctx, e.hit.section))?.current;
      if (e.handled || !p || !vp || !l || l.geometry !== AxisGeometry.linear || !Binding.match(p.zoom, -1, w.modifiers))
        return;
      const f = Math.exp(-(w.deltaY || w.deltaX) * 0.004 * p.speed);
      buf[0] = vp.live.current.rowScale * f;
      buf[1] = w.y - l.scrolledTop;
      vp.setRowScaleFrom(buf);
      e.handled = true;
    });
  },
);
