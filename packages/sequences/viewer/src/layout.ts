/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The layout (seqspec:viewer#state): a section's entries, sizes, data columns, navigator and overview
// geometry and wrapped lines, headless and rebuilt lazily on the first read after an input changes (the view model,
// `ui`, the section's render props, the viewport's size and row scale, `text.version`). The manager also builds each
// section's axis (axis.ts) and hands it to the viewports (`Viewports.update`), and answers hit tests (hit.ts). A
// section's layout is kept across a removal and return within one update, as its viewport is.

import {
  followsTree,
  Height,
  Lanes,
  Layout,
  type ResolvedLanes,
  type ResolvedView,
  RowSizing,
  Sizing,
} from "@molstar/seq-core/format";
import { type Alphabet, Numbering, type SeqSet, type SortedArray } from "@molstar/seq-core/model";
import { Axis, AxisGeometry, SectionAxisView, type ShownIndex, sameList } from "./axis.js";
import { type ConfigManager, SeqViewerConfig } from "./config.js";
import type { ViewModel } from "./graph.js";
import { Hit, hitTest, type StreamLookup, sourceAt } from "./hit.js";
import { type NavigatorEntry, navigatorOf } from "./navigators.js";
import type { PresetEntry } from "./presets.js";
import type { ColumnRenderer, LayerRenderer, SeqViewerRegistry } from "./registry.js";
import {
  GROUP_LABEL_PAD,
  type GroupLabeler,
  groupLabelBand,
  groupLabelPx,
  type LineTable,
  layoutRuns,
  lineTable,
  type PreparedRun,
  prepareRun,
} from "./runs.js";
import type { Scheduler } from "./scheduler.js";
import { type Live, live, type SeqState, state, type Unsubscribable, type WritableLive } from "./state.js";
import type { TextManager } from "./text.js";
import type { UiManager } from "./ui.js";
import type { ViewportProps, ViewportState, Viewports } from "./viewport.js";

/**
 * The layout's and engine's sizes and switches (seqspec:viewer#render-props); a preset's `render` over them per
 * section.
 */
export interface RenderProps extends ViewportProps {
  readonly headerWidth: number;
  readonly rowHeight: number;
  readonly rowBorders: boolean;
  readonly blockGap: number;
  readonly stackGap: number;
  readonly stackTitle: number;
  readonly gutterLeft: number;
  readonly gutterRight: number;
  /** Wrapped: the letters' font px; 0 (the default): `maxLetterPx`'s. */
  readonly wrappedFontPx: number;
  /** Wrapped: the lines' numbers' gutter, px, where there are no group labels. */
  readonly wrappedGutter: number;
  /** Wrapped: positions a group (a line holds whole groups). */
  readonly wrappedGroup: number;
  /** Wrapped: px between groups on a line. */
  readonly wrappedGroupGap: number;
  /** Wrapped: px between lines (a block's first line keeps `stackGap` from the last). */
  readonly wrappedGap: number;
  /** Wrapped: px between a line's letters and its first annotation lane (underlines). */
  readonly wrappedTrackGap: number;
  /** Wrapped: an annotation lane's height under the text, px (its underline 2 px less). */
  readonly wrappedLaneHeight: number;
  /**
   * Wrapped: px between the canvas's edges and the text on every side (the skin's `--seq-canvas-edge` default), the
   * line numbers' gutter inside it (seqspec:viewer#text-measuring).
   */
  readonly wrappedInset: number;
  /**
   * Wrapped: each group's number in the section's numbering: its last position's above the line, right over its last
   * letter, in the ruler's band (`aboveEnd`); its first position's as a superscript before it (`inlineStart`); or
   * none.
   */
  readonly wrappedGroupLabels: "aboveEnd" | "inlineStart" | "none";
  readonly labelFontPx: number;
  /**
   * The largest per-position text (a block's letters, another polymer's residues, three-letter codes, atom symbols),
   * CSS px; 0 (the default): the header labels' size (`labelFontPx`), so it follows whatever sizes those
   * (seqspec:viewer#the-render-backend, per-position text).
   */
  readonly maxLetterPx: number;
  readonly controlWidth: number;
  readonly minLetterCell: number;
  readonly maxCell: number;
  readonly denseRowHeight: number;
  readonly maxRowScale: number;
  readonly ligandCells: number;
  /**
   * A dot plot is never flatter than its `aspect_ratio` stretched this much (seqspec:format#layer-types): its height at
   * least the data area's width ÷ (ratio × this), past the section's visible rows if need be (the rows scroll).
   */
  readonly dotPlotMaxStretch: number;
  readonly reductions: boolean;
  readonly blit: boolean;
  readonly verticalOverview: number;
  /** The rows' vertical navigation: a scrollbar, the overview map, both, or none. */
  readonly verticalNav: "scrollbar" | "overview" | "both" | "none";
  /** Row zoom as a control: the overview window's top and bottom handles, the scrollbar thumb's ends. */
  readonly verticalZoom: boolean;
  /**
   * The circular geometry (seqspec:viewer#the-render-backend): px between rings, the hub's radius (of the radius), the
   * label ring.
   */
  readonly ringGap: number;
  readonly hubRadius: number;
  readonly labelRing: boolean;
  /**
   * The data columns' heads: `shared`, on the last pinned entry's line, a row reserved above the first scrolled entry
   * only when nothing is pinned; `ownRow`, always a row of their own there, never a track's line.
   */
  readonly dataColumnHeads: "shared" | "ownRow";
  /** The least data area, CSS px, data columns leave: they yield from the last below it; 0: never. */
  readonly dataFloor: number;
  /**
   * The hovered element outlined in an ink derived from its color (seqspec:interaction#hover; the style key
   * `hover_emphasis`).
   */
  readonly hoverEmphasis: boolean;
}
export const DefaultRenderProps: RenderProps = {
  headerWidth: 150,
  rowHeight: 22,
  rowBorders: false,
  blockGap: 28,
  stackGap: 14,
  stackTitle: 22,
  gutterLeft: 8,
  gutterRight: 8,
  wrappedFontPx: 0,
  wrappedGutter: 56,
  wrappedGroup: 10,
  wrappedGroupGap: 8,
  wrappedGap: 4,
  wrappedTrackGap: 2,
  wrappedLaneHeight: 5,
  wrappedInset: 8,
  wrappedGroupLabels: "aboveEnd",
  labelFontPx: 12,
  maxLetterPx: 0,
  controlWidth: 0,
  minLetterCell: 8,
  maxCell: 28,
  denseRowHeight: 4,
  maxRowScale: 3,
  ligandCells: 3,
  dotPlotMaxStretch: 2.5,
  reductions: true,
  blit: true,
  verticalOverview: 0,
  verticalNav: "scrollbar",
  verticalZoom: false,
  ringGap: 2,
  hubRadius: 0.35,
  labelRing: true,
  dataColumnHeads: "shared",
  dataFloor: 64,
  hoverEmphasis: true,
};

/** The vertical overview's px where a section's inline style asks for the map and its preset gives none (as `msa`'s). */
const OVERVIEW_WIDTH = 72;
/** The vertical overview's width: its px when `verticalNav` shows the map. */
const overviewWidth = (p: RenderProps) =>
  p.verticalNav === "overview" || p.verticalNav === "both" ? p.verticalOverview : 0;
/** Whether a `box` section's rows get a scrollbar when they overflow. */
export const hasScrollbar = (p: RenderProps) => p.verticalNav === "scrollbar" || p.verticalNav === "both";

/**
 * Render props over the defaults, in one shape: a clone of the defaults assigned over (a spread of two records can
 * leave an object V8 reads slowly, and the viewport's clamp reads it per move).
 */
export function mergeProps(...parts: (Partial<RenderProps> | null | undefined)[]): RenderProps {
  const out = { ...DefaultRenderProps };
  for (const p of parts) if (p) Object.assign(out, p);
  return out;
}

/** The per-position text's cap (`maxLetterPx`, else the header labels' size), CSS px. */
export const letterCap = (p: RenderProps): number => (p.maxLetterPx > 0 ? p.maxLetterPx : p.labelFontPx);

/** The render props that are numbers: px, counts, shares. */
export type NumericRenderProp = {
  [K in keyof RenderProps]: RenderProps[K] extends number ? K : never;
}[keyof RenderProps];
/**
 * A host's bounds for its embedding (`SeqViewerConfig.Render.Clamp`), applied after everything else: a phone's `{ min:
 * { dataFloor: 96 } }` holds whatever a preset or a document's style says.
 */
export interface RenderClamp {
  readonly min?: Partial<Record<NumericRenderProp, number>>;
  readonly max?: Partial<Record<NumericRenderProp, number>>;
}

/**
 * A section's render props, lowest first: the defaults, the host's (`spec.render`, `setRenderProps`), its preset's
 * `render`, its inline style's keys (seqspec:format#sections); then the host's clamp bounds the numbers.
 */
export function sectionProps(
  host: RenderProps,
  preset: Partial<RenderProps> | undefined,
  inline: Partial<RenderProps> | null,
  clamp: RenderClamp | null,
): RenderProps {
  const out = mergeProps(host, preset, inline) as { -readonly [K in keyof RenderProps]: RenderProps[K] };
  // A document asking for the map gets one where its preset gives the map no width.
  if (inline?.verticalNav && inline.verticalNav !== "scrollbar" && out.verticalOverview === 0)
    out.verticalOverview = OVERVIEW_WIDTH;
  // Clamped, then each size that follows another derived from it clamped (the letters' cap from the header labels', the
  // text layout's letters from the cap), then clamped again: a clamp of `labelFontPx` reaches the cap it gives, the
  // cap's reaches the text layout's letters, and one of the derived sizes holds (VL-11; seqspec:viewer#render-props:
  // the host's clamp last). A size left to follow (0) is not clamped before it is derived: a `min` would make it a
  // value of its own, a cap (LC-30).
  bound(out, clamp, true);
  out.maxLetterPx = letterCap(out);
  bound(out, clamp, true);
  if (!(out.wrappedFontPx > 0)) out.wrappedFontPx = out.maxLetterPx;
  bound(out, clamp, false);
  return out;
}
/** The sizes 0 leaves to follow another (`letterCap`, `wrappedFontPx`). */
const FOLLOWING: readonly string[] = ["maxLetterPx", "wrappedFontPx"];
function bound(out: { -readonly [K in keyof RenderProps]: RenderProps[K] }, clamp: RenderClamp | null, auto: boolean) {
  const skip = (k: string) => auto && out[k as NumericRenderProp] === 0 && FOLLOWING.includes(k);
  for (const [k, v] of Object.entries(clamp?.min ?? {}))
    if (!skip(k)) out[k as NumericRenderProp] = Math.max(out[k as NumericRenderProp], v);
  for (const [k, v] of Object.entries(clamp?.max ?? {}))
    if (!skip(k)) out[k as NumericRenderProp] = Math.min(out[k as NumericRenderProp], v);
}

/** `labels`: wrapped, the band above a line holding its group labels (`aboveEnd`). */
export const EntryKind = { track: 0, group: 1, title: 2, labels: 3 } as const;
/** A data column's width when its renderer gives none. */
const COLUMN_WIDTH = 80;
/** Below the rows of a section. */
const BOTTOM = 8;
/** An Int32Array twice as long, its values kept (edge rate). */
const grown = (a: Int32Array): Int32Array<ArrayBuffer> => {
  const b = new Int32Array(2 * a.length);
  b.set(a);
  return b;
};
/** The circular geometry: px outside the outermost ring for its labels (`labelRing`), and the least outer radius. */
const RING_LABELS = 28;
const RING_MIN = 24;
/** A circular section's box without a height: its width, within these (CSS px); a ring at least this wide. */
const RING_BOX = [320, 560] as const;
const RING_FLOOR = 3;
/** Px between a ring label and its tick, and past it; the ruler's ticks (renderers/ruler.ts: 4 px at the bottom). */
const RING_LABEL_GAP = 5;
/** Rings with item labels outside (ring-labels.ts): px kept for them beside the ring, and above and below it. */
const RING_ITEM_SIDE = 96;
const RING_ITEM_TOP = 12;
const RULER_TICK = 4;
/** `setRows`' scratch: the rows' height, the row scale that fits them, the px above them (doubles, not boxed). */
const ROWS = new Float64Array(3);

/** The section's `rows` layer: display rows over its lanes, collapsed clades as wedges. */
export interface RowModel {
  layer: number;
  track: number;
  alignment: number;
  /** Its first entry; a lazy layer has one entry for all rows. */
  entry: number;
  lazy: boolean;
  /** Display rows; each `height` px (rows scale with `rowScale`). */
  count: number;
  height: number;
  laneStart: number;
  laneCount: number;
  /** Display row → lane; null: `laneStart + d`. */
  lane: Int32Array | null;
  /** Display row → a collapsed clade's node (a wedge), -1: a row; null: none collapsed. */
  node: Int32Array | null;
  /** Lane − laneStart → display row, -1 inside a collapsed clade; null: identity. */
  display: Int32Array | null;
  /** Bumped on each build of the rows (the tree lays itself out again). */
  version: number;
}

/** Px right of the data area before the vertical overview: the right gutter, wrapped at least the text's inset. */
const rightOf = (p: RenderProps, wrapped: boolean) =>
  wrapped ? Math.max(p.gutterRight, p.wrappedInset) : p.gutterRight;

/** Struct-of-arrays, one record per section, reused across rebuilds (arrays grow only). */
export interface SectionLayout {
  /** `Geometry`. */
  geometry: number;
  width: number;
  height: number;
  dataLeft: number;
  dataRight: number;
  headerWidth: number;
  entryCount: number;
  pinnedCount: number;
  /** `EntryKind`: a track row, a group header, a stacked block's title. */
  entryKind: Uint8Array;
  entryTrack: Int32Array;
  entryGroup: Int32Array;
  entryLayer: Int32Array;
  /** A rows entry's display row; wrapped: the line. */
  entryRow: Int32Array;
  entryBlock: Int32Array;
  /** Pinned entries: surface y; scrolled ones: y in the scrolled content (shown at `scrolledTop + top − scrollY`). */
  entryTop: Float64Array;
  entryHeight: Float32Array;
  /** A lazy run: its row count. */
  entryRun: Int32Array;
  columnCount: number;
  columnX: Float32Array;
  columnWidth: Float32Array;
  /** The data columns' heads band (surface y), above the first scrolled entry; 0 high: none (`dataColumnHeads`). */
  headsTop: number;
  headsHeight: number;
  navHeight: number;
  /**
   * The `rulerRail` navigator (seqspec:viewer#extension-points): the pinned ruler entry whose bottom border is the
   * rail, its window
   * a thicker segment; -1: none (a section without a ruler takes the slim `rail` row, `navHeight`).
   */
  railEntry: number;
  /**
   * Where the rail's grab zone starts (surface y; it ends at the ruler's bottom): the whole ruler row on a blocks
   * axis, its lower 40% on an alignment's columns, where dragging the ruler above it selects columns.
   */
  railTop: number;
  vnavWidth: number;
  /** Vertical navigation shown: the rows overflow, or (with `verticalZoom`) they are zoomed below their scale, so it stays to zoom back. */
  vscroll: boolean;
  pinnedHeight: number;
  /** Where the scrolled rows start: below the navigator and the pinned entries. */
  scrolledTop: number;
  /** The scrolled rows plus the bottom margin: the viewport's vertical extent. */
  scrolledHeight: number;
  /** The whole section's: navigator, pinned, scrolled. */
  contentHeight: number;
  /** `content` (0) or `box` (1). */
  heightMode: number;
  /** The fitted row scale (box), else 1. */
  minRowScale: number;
  rows: RowModel;
  lines: LineTable;
  /**
   * Per block b, the scrolled entries that draw it, [2b, 2b + 1]: stacked, its title and lane (its marks stay there);
   * wrapped, its lines'; side by side, every scrolled entry.
   */
  blockEntries: Int32Array;
  /**
   * The circular geometry, CSS px: the center, the outer radius, and the radius per px of the entries' content y from
   * `scrolledTop` (rings from the outside in); `ringW` the axis unrolled (`2π·ringR`), the x every draw and hit of the
   * section is in before the polar map (`Ring`).
   */
  ringCx: number;
  ringCy: number;
  ringR: number;
  ringScale: number;
  ringW: number;
  /**
   * The rings' outer edge as drawn, in that frame (CSS px from the section's top): an outermost ruler's tick tops (the
   * room above them is its numbers'), else the first ring's top. Marks span from it to the last ring (CG1b).
   */
  ringTop: number;
}

/** Reused. Rows of the section's `rows` layer in view; tops in the scrolled content. */
export interface VisibleRows {
  count: number;
  entry: Int32Array;
  displayRow: Int32Array;
  dataRow: Int32Array;
  top: Float64Array;
  height: Float32Array;
}

const I = (n: number) => new Int32Array(n);
function emptyLayout(): SectionLayout {
  return {
    geometry: AxisGeometry.linear,
    width: 0,
    height: 0,
    dataLeft: 0,
    dataRight: 0,
    headerWidth: 0,
    entryCount: 0,
    pinnedCount: 0,
    entryKind: new Uint8Array(0),
    entryTrack: I(0),
    entryGroup: I(0),
    entryLayer: I(0),
    entryRow: I(0),
    entryBlock: I(0),
    entryTop: new Float64Array(0),
    entryHeight: new Float32Array(0),
    entryRun: I(0),
    columnCount: 0,
    columnX: new Float32Array(0),
    columnWidth: new Float32Array(0),
    headsTop: 0,
    headsHeight: 0,
    navHeight: 0,
    railEntry: -1,
    railTop: 0,
    vnavWidth: 0,
    vscroll: false,
    pinnedHeight: 0,
    scrolledTop: 0,
    scrolledHeight: 0,
    contentHeight: 0,
    heightMode: Height.box,
    minRowScale: 1,
    ringCx: 0,
    ringCy: 0,
    ringR: 0,
    ringScale: 1,
    ringW: 0,
    ringTop: 0,
    rows: {
      layer: -1,
      track: -1,
      alignment: -1,
      entry: -1,
      lazy: false,
      count: 0,
      height: 0,
      laneStart: 0,
      laneCount: 0,
      lane: null,
      node: null,
      display: null,
      version: 0,
    },
    lines: lineTable(),
    blockEntries: I(0),
  };
}

function grow(l: SectionLayout, c: number) {
  const copy = <T extends Int32Array | Float64Array | Float32Array | Uint8Array>(a: T, b: T) => {
    b.set(a as never);
    return b;
  };
  l.entryKind = copy(l.entryKind, new Uint8Array(c));
  l.entryTrack = copy(l.entryTrack, I(c));
  l.entryGroup = copy(l.entryGroup, I(c));
  l.entryLayer = copy(l.entryLayer, I(c));
  l.entryRow = copy(l.entryRow, I(c));
  l.entryBlock = copy(l.entryBlock, I(c));
  l.entryTop = copy(l.entryTop, new Float64Array(c));
  l.entryHeight = copy(l.entryHeight, new Float32Array(c));
  l.entryRun = copy(l.entryRun, I(c));
}

function push(
  l: SectionLayout,
  kind: number,
  track: number,
  group: number,
  layer: number,
  row: number,
  block: number,
  top: number,
  height: number,
  run: number,
): number {
  const n = l.entryCount;
  if (n >= l.entryTrack.length) grow(l, Math.max(16, n * 2));
  l.entryKind[n] = kind;
  l.entryTrack[n] = track;
  l.entryGroup[n] = group;
  l.entryLayer[n] = layer;
  l.entryRow[n] = row;
  l.entryBlock[n] = block;
  l.entryTop[n] = top;
  l.entryHeight[n] = height;
  l.entryRun[n] = run;
  l.entryCount = n + 1;
  return height;
}

/** The scrolled entry at content y (the last starting at or above it); -1: none. Allocation-free. */
function entryAt(l: SectionLayout, y: number): number {
  let lo = l.pinnedCount;
  let hi = l.entryCount - 1;
  if (hi < lo || y < l.entryTop[lo]) return -1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (l.entryTop[mid] <= y) lo = mid;
    else hi = mid - 1;
  }
  return y < l.entryTop[lo] + l.entryHeight[lo] ? lo : -1;
}
/** The scrolled entries over content [range[0], range[1]): the first and the last, written back into `range`. */
function entriesIn(l: SectionLayout, range: Float64Array) {
  const y0 = range[0];
  const y1 = range[1];
  let lo = l.pinnedCount;
  let hi = l.entryCount - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (l.entryTop[mid] <= y0) lo = mid;
    else hi = mid - 1;
  }
  range[0] = lo;
  hi = l.entryCount - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (l.entryTop[mid] < y1) lo = mid;
    else hi = mid - 1;
  }
  range[1] = Math.min(l.entryCount - 1, lo);
}
/** A display row's lane; -1: a wedge. */
function laneOf(r: RowModel, d: number): number {
  if (r.node && r.node[d] >= 0) return -1;
  return r.lane ? r.lane[d] : r.laneStart + d;
}
/** A lane's display row (a row of a collapsed clade: -1). */
function displayOfLane(r: RowModel, lane: number): number {
  const k = lane - r.laneStart;
  if (k < 0 || k >= r.laneCount) return -1;
  return r.display ? r.display[k] : k;
}
/** The entry a display row draws on (lazy: the one entry, at that row). */
function entryOfDisplay(r: RowModel, d: number): number {
  return d < 0 || d >= r.count ? -1 : r.lazy ? r.entry : r.entry + d;
}

/** An expanded split track's row: its value's row, -1 for any other entry. Allocation-free. */
function splitRow(l: SectionLayout, e: number): number {
  const rows = l.rows.layer >= 0 && l.entryTrack[e] === l.rows.track;
  const track = l.geometry === AxisGeometry.linear && l.entryKind[e] === EntryKind.track;
  return track && !rows ? Math.max(-1, l.entryRow[e]) : -1;
}

/**
 * Whether a track entry draws `rowBorders`' line at its bottom, inside its box: an entry over 8 px, or the last row of
 * an expanded split track (thin split rows share the track's). The line is a device px (two at DPR 2), and
 * the entry's content stands above it, so it centers between the lines: the draw, its overlays, the hit and the
 * keyboard cursor take `max(1, round(dpr)) / dpr` off the entry's height. Allocation-free; no double crosses the call.
 */
export function entryBorder(l: SectionLayout, e: number, rowBorders: boolean): boolean {
  if (!rowBorders) return false;
  return l.entryHeight[e] > 8 || lastSplitRow(l, e);
}

/** Whether entry e is the last row of an expanded split track. Allocation-free. */
function lastSplitRow(l: SectionLayout, e: number): boolean {
  const n = e + 1;
  return splitRow(l, e) >= 0 && (n >= l.entryCount || splitRow(l, n) < 0 || l.entryTrack[n] !== l.entryTrack[e]);
}

/**
 * Whether a split row draws a guide at its bottom edge: any row of an expanded split track but the last, whose
 * bottom is the track's. The guide is a hairline in the grid tone under the data, on or off `rowBorders`, and shifts
 * nothing; the composer draws it only where the row is at least `SPLIT_GUIDE_PX` device px tall.
 */
function splitGuide(l: SectionLayout, e: number): boolean {
  return splitRow(l, e) >= 0 && !lastSplitRow(l, e);
}

/** The least height, device px, a split row has for its guide (below it the rows are too thin for lines to help). */
export const SPLIT_GUIDE_PX = 4;

export const SectionLayouts = {
  entryAt,
  entriesIn,
  laneOf,
  displayOfLane,
  entryOfDisplay,
  splitRow,
  entryBorder,
  splitGuide,
};

/** A layer's `facing` on an entry: on an alignment, a pinned track's residues are the rows below it. */
export function facingOf(l: SectionLayout, ax: Axis, layer: number, track: number, entry: number): number {
  // On rings pairs face the hub: their ends along their ring's outer edge, against the ring outside it, and their
  // chords inward across the hub.
  if (l.geometry === AxisGeometry.circular) return 1;
  if (l.rows.layer >= 0 && l.rows.track !== track && entry < l.pinnedCount) return -1;
  return layer < ax.facing.length ? ax.facing[layer] : 0;
}

/** What hits and the layout read of a split layer's prepared state (LY2's `Features`): its values and their rows. */
export interface SplitFeatures {
  readonly split: {
    readonly labels: readonly string[];
    readonly rows: readonly unknown[];
    readonly rowOf: readonly Int32Array[];
    /** Per row, the one unit its data sits on (-1: none); absent: none. */
    readonly unit?: Int32Array;
  } | null;
}

const TEXT_COLUMN = { name: "text", width: COLUMN_WIDTH } as ColumnRenderer;
/** A data column type's renderer, else `text`'s: layout and prepare resolve it alike. */
export function columnOf(registry: SeqViewerRegistry, type: string): ColumnRenderer {
  return registry.get("columns", type) ?? registry.get("columns", "text") ?? TEXT_COLUMN;
}

export interface LayoutManager {
  get(section: string): Live<SectionLayout> | null;
  visibleRows(section: string): Live<VisibleRows> | null;
  contentHeight(section: string): SeqState<number> | null;
  axis(section: string): SectionAxisView | null;
  /** Headless and allocation-free; the section by number (hot records) or ref. */
  hitTest(section: number | string, x: number, y: number, out: Hit): boolean;
  /** `hitTest` at the point already in `out` (`out.x`, `out.y`): no double crosses a call (pointer-rate callers). */
  hitTestAt(section: number, out: Hit): boolean;
  sectionIndex(ref: string): number;
  sectionRef(i: number): string | null;
  /** The section's render props: `ctx.renderProps` under its preset's. */
  props(section: string): RenderProps;
  /** The current view's section's preset (behaviors' `forSection` reads its overrides). */
  presetOf(section: number): PresetEntry | null;
  /** The section's navigator as drawn (its preset's; a `rulerRail` without a ruler: the `rail`). */
  navigatorOf(section: number): NavigatorEntry | null;
}

export interface LayoutHost {
  readonly registry: SeqViewerRegistry;
  readonly config: ConfigManager;
  readonly renderProps: SeqState<RenderProps>;
  readonly text: TextManager;
  readonly ui: UiManager;
  readonly viewports: Viewports;
  readonly scheduler: Scheduler;
  report(code: string, message: string): void;
  /** After a section's layout was rebuilt (the input dispatcher hits the pointer again). */
  rebuilt?(section: number): void;
  /** A section's viewport changed (the input dispatcher hits a still pointer again, at the next frame). */
  moved?(section: number): void;
  /** A layer's regions (the stream manager's): hits read streamed items. */
  stream?: StreamLookup;
}

/** One section's layout state: its record, axis, render props and the states it publishes. */
export class SectionEntry {
  readonly ref: string;
  /** In the current view; -1: not in it. */
  index = -1;
  model: ViewModel | null = null;
  axis: Axis | null = null;
  props: RenderProps = DefaultRenderProps;
  /** Its preset, for behaviors' `forSection`. */
  preset: PresetEntry | null = null;
  nav = 0;
  /** The preset's navigator, and the one drawn (a `rulerRail` section without a ruler takes the `rail`). */
  navigator: NavigatorEntry | null = null;
  drawnNavigator: NavigatorEntry | null = null;
  /** The preset's navigator lives on the ruler (`rulerRail`). */
  rail = false;
  readonly layout = emptyLayout();
  readonly live: WritableLive<SectionLayout>;
  readonly rows: WritableLive<VisibleRows>;
  readonly content = state(0);
  readonly view = new SectionAxisView();
  runs: PreparedRun[] | null = null;
  /** `text.version` when the runs were prepared. */
  runsText = -1;
  /** What the runs' group labels were prepared for: on, the font px, each block's numbering and its path. */
  runsLabels: readonly unknown[] = [];
  dirty = true;
  /** A row zoom waits: the rows' height and the entries below them change, nothing else (`Layouts.rescale`). */
  rescaled = false;
  /** What it was built for: width, height, row scale. */
  readonly built = new Float64Array(3);
  /** The scrolled entries whose height follows the row scale (a renderer's `zooms`), from the last build. */
  zoomed: Int32Array<ArrayBuffer> = new Int32Array(8);
  zoomedCount = 0;
  /** The section's `rows` layer (its first, VW4), and its rows' key when last laid out. */
  rowsLayer = -1;
  readonly rowsKey = {
    layer: -1,
    laneStart: -1,
    laneCount: -1,
    lanes: null as ResolvedLanes | null,
    clades: null as unknown,
  };
  /** Wrapped: tracks × lines band heights, and one layer's per line (scratch, grown only). */
  bands = new Float32Array(0);
  line = new Float32Array(0);
  /** The visible display rows last published. */
  readonly shownRows = new Int32Array(2);
  /** From the last build, for the scroll: the rows' top in the content, their height (0: none), count, the rows region's top. */
  readonly rowSpan = new Float64Array(4);
  readonly range = new Int32Array(2);
  sub: Unsubscribable | null = null;
  /** Bound by the manager: builds this section's layout if stale. */
  ensure: (() => void) | null = null;
  /** The axis's hook: the scrolled entries over content [y0, y1). */
  readonly entriesIn = (range: Float64Array) => entriesIn(this.layout, range);

  constructor(ref: string) {
    this.ref = ref;
    this.live = live(this.layout);
    this.rows = live<VisibleRows>({
      count: 0,
      entry: I(0),
      displayRow: I(0),
      dataRow: I(0),
      top: new Float64Array(0),
      height: new Float32Array(0),
    });
  }
}

export class Layouts implements LayoutManager {
  /** For tests: layouts built, runs prepared, `layoutRuns` calls. */
  readonly stats = { builds: 0, rescales: 0, prepareRuns: 0, layoutRuns: 0 };
  private readonly map = new Map<string, SectionEntry>();
  /** Sections that left the shown document, closed at the end of the update unless they return in it. */
  private readonly leaving = new Map<string, SectionEntry>();
  private current: SectionEntry[] = [];
  private model: ViewModel | null = null;
  private frame = -1;
  private building = false;
  private readonly reported = new Set<string>();
  private readonly host: LayoutHost;

  constructor(host: LayoutHost) {
    this.host = host;
    const dirty = () => this.dirtyAll();
    host.ui.expanded.subscribe(dirty);
    host.ui.clades.subscribe(dirty);
    // Props, presets and fonts change the axes too (gaps, the data area, the wrapped cell).
    const again = () => this.model && this.update(this.model);
    for (const s of [host.ui.style, host.renderProps, host.text.version, host.registry.version] as SeqState<unknown>[])
      s.subscribe(again);
  }

  sectionIndex(ref: string): number {
    return this.map.get(ref)?.index ?? -1;
  }
  sectionRef(i: number): string | null {
    return this.current[i]?.ref ?? null;
  }
  props(section: string): RenderProps {
    return this.map.get(section)?.props ?? this.host.renderProps.value;
  }
  presetOf(section: number): PresetEntry | null {
    return this.current[section]?.preset ?? null;
  }
  navigatorOf(section: number): NavigatorEntry | null {
    return this.current[section]?.drawnNavigator ?? null;
  }

  get(section: string): Live<SectionLayout> | null {
    const e = this.ready(section);
    return e ? e.live : null;
  }
  visibleRows(section: string): Live<VisibleRows> | null {
    const e = this.ready(section);
    return e ? e.rows : null;
  }
  contentHeight(section: string): SeqState<number> | null {
    const e = this.ready(section);
    return e ? e.content : null;
  }
  axis(section: string): SectionAxisView | null {
    const e = this.ready(section);
    return e?.axis ? this.synced(e) : null;
  }
  hitTest(section: number | string, x: number, y: number, out: Hit): boolean {
    out.x = x;
    out.y = y;
    return this.hitEntry(typeof section === "number" ? this.current[section] : this.map.get(section), out);
  }
  hitTestAt(section: number, out: Hit): boolean {
    return this.hitEntry(this.current[section], out);
  }
  /** A section with no layout hits nothing, and says so: `out` is reset, not left at the last element (IN2). */
  private hitEntry(e: SectionEntry | undefined, out: Hit): boolean {
    if (e) this.ensure(e);
    if (!e?.axis || !e.model) {
      Hit.reset(out);
      out.section = e ? e.index : -1;
      return false;
    }
    const residue = this.host.config.get(SeqViewerConfig.UI.OwnRowHover) === "residue";
    return hitTest(e, this.synced(e), this.host.registry, out, this.host.stream ?? null, residue);
  }

  private ready(section: string): SectionEntry | null {
    const e = this.map.get(section) ?? null;
    if (e) this.ensure(e);
    return e;
  }
  /** The axis view under the current layout and viewport; re-synced only when either version moved. */
  private synced(e: SectionEntry): SectionAxisView {
    const vp = this.host.viewports.get(e.ref);
    const v = e.view;
    if (!vp || !e.axis) return v;
    const lv = e.live.version.value;
    const vv = vp.live.version.value;
    if (v.axis !== e.axis || v.at[0] !== lv || v.at[1] !== vv) {
      v.minLetterCell = e.props.minLetterCell;
      v.sync(e.axis, e.layout, vp.live.current);
      v.runs = e.runs ?? NO_RUNS;
      v.gap = e.props.wrappedGroupGap;
      v.fontPx = e.props.wrappedFontPx;
      v.laneHeight = e.props.wrappedLaneHeight;
      v.trackGap = e.props.wrappedTrackGap;
      v.gutterLeft = e.props.gutterLeft;
      v.gutterRight = e.props.gutterRight;
      v.reductions = e.props.reductions;
      v.at[0] = lv;
      v.at[1] = vv;
    }
    return v;
  }

  /**
   * Lifecycle: one entry per section of every view of the shown document, beside its viewport. A ref that leaves and
   * returns within one update keeps its entry, as `Viewports.sync` keeps its viewport.
   */
  sync(sections: readonly string[]) {
    for (const ref of sections) {
      if (this.map.has(ref)) continue;
      const back = this.leaving.get(ref);
      this.leaving.delete(ref);
      if (back) {
        this.map.set(ref, back);
        continue;
      }
      const e = new SectionEntry(ref);
      e.ensure = () => this.ensure(e);
      this.map.set(ref, e);
      const s = this.host.viewports.get(ref)?.live;
      if (s) e.sub = s.version.subscribe(() => this.followViewport(e, s.current));
    }
    for (const [ref, e] of this.map) {
      if (sections.includes(ref)) continue;
      this.map.delete(ref);
      if (this.leaving.size === 0) void Promise.resolve().then(this.drop);
      this.leaving.set(ref, e);
    }
  }
  private readonly drop = () => {
    for (const e of this.leaving.values()) this.close(e);
    this.leaving.clear();
  };

  /** The current view on every graph emission: render props, presets and axes per section, then the viewports. */
  update(model: ViewModel | null) {
    this.model = model;
    const v = model?.resolved ?? null;
    const set = model?.set ?? null;
    for (const e of this.map.values()) e.index = -1;
    this.current = [];
    for (let i = 0; v && set && i < v.sections.count; i++) {
      const ref = set.strings.get(v.sections.ref[i]);
      const e = this.map.get(ref);
      if (!e) continue;
      e.index = i;
      this.current[i] = e;
      this.configure(e, model as ViewModel, i);
    }
    for (const e of this.map.values()) {
      if (e.index >= 0) continue;
      e.axis = null;
      e.model = null;
    }
    this.host.viewports.update(v, set, (i) => this.current[i]?.axis ?? null);
    this.dirtyAll();
  }

  private configure(e: SectionEntry, model: ViewModel, i: number) {
    const { resolved: v, set } = model;
    const preset = this.preset(e.ref, v, set, i);
    e.preset = preset;
    const clamp = this.host.config.get(SeqViewerConfig.Render.Clamp);
    const props = sectionProps(this.host.renderProps.value, preset?.render, v.sections.render[i], clamp);
    const wrapped = preset?.geometry ? preset.geometry === "wrapped" : v.layout === Layout.text;
    const circular = preset?.geometry === "circular" && !wrapped;
    e.rowsLayer = wrapped ? -1 : rowsLayer(v, set, i);
    const nav = navigatorOf(this.host.registry, preset?.navigator);
    e.navigator = e.drawnNavigator = nav;
    e.rail = nav?.parts.onRuler === true;
    e.nav = nav ? nav.height(props) : 0;
    const old = e.axis;
    const group = e.props.wrappedGroup;
    e.model = model;
    e.props = props;
    const vp = this.host.viewports.get(e.ref);
    if (vp) vp.props = props;
    const ax = new Axis(v, set, i, props.blockGap, props.ligandCells);
    const l = e.layout;
    l.geometry = wrapped ? AxisGeometry.wrapped : circular ? AxisGeometry.circular : AxisGeometry.linear;
    l.headerWidth = props.headerWidth;
    l.navHeight = e.nav;
    l.vnavWidth = overviewWidth(props);
    // Data columns between the header and the data area (the section's `rows` layer's), never in wrapped.
    const rows = e.rowsLayer;
    const columns = (rows >= 0 ? (v.layers.params[rows].data_columns as { type?: string }[] | null) : null) ?? [];
    if (l.columnX.length < columns.length) {
      l.columnX = new Float32Array(columns.length);
      l.columnWidth = new Float32Array(columns.length);
    }
    // The data area's left with the first k columns shown (they yield to a narrow width: `Axis.fitWidth`).
    const lefts = new Float64Array(columns.length + 1);
    let x = props.headerWidth + props.controlWidth;
    lefts[0] = x + props.gutterLeft;
    columns.forEach((c, k) => {
      l.columnX[k] = x;
      x += l.columnWidth[k] = this.column(e, String(c.type ?? "text")).width;
      lefts[k + 1] = x + props.gutterLeft;
    });
    // Wrapped, the gutter holds the lines' numbers, drawn only without group labels: else no width.
    const gutter = props.wrappedGroupLabels === "none" ? props.wrappedGutter : 0;
    if (wrapped) lefts[0] = props.headerWidth + props.wrappedInset + gutter;
    ax.lefts = lefts;
    ax.reserve = rightOf(props, wrapped) + overviewWidth(props);
    ax.floor = props.dataFloor;
    ax.left = lefts[columns.length];
    ax.insetX = ax.left + ax.reserve;
    ax.columnsShown = columns.length;
    const vp0 = this.host.viewports.get(e.ref)?.live.current;
    if (vp0 && vp0.width > 0) ax.fitWidth(vp0.width);
    l.columnCount = ax.columnsShown;
    if (wrapped)
      ax.fixedCell = this.host.text.cell(this.host.text.font("mono", props.wrappedFontPx), alphabetOf(ax, 0));
    ax.entriesIn = e.entriesIn;
    ax.rows = l.rows;
    ax.ensure = e.ensure;
    e.axis = ax;
    // Runs survive an update that keeps the axis, the cell, the groups and the fonts (a color, a track moved).
    const same = old?.equals(ax) && old.fixedCell === ax.fixedCell && group === props.wrappedGroup;
    // ...and the group labels' numbering (the Numbering switch re-resolves the view, keeping the axis).
    const b0 = v.sections.blockStart[i];
    const labels: unknown[] = [props.wrappedGroupLabels, props.wrappedFontPx];
    for (let b = 0; wrapped && b < ax.count; b++)
      labels.push(v.blocks.numbering[b0 + b], v.blocks.numberingPath[b0 + b]);
    if (!same || e.runsText !== this.host.text.version.value || !sameList(labels, e.runsLabels)) e.runs = null;
    e.runsLabels = labels;
  }

  /** A data column's renderer; an unknown type is drawn as `text`, reported once per section and type. */
  private column(e: SectionEntry, type: string): ColumnRenderer {
    const r = columnOf(this.host.registry, type);
    if (r.name !== type && !this.reported.has(`${e.ref}|column ${type}`)) {
      this.reported.add(`${e.ref}|column ${type}`);
      this.host.report("svs.render.fallback", `section "${e.ref}": no data column type "${type}"; drawn as "text"`);
    }
    return r;
  }

  /** `ui.style`, else the section's `style` through `Presets.Alias`, else `Presets.Default` by layout and axis. */
  private preset(ref: string, v: ResolvedView, set: SeqSet, i: number) {
    const c = this.host.config;
    const style = v.sections.style[i] >= 0 ? set.strings.get(v.sections.style[i]) : null;
    const defaults = c.get(SeqViewerConfig.Presets.Default);
    const columns = v.blocks.alignment[v.sections.blockStart[i]] >= 0;
    const fallback = v.layout === Layout.text ? defaults.text : columns ? defaults.columns : defaults.blocks;
    const key =
      this.host.ui.style.value.get(ref) ?? (style ? (c.get(SeqViewerConfig.Presets.Alias)[style] ?? style) : fallback);
    const found = this.host.registry.get("presets", key);
    if (!found && key !== fallback && !this.reported.has(`${ref}|${key}`)) {
      this.reported.add(`${ref}|${key}`);
      this.host.report("svs.render.preset", `section "${ref}": no preset "${key}"; drawn with "${fallback}"`);
    }
    return found ?? this.host.registry.get("presets", fallback) ?? null;
  }

  private dirtyAll() {
    for (const e of this.map.values()) e.dirty = true;
    if (this.map.size > 0) this.schedule();
  }
  /** A frame reads the layouts if nobody has. No closure here: V8 would allocate its context per call. */
  private schedule() {
    if (this.frame < 0) this.frame = this.host.scheduler.requestFrame(this.onFrame);
  }
  private readonly onFrame = () => {
    this.frame = -1;
    for (const e of this.map.values()) this.ensure(e);
  };

  /**
   * A viewport change: a new size rebuilds; a row scale alone rescales the rows (`rescale`), else rebuilds; a scroll
   * moves the visible rows only.
   */
  private followViewport(e: SectionEntry, s: ViewportState) {
    if (e.index >= 0) this.host.moved?.(e.index);
    if (s.width !== e.built[0] || s.height !== e.built[1]) {
      e.dirty = true;
      this.schedule();
    } else if (s.rowScale !== e.built[2]) {
      e.rescaled = true;
      this.schedule();
    } else if (!e.dirty && !e.rescaled) this.followScroll(e, s, false);
  }

  private ensure(e: SectionEntry) {
    if (this.building) return;
    if (e.dirty) this.rebuild(e);
    else if (e.rescaled) this.rescale(e);
  }
  /**
   * A row zoom without a rebuild (allocation-free): the linear, unstacked rows' height to the row scale, and the
   * entries whose height follows it (`zoomed`: a renderer's `zooms`, a dense matrix's) to theirs, the scrolled entries
   * laid out again from the first, the content's height; a change that shows or hides the vertical overview (the data
   * area's width) rebuilds. Every other entry's height ignores the row scale.
   */
  private rescale(e: SectionEntry) {
    e.rescaled = false;
    const l = e.layout;
    const r = l.rows;
    const vp = this.host.viewports.get(e.ref);
    const s = vp?.live.current;
    const linear = l.geometry === AxisGeometry.linear && e.axis !== null && !e.axis.stacked;
    const rows = r.layer >= 0 && r.count > 0;
    if (!s || !linear || (!rows && e.zoomedCount === 0) || (rows && (r.lazy ? 1 : r.count) + r.entry > l.entryCount)) {
      this.rebuild(e);
      return;
    }
    const p = e.props;
    const was = e.built[2];
    e.built[2] = s.rowScale;
    const height = p.rowHeight * s.rowScale;
    // The new heights' sum against the old, to see the overview before anything moves.
    let delta = rows ? r.count * (height - r.height) : 0;
    for (let z = 0; z < e.zoomedCount; z++) {
      const k = e.zoomed[z];
      delta += this.trackHeight(e, l.entryTrack[k]) - l.entryHeight[k];
    }
    const scrolled = l.scrolledHeight + delta;
    const over = l.heightMode === Height.box && scrolled > s.height - l.scrolledTop + 0.5;
    const vscroll = over || (l.heightMode === Height.box && p.verticalZoom && s.rowScale < 1 - 1e-6);
    if (vscroll !== l.vscroll) {
      e.built[2] = was;
      this.rebuild(e);
      return;
    }
    this.building = true;
    try {
      if (rows) {
        r.height = height;
        const last = r.lazy ? r.entry + 1 : r.entry + r.count;
        for (let k = r.entry; k < last; k++) l.entryHeight[k] = r.lazy ? r.count * height : height;
      }
      for (let z = 0; z < e.zoomedCount; z++)
        l.entryHeight[e.zoomed[z]] = this.trackHeight(e, l.entryTrack[e.zoomed[z]]);
      // As the build lays them out: each scrolled entry's top the sum of the heights above it.
      // The rows' heights in doubles, as the build adds them (the entries keep floats).
      const last = rows ? (r.lazy ? r.entry + 1 : r.entry + r.count) : -1;
      let y = l.pinnedCount < l.entryCount ? l.entryTop[l.pinnedCount] : 0;
      for (let k = l.pinnedCount; k < l.entryCount; k++) {
        l.entryTop[k] = y;
        y += k >= r.entry && k < last ? (r.lazy ? r.count * height : height) : l.entryHeight[k];
      }
      l.scrolledHeight = scrolled;
      l.contentHeight = l.scrolledTop + scrolled;
      this.stats.rescales++;
      // What `rebuild` publishes: a box with `height_rows` keeps its preferred height (its rows' at scale 1); any other
      // section's is its content, in whole px (an int, so setting it boxes nothing), set only when it changed.
      const fixed = e.model && e.index >= 0 ? e.model.resolved.sections.heightRows[e.index] : -1;
      const content = Math.ceil(l.contentHeight) | 0;
      if ((l.heightMode === Height.content || fixed <= 0) && e.content.value !== content) e.content.set(content);
      this.published(e, vp, s);
    } finally {
      this.building = false;
    }
    if (e.index >= 0) this.host.rebuilt?.(e.index);
  }
  private rebuild(e: SectionEntry) {
    this.building = true;
    try {
      e.dirty = false;
      e.rescaled = false;
      const vp = this.host.viewports.get(e.ref);
      const s = vp?.live.current;
      e.built[0] = s?.width ?? 0;
      e.built[1] = s?.height ?? 0;
      e.built[2] = s?.rowScale ?? 1;
      const inset = e.axis?.insetX ?? 0;
      this.build(e, e.built[0], e.built[1], e.built[2]);
      // The overview came or went: the data area's width changed, so the viewport's limits did.
      if (s && e.axis && e.axis.insetX !== inset) vp?.setSize(s.width, s.height, s.dpr);
      this.stats.builds++;
      // Whole px: a small integer is no allocation when the row zoom sets it (`rescale`).
      e.content.set(Math.ceil(e.layout.heightMode === Height.content ? e.layout.contentHeight : this.preferred(e)) | 0);
      this.published(e, vp, s);
    } finally {
      this.building = false;
    }
    if (e.index >= 0) this.host.rebuilt?.(e.index);
  }
  /** A built or rescaled layout given out: the viewport's rows, the live layout, the visible rows. */
  private published(e: SectionEntry, vp: ReturnType<Viewports["get"]>, s: ViewportState | undefined) {
    const l = e.layout;
    if (vp) {
      ROWS[0] = l.scrolledHeight;
      ROWS[1] = l.minRowScale;
      ROWS[2] = l.scrolledTop;
      vp.setRows(ROWS);
    }
    e.live.bump();
    const r = l.rows;
    const on = r.layer >= 0 && r.height > 0 && r.count > 0;
    e.rowSpan[0] = on ? l.entryTop[r.entry] : 0;
    e.rowSpan[1] = on ? r.height : 0;
    e.rowSpan[2] = on ? r.count : 0;
    e.rowSpan[3] = l.scrolledTop;
    if (s) this.followScroll(e, s, true);
  }

  /** `box`: `height_rows` × row height as the preferred height, else as tall as its rows. */
  private preferred(e: SectionEntry): number {
    const l = e.layout;
    const rows = e.model && e.index >= 0 ? e.model.resolved.sections.heightRows[e.index] : -1;
    return rows > 0 ? l.navHeight + l.pinnedHeight + rows * e.props.rowHeight : l.contentHeight;
  }

  private build(e: SectionEntry, width: number, height: number, rowScale: number) {
    const l = e.layout;
    l.entryCount = 0;
    l.pinnedCount = 0;
    l.width = width;
    l.height = height;
    l.rows.layer = -1;
    l.lines.count = 0;
    e.zoomedCount = 0;
    const ax = e.axis;
    const model = e.model;
    if (!ax || !model) {
      l.pinnedHeight = l.scrolledHeight = l.contentHeight = 0;
      l.scrolledTop = l.navHeight = 0;
      return;
    }
    const v = model.resolved;
    const p = e.props;
    ax.fitWidth(width);
    l.columnCount = ax.columnsShown;
    l.dataLeft = ax.left;
    const wrapped = l.geometry === AxisGeometry.wrapped;
    l.dataRight = Math.max(ax.left, width - rightOf(p, wrapped) - overviewWidth(p));
    l.heightMode = v.sections.height[e.index];
    const items = this.items(e, v, model.set);
    // The rulerRail lives on the first pinned ruler; without one, the slim rail row.
    const linearSide = l.geometry === AxisGeometry.linear && !ax.stacked;
    const ruler =
      e.rail && linearSide
        ? items.find((t) => t >= 0 && isPinned(v, model.set, t) && isRuler(v, model.set, t))
        : undefined;
    const fallback = e.rail && ruler === undefined ? (this.host.registry.get("navigators", "rail") ?? null) : null;
    e.drawnNavigator = fallback ?? e.navigator;
    l.navHeight = e.rail ? (fallback?.height(p) ?? 0) : e.nav;
    l.railEntry = -1;
    let pinned = 0;
    let y = l.navHeight;
    if (linearSide)
      for (let k = 0; k < items.length; k++) {
        const t = items[k];
        if (t < 0 || !isPinned(v, model.set, t)) continue;
        if (t === ruler) {
          l.railEntry = l.entryCount;
          l.railTop = y + (ax.alignment[0] >= 0 ? 0.6 * this.trackHeight(e, t) : 0);
        }
        y += push(l, EntryKind.track, t, v.tracks.group[t], firstLayer(v, t), -1, -1, y, this.trackHeight(e, t), 0);
        items[k] = NONE;
        pinned++;
      }
    l.pinnedCount = ax.pinnedCount = pinned;
    // Data columns' heads are a band directly above the first scrolled entry: the last pinned entry's line when they
    // share it, else a row of their own.
    l.headsTop = l.headsHeight = 0;
    if (l.columnCount > 0 && pinned > 0 && p.dataColumnHeads === "shared") {
      l.headsTop = l.entryTop[pinned - 1];
      l.headsHeight = l.entryHeight[pinned - 1];
    } else if (l.columnCount > 0) {
      l.headsTop = y;
      l.headsHeight = p.rowHeight;
      y += p.rowHeight;
    }
    l.pinnedHeight = y - l.navHeight;
    l.scrolledTop = y;
    let h = 0;
    if (l.geometry === AxisGeometry.wrapped) h = this.wrapped(e, items);
    else if (ax.stacked)
      for (let b = 0; b < ax.count; b++) {
        h += b > 0 ? p.stackGap : 0;
        h += push(l, EntryKind.title, -1, -1, -1, -1, b, h, p.stackTitle, 0);
        h = this.scrolled(e, items, b, h, rowScale);
      }
    else h = this.scrolled(e, items, -1, 0, rowScale);
    if (l.geometry === AxisGeometry.circular) {
      this.ring(e, h, width, height);
      this.blockEntries(l, ax);
      return;
    }
    // `content` sections are their rows plus a margin; a `box` section's rows scroll exactly to their last; wrapped
    // text has its inset below.
    const margin = l.heightMode === Height.content && !wrapped ? BOTTOM : 0;
    l.scrolledHeight = h + margin;
    l.contentHeight = l.scrolledTop + l.scrolledHeight;
    const r = l.rows;
    const scaled = r.layer >= 0 ? r.count * r.height : 0;
    const fit = scaled > 0 ? (height - l.scrolledTop - margin - (h - scaled)) / (r.count * p.rowHeight) : 1;
    l.minRowScale = l.heightMode === Height.content ? 1 : Math.max(1e-4, Math.min(1, fit));
    // The vertical overview only when the rows scroll: never an empty box beside rows that all fit.
    const over = l.heightMode === Height.box && l.scrolledHeight > height - l.scrolledTop + 0.5;
    l.vscroll = over || (l.heightMode === Height.box && p.verticalZoom && rowScale < 1 - 1e-6);
    l.vnavWidth = l.vscroll ? overviewWidth(p) : 0;
    l.dataRight = Math.max(ax.left, width - rightOf(p, wrapped) - l.vnavWidth);
    ax.insetX = ax.left + rightOf(p, wrapped) + l.vnavWidth;
    this.blockEntries(l, ax);
  }

  /**
   * Rings: the entries from the outside in, `ringGap` apart, filling the ring from the label ring (`labelRing`) down to
   * the hub (`hubRadius` of the radius): each as wide as its linear height scaled to the radius, at least `RING_FLOOR`
   * px. The radius fills the box, the smaller of its width and height less the label ring; the box is the section's
   * height (`box`), else `height_rows` rows, else as tall as it is wide within [320, 560] px. Nothing scrolls.
   */
  private ring(e: SectionEntry, h: number, width: number, height: number) {
    const l = e.layout;
    const p = e.props;
    let thin = Number.POSITIVE_INFINITY;
    for (let k = l.pinnedCount; k < l.entryCount; k++) {
      l.entryTop[k] += (k - l.pinnedCount) * p.ringGap;
      thin = Math.min(thin, l.entryHeight[k]);
    }
    const bands = h + Math.max(0, l.entryCount - l.pinnedCount - 1) * p.ringGap;
    const model = e.model as ViewModel;
    const rows = e.index >= 0 ? model.resolved.sections.heightRows[e.index] : -1;
    // Item labels outside the ring (a renderer's `ringLabels`, with `labels`): room beside it (a sixth of the width,
    // within [64, 96] px) and above and below it; the section no taller than the ring and that room need.
    const items = p.labelRing && this.ringLabeled(l, model);
    const side = items ? Math.min(RING_ITEM_SIDE, Math.max(64, width / 6)) : 0;
    const square = Math.min(RING_BOX[1], Math.max(RING_BOX[0], width));
    const want =
      rows > 0 ? rows * p.rowHeight : items ? Math.min(square, width - 2 * side + 2 * RING_ITEM_TOP) : square;
    const box = l.heightMode === Height.box && height > l.scrolledTop ? height - l.scrolledTop : want;
    // An outermost ruler's numbers sit just past its ticks (RingBackend): the margin is the widest number less the room
    // above the ticks, found in two passes (that room scales with the radius). Else the label ring's allowance.
    const first = l.pinnedCount < l.entryCount ? l.entryTrack[l.pinnedCount] : -1;
    const ruler = first >= 0 && isRuler(model.resolved, model.set, first);
    const t = this.host.text;
    const widest = ruler ? Math.ceil(t.numberWidth(t.font("ui", 11), Math.max(1, e.axis?.cells ?? 1))) : 0;
    let margin = !p.labelRing ? BOTTOM : ruler ? widest + 2 * RING_LABEL_GAP : RING_LABELS;
    const across = items ? Math.min(width - 2 * side, box - 2 * RING_ITEM_TOP) : Math.min(width, box);
    for (let pass = 0; pass < (ruler && p.labelRing ? 2 : 1); pass++) {
      l.ringR = Math.max(RING_MIN, (across - 2 * margin) / 2);
      const fill = (l.ringR * (1 - p.hubRadius)) / Math.max(1, bands);
      l.ringScale = Math.min(l.ringR / Math.max(1, bands), Math.max(fill, thin > 0 ? RING_FLOOR / thin : fill));
      const above = (l.entryHeight[l.pinnedCount] - RULER_TICK) * l.ringScale;
      if (ruler && p.labelRing) margin = Math.max(BOTTOM, widest + 2 * RING_LABEL_GAP - above);
    }
    l.ringCx = width / 2;
    l.ringCy = l.scrolledTop + Math.max(box, 2 * (l.ringR + margin)) / 2;
    l.ringW = 2 * Math.PI * l.ringR;
    const k0 = l.pinnedCount;
    l.ringTop = l.scrolledTop + (k0 < l.entryCount ? l.entryTop[k0] + (ruler ? l.entryHeight[k0] - RULER_TICK : 0) : 0);
    // Whole px: the composer passes it to its clip (a fractional double crossing a call is boxed).
    l.dataRight = Math.ceil(l.dataLeft + l.ringW);
    l.contentHeight = l.scrolledTop + Math.max(want, 2 * (l.ringR + margin));
    l.scrolledHeight = l.contentHeight - l.scrolledTop;
    l.minRowScale = 1;
    l.vscroll = false;
    l.vnavWidth = 0;
  }

  /** Whether a layer on the rings lists item labels (its renderer's `ringLabels`, its `labels` on). */
  private ringLabeled(l: SectionLayout, model: ViewModel): boolean {
    const rv = model.resolved;
    for (let e = 0; e < l.entryCount; e++) {
      const t = l.entryTrack[e];
      for (let k = l.entryLayer[e]; t >= 0 && k >= 0 && k < rv.layers.count && rv.layers.track[k] === t; k++)
        if (
          rv.layers.params[k].labels === true &&
          this.host.registry.get("renderers", model.prepared[k]?.renderer ?? "")?.ringLabels
        )
          return true;
    }
    return false;
  }

  /** Per block, its scrolled entries (`SectionLayout.blockEntries`): its own when stacked or wrapped. */
  private blockEntries(l: SectionLayout, ax: Axis) {
    if (l.blockEntries.length < 2 * ax.count) l.blockEntries = I(2 * ax.count);
    const own = ax.stacked || l.geometry === AxisGeometry.wrapped;
    for (let b = 0; b < ax.count; b++) {
      l.blockEntries[2 * b] = own ? l.entryCount : l.pinnedCount;
      l.blockEntries[2 * b + 1] = own ? l.pinnedCount : l.entryCount;
    }
    for (let e = l.pinnedCount; own && e < l.entryCount; e++) {
      const b = l.entryBlock[e];
      if (b < 0 || b >= ax.count) continue;
      l.blockEntries[2 * b] = Math.min(l.blockEntries[2 * b], e);
      l.blockEntries[2 * b + 1] = Math.max(l.blockEntries[2 * b + 1], e + 1);
    }
  }

  /**
   * Visible groups and tracks in order: a track as its index, a group header as −2 − g; collapsed groups' children left
   * out but their `summary` tracks, which show only then (seqspec:format#placement).
   */
  private items(e: SectionEntry, v: ResolvedView, set: SeqSet): number[] {
    const out: number[] = [];
    const tr = v.tracks;
    const g = v.groups;
    const emitted = new Set<number>();
    const expanded = this.host.ui.expanded.value;
    const open = (k: number): boolean => {
      if (k < 0) return true;
      if (!open(g.parent[k]) || g.hidden[k]) return false;
      if (!emitted.has(k)) {
        emitted.add(k);
        out.push(-2 - k);
      }
      return expanded.has(set.strings.get(g.ref[k]));
    };
    for (let t = 0; t < tr.count; t++) {
      if (tr.section[t] !== e.index || tr.hidden[t]) continue;
      const k = tr.group[t];
      if (tr.params[t]?.summary !== true) {
        if (open(k)) out.push(t);
      } else if (k >= 0 && open(g.parent[k]) && !g.hidden[k] && !open(k)) out.push(t);
    }
    return out;
  }

  /** Scrolled entries of the linear geometry (one block when stacked), from content y; returns the y after them. */
  private scrolled(e: SectionEntry, items: number[], block: number, y: number, rowScale: number): number {
    const l = e.layout;
    const model = e.model as ViewModel;
    const v = model.resolved;
    const p = e.props;
    for (const t of items) {
      if (t === NONE) continue;
      if (t < 0) {
        y += push(l, EntryKind.group, -1, -2 - t, -1, -1, block, y, p.rowHeight, 0);
        continue;
      }
      // Lanes come from the section's first `rows` layer only (VW4).
      const rows = e.rowsLayer >= 0 && v.layers.track[e.rowsLayer] === t ? e.rowsLayer : -1;
      // An expanded split track: its label's entry (row -2, its other layers), then a row per value.
      const k = rows < 0 ? v.tracks.split[t] : -1;
      const split = k >= 0 && this.host.ui.expanded.value.has(model.set.strings.get(v.tracks.ref[t]));
      const values = split ? ((model.prepared[k]?.state as SplitFeatures | undefined)?.split?.labels.length ?? 0) : 0;
      if (values > 0) {
        y += push(l, EntryKind.track, t, v.tracks.group[t], firstLayer(v, t), -2, block, y, p.rowHeight, 0);
        const h = this.renderer(v, model.set, k)?.splitHeight ?? this.trackHeight(e, t);
        for (let r = 0; r < values; r++) y += push(l, EntryKind.track, t, v.tracks.group[t], k, r, block, y, h, 0);
        continue;
      }
      if (rows < 0 || block >= 0) {
        if (block < 0 && this.zooms(e, t)) {
          if (e.zoomedCount === e.zoomed.length) e.zoomed = grown(e.zoomed);
          e.zoomed[e.zoomedCount++] = l.entryCount;
        }
        y += push(l, EntryKind.track, t, v.tracks.group[t], firstLayer(v, t), -1, block, y, this.trackHeight(e, t), 0);
        continue;
      }
      const r = this.rowModel(e, t, rows, p.rowHeight * rowScale);
      r.entry = l.entryCount;
      if (r.lazy) y += push(l, EntryKind.track, t, v.tracks.group[t], rows, 0, 0, y, r.count * r.height, r.count);
      else
        for (let d = 0; d < r.count; d++)
          y += push(l, EntryKind.track, t, v.tracks.group[t], rows, d, 0, y, r.height, 0);
    }
    return y;
  }

  /**
   * A track's rows over its lanes (tree order and `row_subset` are the format's), collapsed clades folded. A build with
   * the same rows (layer, lanes, clades) keeps them and their version, so a resize or a row zoom lays no tree out again
   * (the key is the entry's, captured when they were made, never the record a build resets).
   */
  private rowModel(e: SectionEntry, t: number, layer: number, height: number): RowModel {
    const model = e.model as ViewModel;
    const { resolved: v, set } = model;
    const r = e.layout.rows;
    const a = (e.axis as Axis).alignment[0];
    const aln = set.alignments[a];
    Object.assign(r, { layer, track: t, alignment: a, height });
    r.laneStart = v.tracks.laneStart[t];
    r.laneCount = v.tracks.laneCount[t];
    const renderer = this.renderer(v, set, layer);
    r.lazy = Boolean(aln?.lazy || renderer?.lazy);
    const clades = (aln && this.host.ui.clades.value.get(set.strings.get(aln.key))) ?? null;
    const key = e.rowsKey;
    const same = key.layer === layer && key.laneStart === r.laneStart && key.laneCount === r.laneCount;
    if (same && key.lanes === v.lanes && key.clades === clades) return r;
    Object.assign(key, { layer, laneStart: r.laneStart, laneCount: r.laneCount, lanes: v.lanes, clades });
    r.version++;
    Object.assign(r, { lane: null, node: null, display: null, count: r.laneCount });
    const tree = aln?.tree;
    if (!clades?.length || !tree || !followsTree(v.layers.params[layer].order, true)) return r;
    const leafOf = new Int32Array(aln.rowCount).fill(-1);
    for (let n = 0; n < tree.count; n++) if (tree.leafRow[n] >= 0) leafOf[tree.leafRow[n]] = n;
    const lane: number[] = [];
    const node: number[] = [];
    const display = new Int32Array(r.laneCount).fill(-1);
    let open = -1;
    for (let k = 0; k < r.laneCount; k++) {
      const leaf = leafOf[Lanes.rowOf(v.lanes, r.laneStart + k)];
      let wedge = -1;
      for (const c of clades as SortedArray)
        if (c <= leaf && leaf < tree.cladeEnd[c] && (wedge < 0 || c < wedge)) wedge = c;
      if (wedge >= 0 && wedge === open) continue;
      open = wedge;
      display[k] = lane.length;
      lane.push(wedge >= 0 ? -1 : r.laneStart + k);
      node.push(wedge);
    }
    Object.assign(r, { lane: Int32Array.from(lane), node: Int32Array.from(node), display, count: lane.length });
    return r;
  }

  /** Wrapped: each block a prepared run laid out to the data width; a band per visible track per line. */
  private wrapped(e: SectionEntry, items: number[]): number {
    const l = e.layout;
    const model = e.model as ViewModel;
    const v = model.resolved;
    const ax = e.axis as Axis;
    const p = e.props;
    if (!e.runs) {
      const runs: PreparedRun[] = [];
      const labeler = p.wrappedGroupLabels !== "none" ? this.labeler(e) : null;
      for (let b = 0; b < ax.count; b++) runs.push(prepareRun(b, ax, ax.fixedCell, p.wrappedGroup, labeler));
      e.runs = runs;
      e.runsText = this.host.text.version.value;
      this.stats.prepareRuns += runs.length;
    }
    // No lines before the viewport has a width (one position per line would be a line per position).
    const width = l.dataRight - l.dataLeft;
    const n = width > 0 ? layoutRuns(e.runs, width, p.wrappedGroupGap, l.lines) : 0;
    if (width > 0) this.stats.layoutRuns++;
    // A track that draws nothing wrapped (a ruler, omitted layers) takes no entries; `aboveEnd` puts a band for the
    // group labels above each line.
    const all = items.filter((t) => t >= 0);
    const above = p.wrappedGroupLabels === "aboveEnd" ? groupLabelBand(p.wrappedFontPx) : 0;
    if (e.bands.length < all.length * n) e.bands = new Float32Array(2 * all.length * n);
    e.view.axis = ax;
    e.view.runs = e.runs;
    // The lanes' bands read these (`laneLineHeights`).
    e.view.laneHeight = p.wrappedLaneHeight;
    e.view.trackGap = p.wrappedTrackGap;
    const tracks: number[] = [];
    for (const t of all) if (this.bands(e, t, n, tracks.length * n)) tracks.push(t);
    let y = p.wrappedInset;
    for (let k = 0; k < n; k++) {
      if (k > 0) y += l.lines.block[k] !== l.lines.block[k - 1] ? p.stackGap : p.wrappedGap;
      l.lines.top[k] = y;
      if (above > 0) y += push(l, EntryKind.labels, -1, -1, -1, k, l.lines.block[k], y, above, 0);
      for (let j = 0; j < tracks.length; j++) {
        const t = tracks[j];
        y += push(
          l,
          EntryKind.track,
          t,
          v.tracks.group[t],
          firstLayer(v, t),
          k,
          l.lines.block[k],
          y,
          e.bands[j * n + k],
          0,
        );
      }
    }
    return y + p.wrappedInset;
  }

  /**
   * The group labels: a group's last position (above the line) or first (inline) in the block's numbering through its
   * path, as the hover names it (insertion codes as they are; none where the numbering has no label, and none on a
   * ligand); ordinal without a numbering. Widths from the mono font's advances, measured once per font.
   */
  private labeler(e: SectionEntry): GroupLabeler {
    const { resolved: v, set } = e.model as ViewModel;
    const ax = e.axis as Axis;
    const b0 = v.sections.blockStart[e.index];
    const tm = this.host.text;
    const font = tm.font("mono", groupLabelPx(e.props.wrappedFontPx));
    const pad = e.props.wrappedGroupLabels === "inlineStart" ? GROUP_LABEL_PAD : 0;
    return {
      inline: pad > 0,
      text: (b, k) => positionLabel(v, set, ax, b0, b, k),
      width(text) {
        const a = tm.advances(font);
        let w = 0;
        for (let c = 0; c < text.length; c++) {
          const code = text.charCodeAt(c);
          if (code >= 128) return tm.label(text, font).width + pad;
          w += a[code];
        }
        return w + pad;
      },
    };
  }

  /**
   * A track's band per line into `e.bands` from `at`: its layers' `lineHeights` (the tallest, from each line's block
   * coordinate system runs, through the axis view), else letters at the font's size, else a row; false for a track of
   * rulers and omitted layers (they draw nothing wrapped). Into the entry's scratch: nothing per call.
   */
  private bands(e: SectionEntry, t: number, lines: number, at: number): boolean {
    const model = e.model as ViewModel;
    const { resolved: v, set } = model;
    const p = e.props;
    const out = e.bands;
    let fixed = 0;
    let nothing = true;
    out.fill(0, at, at + lines);
    for (let k = firstLayer(v, t); k >= 0 && k < v.layers.count && v.layers.track[k] === t; k++) {
      const r = this.renderer(v, set, k);
      if (r?.lineHeights && lines > 0) {
        if (e.line.length < lines) e.line = new Float32Array(2 * lines);
        e.line.fill(0);
        r.lineHeights(model.prepared[k]?.state, v.layers.params[k], e.view, e.layout.lines, e.line);
        for (let i = 0; i < lines; i++) out[at + i] = Math.max(out[at + i], e.line[i]);
        nothing = false;
        continue;
      }
      const type = v.layers.type[k] >= 0 ? set.strings.get(v.layers.type[k]) : "";
      if (type === "" || type === "ruler") continue;
      nothing = false;
      fixed = Math.max(fixed, type === "letters" ? Math.ceil(p.wrappedFontPx * 1.6) : p.rowHeight);
    }
    // A layer with per-line heights shares the band with the fixed ones (letters and their emphasis).
    for (let i = 0; i < lines; i++) out[at + i] = Math.max(out[at + i], fixed);
    return !nothing;
  }

  /**
   * A track's `height` hint (in rows) wins, raised to the floor a layer's `minTrackHeight` names (a row is never under
   * a device px); else its tallest layer: row sizing under `uniform`, the type's own otherwise (given the data area's
   * width and a box's room, from this build: a resize rebuilds once).
   */
  private trackHeight(e: SectionEntry, t: number): number {
    const model = e.model as ViewModel;
    const { resolved: v, set } = model;
    const p = e.props;
    const hint = v.tracks.height[t];
    if (hint > 0) {
      let floor = hint * p.rowHeight;
      for (let k = firstLayer(v, t); k >= 0 && k < v.layers.count && v.layers.track[k] === t; k++) {
        const fn = v.layers.type[k] >= 0 ? this.renderer(v, set, k)?.minTrackHeight : undefined;
        if (fn) floor = Math.max(floor, fn(model.prepared[k]?.state, v.layers.params[k]));
      }
      return floor;
    }
    const natural = v.sections.rowSizing[e.index] === RowSizing.natural;
    // The data area's width and a box's room for its rows, whole px (ints cross the call unboxed): a plot keeps its aspect.
    const l = e.layout;
    const width = Math.max(0, Math.floor(l.dataRight - l.dataLeft)) | 0;
    const room = l.heightMode === Height.box ? Math.max(0, Math.floor(e.built[1] - l.scrolledTop)) | 0 : 0;
    let h = 0;
    for (let k = firstLayer(v, t); k >= 0 && k < v.layers.count && v.layers.track[k] === t; k++) {
      if (v.layers.type[k] < 0) continue;
      const r = this.renderer(v, set, k);
      const own = natural || v.layers.sizing[k] === Sizing.free;
      const fn = own ? r?.height : undefined;
      h = Math.max(h, fn ? fn(model.prepared[k]?.state, v.layers.params[k], p, e.built[2], width, room) : p.rowHeight);
    }
    return h || p.rowHeight;
  }

  /** Whether a track's height follows the row scale: a layer of it drawn by a renderer that `zooms`, with no height hint. */
  private zooms(e: SectionEntry, t: number): boolean {
    const model = e.model as ViewModel;
    const { resolved: v, set } = model;
    if (v.tracks.height[t] > 0) return false;
    for (let k = firstLayer(v, t); k >= 0 && k < v.layers.count && v.layers.track[k] === t; k++)
      if (v.layers.type[k] >= 0 && this.renderer(v, set, k)?.zooms) return true;
    return false;
  }

  private renderer(v: ResolvedView, set: SeqSet, layer: number): LayerRenderer | undefined {
    const type = v.layers.type[layer];
    return type >= 0 ? this.host.registry.get("renderers", set.strings.get(type)) : undefined;
  }

  /**
   * The visible rows of the `rows` layer, republished when the first or last changes (or after a build). Per move it
   * reads only numbers the build cached on the entry, so it stays one shape however many sections there are.
   */
  private followScroll(e: SectionEntry, s: ViewportState, force: boolean) {
    const h = e.rowSpan[1];
    let d0 = 0;
    let d1 = -1;
    if (h > 0) {
      d0 = Math.max(0, Math.floor((s.scrollY - e.rowSpan[0]) / h));
      d1 = Math.min(e.rowSpan[2] - 1, Math.floor((s.scrollY + s.height - e.rowSpan[3] - e.rowSpan[0]) / h));
    }
    if (!force && d0 === e.shownRows[0] && d1 === e.shownRows[1]) return;
    e.shownRows[0] = d0;
    e.shownRows[1] = d1;
    this.fillRows(e, d0, d1);
  }

  private fillRows(e: SectionEntry, d0: number, d1: number) {
    const r = e.layout.rows;
    const out = e.rows.current;
    const n = Math.max(0, d1 - d0 + 1);
    if (out.entry.length < n) {
      const c = Math.max(n, out.entry.length * 2);
      Object.assign(out, { entry: I(c), displayRow: I(c), dataRow: I(c), top: new Float64Array(c) });
      out.height = new Float32Array(c);
    }
    const lanes = n > 0 ? (e.model as ViewModel).resolved.lanes : null;
    for (let k = 0; lanes && k < n; k++) {
      const d = d0 + k;
      const lane = laneOf(r, d);
      out.entry[k] = entryOfDisplay(r, d);
      out.displayRow[k] = d;
      out.dataRow[k] = lane >= 0 ? Lanes.rowOf(lanes, lane) : -1;
      out.top[k] = e.rowSpan[0] + d * r.height;
      out.height[k] = r.height;
    }
    out.count = n;
    e.rows.bump();
  }

  private close(e: SectionEntry) {
    e.sub?.unsubscribe();
    e.live.complete();
    e.rows.complete();
    e.content.complete();
  }
  dispose() {
    if (this.frame >= 0) this.host.scheduler.cancelFrame(this.frame);
    for (const e of this.map.values()) this.close(e);
    this.drop();
    this.map.clear();
  }
}

/** A pinned item, taken out of the scrolled order. */
const NONE = -1;
const NO_RUNS: readonly PreparedRun[] = [];

/**
 * Shown position k of block b's label (b0: the section's first block): its numbering's through the numbering path, as
 * the hover names it, insertion codes as they are; none on a ligand or where the numbering has none; ordinal without a
 * numbering. The group labels' and the wrapped gutter's, formatted once per group or line.
 */
export function positionLabel(v: ResolvedView, set: SeqSet, ax: Axis, b0: number, b: number, k: number): string {
  if (v.blocks.ligand[b0 + b]) return "";
  const i = ax.indexOf(b, k as ShownIndex);
  const n = v.blocks.numbering[b0 + b];
  if (!n) return String(i + 1);
  const path = v.paths[v.blocks.numberingPath[b0 + b]];
  const j = path ? sourceAt(path, i) : i;
  return j < 0 ? "" : Numbering.format(n, j as never, set.strings);
}

/**
 * Per track its first layer (-1: none), made once per resolved view's `layers.track` (a track's layers are
 * contiguous). A scan per call made a build O(tracks × layers): 1,100 instance tracks cost 13 ms a rebuild.
 */
let firstOf = new Int32Array(0);
let firstOfTrack: Int32Array | null = null;
function firstLayer(v: ResolvedView, t: number): number {
  const track = v.layers.track;
  if (track !== firstOfTrack) {
    firstOfTrack = track;
    if (firstOf.length < v.tracks.count) firstOf = new Int32Array(Math.max(v.tracks.count, 2 * firstOf.length));
    firstOf.fill(-1);
    for (let k = v.layers.count - 1; k >= 0; k--) if (track[k] >= 0 && track[k] < firstOf.length) firstOf[track[k]] = k;
  }
  return t >= 0 && t < v.tracks.count ? firstOf[t] : -1;
}
const typeIs = (v: ResolvedView, set: SeqSet, k: number, name: string) =>
  v.layers.type[k] >= 0 && set.strings.get(v.layers.type[k]) === name;
/** The track's `rows` layer; -1: none. */
function rowsOf(v: ResolvedView, set: SeqSet, t: number): number {
  for (let k = firstLayer(v, t); k >= 0 && k < v.layers.count && v.layers.track[k] === t; k++)
    if (typeIs(v, set, k, "rows")) return k;
  return -1;
}
/** The section's first `rows` layer on a shown track. */
function rowsLayer(v: ResolvedView, set: SeqSet, s: number): number {
  for (let t = 0; t < v.tracks.count; t++) {
    if (v.tracks.section[t] !== s || v.tracks.hidden[t]) continue;
    const k = rowsOf(v, set, t);
    if (k >= 0) return k;
  }
  return -1;
}
/** Pinned: on an alignment axis every track without rows (column data); on blocks, rulers. */
function isPinned(v: ResolvedView, set: SeqSet, t: number): boolean {
  if (v.tracks.group[t] >= 0) return false;
  const columns = v.blocks.alignment[v.sections.blockStart[v.tracks.section[t]]] >= 0;
  return columns ? rowsOf(v, set, t) < 0 : isRuler(v, set, t);
}
/** A track of rulers only. */
function isRuler(v: ResolvedView, set: SeqSet, t: number): boolean {
  let any = false;
  for (let k = firstLayer(v, t); k >= 0 && k < v.layers.count && v.layers.track[k] === t; k++) {
    if (v.layers.type[k] < 0) continue;
    if (!typeIs(v, set, k, "ruler")) return false;
    any = true;
  }
  return any;
}
/** A block's alphabet: its unit's molecule's or its alignment's. */
function alphabetOf(ax: Axis, b: number): Alphabet | null {
  const set = ax.set;
  const a = set.alignments[ax.alignment[b]];
  if (a) return a.alphabet;
  return ax.unit[b] >= 0 ? (set.molecules[set.units.molecule[ax.unit[b]]]?.alphabet ?? null) : null;
}
