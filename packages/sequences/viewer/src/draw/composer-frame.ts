/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The section composer's state and the geometry its passes share (draw/composer.ts): `ComposerFrame` holds every
// field (`SectionComposer` extends it, so one object, one hidden class), and its functions are what the main pass, the
// entries, the headers, the marks and the overlay all use: the entries in view, the clip, an entry's box and its
// clip opened, the view placed on a block, line and row, a layer primed and released, failures, the rows' band, the
// gutter's start, names, a mark's tint and paint, a block's lane, the shown runs of a run, the highlight bands behind
// the entries and the break hairlines. Per frame nothing is allocated: numbers in typed arrays and fields, no double
// passed to a helper (RN1's draw pattern).

import { Lanes } from "@molstar/seq-core/format";
import { type Axis, AxisGeometry, type SectionAxisView, type ShownIndex } from "../axis.js";
import type { ChannelMarks } from "../channels.js";
import type { SeqViewerContext } from "../context.js";
import type { Skin, ViewModel } from "../graph.js";
import { EntryKind, entryBorder, facingOf, type SectionLayout, SectionLayouts, type SplitFeatures } from "../layout.js";
import { MarkBit, type SectionMarks } from "../marks.js";
import type { LayerRenderer } from "../registry.js";
import { Ring } from "../ring.js";
import { RING_WINDOW } from "../ring-window.js";
import type { StreamView } from "../streams.js";
import type { ViewportState } from "../viewport.js";
import { WINDOW_FRAME } from "../windows.js";
import { PathBatch, RectBatch, type RenderBackend, SkinToken, TextBatch } from "./backend.js";
import { Miniature } from "./miniature.js";
import { NavView } from "./nav.js";
import { RingBackend } from "./ring.js";
import { RingLabels } from "./ring-labels.js";
import { DrawView } from "./view.js";

/** An entry's parts; `clear` clears its data box first. */
export const Part = { header: 1, data: 2, gutter: 4, clear: 8 } as const;
export const ALL = Part.header | Part.data | Part.gutter;

// Enum-like members read once (a read per frame allocated under the test runner's module wrapper).
const CIRCULAR = AxisGeometry.circular;
const LINEAR = AxisGeometry.linear;
const WRAPPED = AxisGeometry.wrapped;
const TRACK = EntryKind.track;
const BIT_ACCENT = MarkBit.accent;
const BIT_HOVER = MarkBit.hover;
const BIT_SELECTION = MarkBit.selection;
const PART_CLEAR = Part.clear;
const TOKEN_ACCENT = SkinToken.accent;
const TOKEN_HIGHLIGHT = SkinToken.highlight;
const TOKEN_HOVER = SkinToken.hover;
const TOKEN_INK_FAINT = SkinToken.inkFaint;
const TOKEN_SELECTION = SkinToken.selection;
const TOKEN_SELECTION_BORDER = SkinToken.selectionBorder;

export interface ComposerStats {
  frames: number;
  fullDraws: number;
  copies: number;
  entries: number;
  cacheBuilds: number;
  /** Dense or lazy rows drawn in one `rows` call. */
  rowsCalls: number;
}

/**
 * The composer's state: every field `SectionComposer` and its parts read, initialized here so the composer is one
 * object of one hidden class. Not an API: `/draw` exports `SectionComposer`, and a host reads its documented
 * members only (`stats`, `section`, `live`, `animating`, `again`, `animStep`).
 */
export class ComposerFrame {
  readonly stats: ComposerStats = { frames: 0, fullDraws: 0, copies: 0, entries: 0, cacheBuilds: 0, rowsCalls: 0 };
  readonly section: string;
  /** A surface's (frames follow): skeletons shimmer and arrivals fade in; else (an export) they hold still. */
  live = false;
  /** The last overlay drew a shimmer or a fade: the next overlay frame is needed. */
  animating = false;
  /** The last tracks pass left work for the next frame (`DrawView.more`). */
  again = false;
  /**
   * When `animating`, the ms to the next overlay frame: 0 while a region fades, else `Stream.ShimmerFps`'s period (in an
   * array: a double stored in a field was a number allocated a frame).
   */
  readonly animStep = new Float64Array(1);
  /** The overlay's time (ms), for the progress draws. */
  readonly now = new Float64Array(1);
  readonly ctx: SeqViewerContext;
  /** The pass's view: the main pass's, or the overlay's while it draws (`prime`, `place` and the layers read it). */
  v: DrawView;
  readonly main = new DrawView();
  readonly ov = new DrawView();
  /** The surface's channels (one made on first use unless given) and this section's marks, read once a frame. */
  channels: ChannelMarks | null;
  readonly own: boolean;
  marks: SectionMarks | null = null;
  lay: SectionLayout | null = null;
  primaryAt: Axis | null = null;
  /** The highlight rebuild the main pass's bands were drawn at; whether it drew any. */
  bandAt = -1;
  bandsDrawn = false;
  /** The navigator's view; the window's frame (windows.ts); the axis its block names are of. */
  readonly nview = new NavView();
  readonly g = new Float64Array(WINDOW_FRAME);
  namesAt: Axis | null = null;
  /** The tracks in miniature (`seqv:map`, the navigator's map). */
  readonly mini: Miniature;
  model: ViewModel | null = null;
  skin: Skin | null = null;
  /** By layer: its entry, its stream, the stream version drawn and whether it changed since. */
  renderers: (LayerRenderer | undefined)[] = [];
  refs: string[] = [];
  streams: (StreamView | null)[] = [];
  streamAt = new Float64Array(0);
  fresh = new Uint8Array(0);
  /** Per layer, its last draw left work (`DrawView.more`): its entries drawn again whole next frame. */
  pend = new Uint8Array(0);
  /** Per layer, the columns [lo, hi) its stream's changes cover (-1: all); `changes` takes a stream's pairs. */
  cols = new Int32Array(0);
  readonly changes = new Int32Array(64);
  readonly failed = new Set<number>();
  /** Layers reported as declaring a geometry and having no `draw`, once each until the model changes. */
  readonly undrawn = new Set<number>();
  /** Layers drawn as bars on rings (`ringless`), reported once each until the model changes. */
  readonly unringed = new Set<number>();
  readonly caches = new Map<unknown, { r: LayerRenderer; c: unknown }>();
  /** At the last draw: scrollX, scrollY, cell, rowScale, layout version, dpr, width, height (device px). */
  readonly drawn = new Float64Array(8).fill(Number.NaN);
  /** The navigator as drawn: layout version, width, dpr. */
  readonly navAt = new Float64Array(3).fill(Number.NaN);
  /**
   * The clip of the parts drawn now (device px: x0, x1, y0, y1), and the entry's top, height, visible y0, y1 and its
   * content's height: the height less the row border drawn inside its box (`SectionLayouts.entryBorder`), so what is
   * drawn on the row centers between its lines.
   */
  readonly clip = new Int32Array(4);
  readonly box = new Float64Array(5);
  /** The last clip `open` pushed: x0, x1 (device px). */
  readonly opened = new Int32Array(2);
  /** A label's x, center y and width (CSS px); the fit's width and ellipsis; the rows band's row-0 top and height. */
  readonly lab = new Float64Array(3);
  readonly fit = new Float64Array(2);
  readonly rg = new Float64Array(2);
  /** A horizontal strip's shown range per block, in place of the viewport's. */
  stripLo = new Int32Array(0);
  stripHi = new Int32Array(0);
  runs = new Int32Array(16);
  /** A run's shown runs (at most two: a cyclic origin); a lane's rows and the data area's x (device px); a span. */
  readonly shown = new Int32Array(4);
  shownCount = 0;
  /** `litCarets`' index runs of the shown range (start, stop pairs), grown as needed. */
  caretRuns = new Int32Array(16);
  readonly ext = new Int32Array(4);
  readonly span = new Int32Array(2);
  /** A renderer's `rowHover` band (device rows). */
  readonly bandRows = new Int32Array(2);
  /** Side by side: where the entries drawn end (CSS px), the marks' bottom. */
  bottom = 0;
  /** The rows layer's data row → display row (-1: not shown), at layout version `displayAt`. */
  display = new Int32Array(0);
  displayAt = -1;
  /** The cursor's rect (`cursorRect`), and a hovered tree node (-1: none) whose clade the overlay draws. */
  readonly cursor = new Float64Array(6);
  /** The circular geometry's map and backend. */
  readonly ring = new Ring();
  readonly ringBe = new RingBackend();
  readonly ringAt = new Float64Array(2).fill(Number.NaN);
  /** The ring's window on its linked linear section (ring-window.ts). */
  readonly rw = new Float64Array(RING_WINDOW);
  readonly rb = new Float64Array(2);
  /** The ring's item labels (ring-labels.ts), placed at a layout version, DPR and model. */
  readonly ringLabels = new RingLabels();
  readonly labelsAt = new Float64Array(2).fill(Number.NaN);
  labelsOf: ViewModel | null = null;
  labelsHidden = false;
  clade = -1;
  /** The dense rows' band was drawn in this region (one `rows` call covers every dense entry). */
  dense = false;
  version = -1;
  lineText: (string | null)[] = [];
  readonly text = new TextBatch(16);
  readonly rects = new RectBatch(16);
  readonly tri = new PathBatch(4);
  /** The overlay's fills and outlines (outlines over every fill). */
  readonly fills = new RectBatch(64);
  readonly edges = new RectBatch(16);

  constructor(ctx: SeqViewerContext, section: string, channels?: ChannelMarks) {
    this.ctx = ctx;
    this.section = section;
    this.v = this.main;
    this.mini = new Miniature(ctx, section);
    this.main.tracks = this.ov.tracks = this.mini;
    // The overlay reads the density mode the main pass drew (`densityShown`, seqspec:interaction#hover).
    this.ov.densityModes = this.main.densityModes;
    // A dot plot's overlay draws the rectangle being dragged (seqspec:interaction#selection).
    this.ov.pending = ctx.selection.pending.current;
    this.channels = channels ?? null;
    this.own = !channels;
  }
}

export function view(f: ComposerFrame, v: DrawView, be: RenderBackend, av: SectionAxisView) {
  v.from(av);
  // On rings a draw is a linear one along the unrolled axis; `ring` says where it lands.
  v.ring = av.geometry === CIRCULAR ? f.ring : null;
  if (v.ring) {
    v.geometry = LINEAR;
    f.ring.drawn = v;
  }
  v.backend = be;
  v.skin = f.skin as Skin;
  v.props = f.ctx.layout.props(f.section);
  v.text = f.ctx.text;
}
/** The entries in view: pinned, then scrolled; -1 past the last. */
export function next(l: SectionLayout, s: ViewportState, e: number): number {
  let n = e + 1;
  if (n >= l.pinnedCount && n < s.firstEntry) n = s.firstEntry;
  return n < l.entryCount && (n < l.pinnedCount || n <= s.lastEntry) ? n : -1;
}
export function setClip(f: ComposerFrame, x0: number, x1: number, y0: number, y1: number) {
  const c = f.clip;
  c[0] = x0;
  c[1] = x1;
  c[2] = y0;
  c[3] = y1;
}
/** Entry e's top and height, and its visible band within the clip, into `box` (CSS px); false: none of it shows. */
export function boxOf(f: ComposerFrame, l: SectionLayout, e: number): boolean {
  const v = f.v;
  const b = f.box;
  const c = f.clip;
  const pinned = e < l.pinnedCount;
  const top = pinned ? l.entryTop[e] : l.scrolledTop + l.entryTop[e] - v.scrollY;
  b[0] = top;
  b[1] = l.entryHeight[e];
  b[2] = Math.max(top, pinned ? l.navHeight : l.scrolledTop, c[2] / v.dpr);
  b[3] = Math.min(top + b[1], pinned ? l.scrolledTop : l.height, c[3] / v.dpr);
  // The border's device px, inline (a double passed to a call is boxed).
  const w = Math.max(1, Math.round(v.dpr)) | 0;
  b[4] = l.entryKind[e] === TRACK && entryBorder(l, e, v.props.rowBorders) ? b[1] - w / v.dpr : b[1];
  return b[3] > b[2];
}
/** A clip of CSS x [x0, x1) over the entry's visible band, within the current clip; false: empty, nothing pushed. */
export function open(f: ComposerFrame, x0: number, x1: number, parts: number): boolean {
  const c = f.clip;
  const dpr = f.v.dpr;
  const a = Math.max(c[0], Math.round(x0 * dpr)) | 0;
  const z = Math.min(c[1], Math.round(x1 * dpr)) | 0;
  const y0 = Math.max(c[2], Math.round(f.box[2] * dpr)) | 0;
  const y1 = Math.min(c[3], Math.round(f.box[3] * dpr)) | 0;
  if (z <= a || y1 <= y0) return false;
  const o = f.opened;
  o[0] = a;
  o[1] = z;
  f.v.backend.pushClip(a, y0, z - a, y1 - y0);
  if (parts & PART_CLEAR) {
    f.v.backend.clear(a, y0, z - a, y1 - y0);
    bands(f, f.lay as SectionLayout);
  }
  return true;
}
/** The view's block, line and row for entry e. */
export function place(f: ComposerFrame, l: SectionLayout, e: number) {
  const v = f.v;
  const r = l.rows;
  v.dataRow = v.row = -1;
  if (r.layer >= 0 && l.entryTrack[e] === r.track && e >= r.entry && !r.lazy) onRow(f, l, l.entryRow[e]);
  if (l.geometry === WRAPPED) {
    const line = l.entryRow[e];
    v.line = line;
    v.onlyBlock = l.lines.block[line];
    v.lineStart = l.lines.start[line];
    v.lineStop = l.lines.stop[line];
  } else v.onlyBlock = v.axis?.stacked ? l.entryBlock[e] : -1;
}
/** The view on display row d of the rows layer: its lane offset and data row (-1: a wedge). */
export function onRow(f: ComposerFrame, l: SectionLayout, d: number) {
  const lane = SectionLayouts.laneOf(l.rows, d);
  f.v.row = lane >= 0 ? lane - l.rows.laneStart : -1;
  f.v.dataRow = lane >= 0 ? Lanes.rowOf((f.model as ViewModel).resolved.lanes, lane) : -1;
}
/** Layer k's prepared state on entry e: an expanded split row's own. */
export function stateOf(f: ComposerFrame, k: number, l: SectionLayout, e: number): unknown {
  const at = SectionLayouts.splitRow(l, e);
  const s0 = (f.model as ViewModel).prepared[k].state;
  return at >= 0 ? ((s0 as SplitFeatures).split?.rows[at] ?? s0) : s0;
}
/** The view's layer state: its index, themes, cache, stream and facing. */
export function prime(f: ComposerFrame, k: number, l: SectionLayout, e: number) {
  const v = f.v;
  const m = f.model as ViewModel;
  const r = f.renderers[k] as LayerRenderer;
  v.layer = k;
  v.setThemes(m.prepared[k].themes);
  v.cache = r.cache ? cacheOf(f, r, m.prepared[k].state) : null;
  v.stream = f.streams[k];
  v.facing = v.axis ? facingOf(l, v.axis, k, m.resolved.layers.track[k], e) : 0;
}
export function release(f: ComposerFrame) {
  f.v.cache = f.v.stream = null;
}
function cacheOf(f: ComposerFrame, r: LayerRenderer, s: unknown): unknown {
  const held = f.caches.get(s);
  if (held !== undefined) return held.c;
  const c = r.cache?.create(s, f.skin as Skin);
  f.caches.set(s, { r, c });
  f.stats.cacheBuilds++;
  return c;
}
export function fail(f: ComposerFrame, k: number, err: unknown) {
  f.failed.add(k);
  const m = f.model as ViewModel;
  const ref = m.set.strings.get(m.resolved.layers.ref[k]);
  const message = err instanceof Error ? err.message : String(err);
  f.ctx.report("seqv.renderer.error", `layer "${ref}" threw while drawing and is omitted: ${message}`);
}
/** A rows entry drawn as one band: lazy rows, or rows thinner than `denseRowHeight`. */
export function banded(f: ComposerFrame, l: SectionLayout, e: number): boolean {
  const r = l.rows;
  return r.layer >= 0 && l.entryTrack[e] === r.track && e >= r.entry && (r.lazy || r.height < f.v.props.denseRowHeight);
}
/** Where the gutter starts: the header's and data columns' end. */
export function gutterWidth(f: ComposerFrame, l: SectionLayout): number {
  return l.geometry === WRAPPED ? l.dataLeft - l.headerWidth : f.v.props.gutterLeft;
}
export function name(f: ComposerFrame, id: number): string {
  return id >= 0 ? (f.model as ViewModel).set.strings.get(id) : "";
}
export function depth(f: ComposerFrame, group: number): number {
  const parent = (f.model as ViewModel).resolved.groups.parent;
  let d = 0;
  for (let k = group; k >= 0 && d < 16; k = parent[k]) d++;
  return d;
}
export function index(f: ComposerFrame): number {
  return f.ctx.layout.sectionIndex(f.section);
}
/** Where a filter skips positions or a cyclic origin falls: a full-height hairline, from the shown runs. */
export function breaks(f: ComposerFrame, l: SectionLayout) {
  const v = f.v;
  const ax = v.axis;
  const c = f.clip;
  if (!ax || l.geometry !== LINEAR) return;
  const r = f.rects;
  const w = Math.max(1, Math.round(v.dpr)) | 0;
  r.reset();
  r.begin(TOKEN_INK_FAINT);
  for (let b = 0; b < ax.count; b++) {
    const lo = v.visLo[b];
    let n = v.shownToRuns(b, lo, v.visHi[b], f.runs);
    if (n < 2) continue;
    if (2 * n > f.runs.length) {
      f.runs = new Int32Array(4 * n);
      n = v.shownToRuns(b, lo, v.visHi[b], f.runs);
    }
    // Stacked: over the block's own lane only.
    const e0 = l.blockEntries[2 * b];
    const e1 = l.blockEntries[2 * b + 1] - 1;
    const stacked = ax.stacked && e1 >= e0;
    const top = stacked ? (l.scrolledTop + l.entryTop[e0] - v.scrollY) * v.dpr : l.navHeight * v.dpr;
    const bottom = stacked ? (l.scrolledTop + l.entryTop[e1] + l.entryHeight[e1] - v.scrollY) * v.dpr : c[3];
    const y0 = Math.max(c[2], Math.round(top)) | 0;
    const y1 = Math.min(c[3], Math.round(bottom)) | 0;
    let k = lo;
    for (let q = 0; q < n && y1 > y0; q++) {
      const x = (v.xd(b, k) - (w >> 1)) | 0;
      if (q > 0 && v.breakAt(b, k as ShownIndex) && x >= c[0] && x + w <= c[1]) r.push(x, y0, w, y1 - y0);
      k += Math.abs(f.runs[2 * q + 1] - f.runs[2 * q]);
    }
  }
  if (r.n > 0) v.backend.rects(r);
}
/** The fills' style of channel `bit` (seqspec:interaction#the-skin): the accent's color as faint as a highlight. */
export function tint(f: ComposerFrame, bit: number, alpha: number) {
  const t =
    bit === BIT_SELECTION
      ? TOKEN_SELECTION
      : bit === BIT_HOVER
        ? TOKEN_HOVER
        : bit === BIT_ACCENT
          ? TOKEN_ACCENT
          : TOKEN_HIGHLIGHT;
  f.fills.begin(t, bit === BIT_ACCENT ? alpha >> 2 : alpha);
}
/**
 * One mark: shown [k0, k1) of block b over device rows [y0, y1), the cells' pixels (a wide cell's gap left out, as
 * `cellSpan`) within the data area, in channel `bit`'s tint; a selection's outline over every fill: its sides (a
 * column band) or, `box`, all four (cells on rows), one device px wide at every geometry.
 */
export function paintMark(
  f: ComposerFrame,
  b: number,
  k0: number,
  k1: number,
  y0: number,
  y1: number,
  bit: number,
  box: boolean,
  fill = true,
) {
  const v = f.v;
  const L = f.ext[2];
  const R = f.ext[3];
  const xl = v.xd(b, k0);
  const xr = v.geometry === WRAPPED ? v.xe(b, k1 - 1) : v.xd(b, k1) - (v.cell * v.span[b] > 6 ? 1 : 0);
  if (xr <= L || xl >= R || y1 <= y0) return;
  const x0 = Math.max(xl, L);
  const x1 = Math.max(x0 + 1, Math.min(xr, R));
  if (fill) {
    tint(f, bit, 255);
    f.fills.push(x0, y0, x1 - x0, y1 - y0);
  }
  const w = Math.max(1, Math.round(v.dpr)) | 0;
  const g = f.edges;
  if (bit !== BIT_SELECTION) return;
  if (box && y1 - y0 > 2 * w) {
    g.begin(TOKEN_SELECTION_BORDER, 255, w);
    g.push(x0, y0, x1 - x0, y1 - y0);
  } else if (!box && x1 - x0 > 2 * w) {
    g.begin(TOKEN_SELECTION_BORDER);
    g.push(x0, y0, w, y1 - y0);
    g.push(x1 - w, y0, w, y1 - y0);
  }
}
/** Block b's marks' rows (device px) into `ext`: stacked, its own lane's entries only; else every entry drawn. */
export function laneOf(f: ComposerFrame, l: SectionLayout, b: number): boolean {
  const v = f.v;
  // On rings, from the drawn outer edge (an outermost ruler's ticks, not its numbers' room) to the last ring.
  let y0 = l.geometry === CIRCULAR ? l.ringTop : l.navHeight;
  let y1 = f.bottom;
  if (v.axis?.stacked) {
    const e0 = l.blockEntries[2 * b];
    const e1 = l.blockEntries[2 * b + 1] - 1;
    if (e1 < e0) return false;
    y0 = Math.max(l.scrolledTop, l.scrolledTop + l.entryTop[e0] - v.scrollY);
    y1 = Math.min(l.height, l.scrolledTop + l.entryTop[e1] + l.entryHeight[e1] - v.scrollY);
  }
  f.ext[0] = Math.round(y0 * v.dpr) | 0;
  f.ext[1] = Math.round(y1 * v.dpr) | 0;
  return f.ext[1] > f.ext[0];
}
/** The `q`th shown run of block-system run [s, t) (q = 0 converts it), clipped to [lo, hi), into `span`; false: none. */
export function shownRun(
  f: ComposerFrame,
  b: number,
  s: number,
  t: number,
  lo: number,
  hi: number,
  q: number,
): boolean {
  if (q === 0) f.shownCount = f.v.runToShown(b, s, t, f.shown);
  if (q >= f.shownCount) return false;
  f.span[0] = Math.max(lo, f.shown[2 * q]);
  f.span[1] = Math.min(hi, f.shown[2 * q + 1]);
  return true;
}
/** Highlights of style `band`: a full-height band per run behind the entries, in the clip (linear). */
export function bands(f: ComposerFrame, l: SectionLayout) {
  const m = f.marks;
  const v = f.v;
  const c = f.clip;
  if (!m || m.band.count === 0 || l.geometry === WRAPPED || !v.axis) return;
  f.bottom = l.height;
  const r = f.rects;
  r.reset();
  r.begin(TOKEN_HIGHLIGHT);
  for (let b = 0; b < v.axis.count; b++) {
    if (!laneOf(f, l, b)) continue;
    const y0 = Math.max(c[2], f.ext[0]) | 0;
    const y1 = Math.min(c[3], f.ext[1]) | 0;
    for (let q = m.band.from(b); y1 > y0 && q < m.band.to(b); q++)
      for (let k = 0; shownRun(f, b, m.band.start[q], m.band.stop[q], v.visLo[b], v.visHi[b], k); k++) {
        const x0 = Math.max(c[0], f.ext[2], v.xd(b, f.span[0])) | 0;
        const x1 = Math.min(c[1], f.ext[3], v.xd(b, f.span[1])) | 0;
        if (x1 > x0) r.push(x0, y0, x1 - x0, y1 - y0);
      }
  }
  if (r.n > 0) v.backend.rects(r);
}
