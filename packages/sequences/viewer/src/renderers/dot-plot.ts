/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `dot_plot` renderer (seqspec:format#layer-types, seqspec:interaction#hover): positions against positions on the
// section's axis across and its window down the track (draw/dot-plot.ts). Its hit is a cell, the position under x as
// every track's and its partner under y (`Hit.pairUnit`, `pairIndex`), which the hover lights with it; its label adds
// the partner and the cell's value. Linear only: on rings the engine draws its extents as bars and says so.
import {
  type AxisView,
  type DotPlotState,
  DotPlotTiles,
  dotPlotAt,
  dotPlotValue,
  drawDotPlot,
  gridAt,
  type Hit,
  type HoverState,
  kernelEntry,
  type LabelContext,
  type LabelSink,
  overlayDotPlot,
  prepareDotPlot,
  type RenderProps,
} from '../draw/api.js';

const io = new Float64Array(2);
const at = new Int32Array(2);
/** The cell under the point: the hit's position, and the y axis's at the point's height as its partner. */
function hit(v: AxisView, _s: DotPlotState, _p: unknown, out: Hit): boolean {
  const ax = v.axis;
  if (!ax || out.block < 0 || out.index < 0) return false;
  io[1] = out.point[1];
  // On a stacked section y is the lane's block (`onlyBlock`).
  if (dotPlotAt(v, v.top, v.height, io, false, ax.stacked ? v.onlyBlock : -1, at) && ax.unit[at[0]] >= 0) {
    out.pairBlock = at[0];
    out.pairUnit = ax.unit[at[0]];
    out.pairIndex = ax.indexOf(at[0], at[1] as never);
  }
  return true;
}

/** The partner in its block's numbering, and the cell's value under its annotation's label. Edge rate. */
function label(h: HoverState, lc: LabelContext, out: LabelSink, s: DotPlotState): void {
  const set = lc.set;
  const u = h.pairUnit;
  if (u < 0 || !s) return;
  const rv = lc.view.resolved;
  const b0 = rv.sections.blockStart[lc.section];
  const by = h.pairBlock;
  const mol = set.molecules[set.units.molecule[u]];
  if (!mol || by < 0) return;
  const n = rv.blocks.numbering[b0 + by] ?? set.units.numberings[u].items[0] ?? mol.numberings.items[0];
  const i = h.pairIndex;
  const res = mol.element ? '' : String.fromCharCode(mol.code[i]);
  const name = set.strings.get(set.units.label[u]);
  out.line('with', `${name} ${res}${n && i < n.length ? lc.format(n, i as never) : i + 1}`);
  // The cell's value in each grid of the track's dot plots (the hit names the topmost layer).
  for (let k = 0; k < rv.layers.count; k++) {
    const st = rv.layers.track[k] === h.track ? (lc.view.prepared[k]?.state as DotPlotState | undefined) : undefined;
    const at = st?.grids ? gridAt(st, h.block, by) : -1;
    const value = st && at >= 0 ? dotPlotValue(st, h.block, h.index, by, i) : Number.NaN;
    const key = st && at >= 0 ? set.strings.get(set.annotations[st.grids[at].annotation]?.label ?? -1) : '';
    if (!Number.isNaN(value)) out.line(key || 'value', `${Number(value.toPrecision(4))}`);
  }
}

/** Without an `aspect_ratio`: 320 px at a row zoom of 1, following the row zoom. */
const FIXED = 320;
/** With one: never under this, and in a section as tall as its content never over `MAX` (a box's room caps it). */
const MIN = 120;
const MAX = 720;
/** A shrink held (hysteresis) is smaller than this: a page scrollbar's width coming and going. */
const SETTLE = 32;
/**
 * The plot's height, whole px (seqspec:format#layer-types): the data area's width ÷ `aspect_ratio` (1 by default),
 * capped by the section's room (a box's visible rows under its pinned ones, less a margin; else `MAX`) and raised to
 * `MIN`; then at least width ÷ (ratio × `dotPlotMaxStretch`), which wins over both, so a wide screen never flattens the
 * cells past that stretch (the rows scroll). A ratio overrides the row zoom. It follows its target (a resize lands
 * exactly), except a bounce: a shrink back to the target it just grew from, by less than `SETTLE` px and within the
 * cap, keeps the taller height while the target stays there. So a section sized to its content whose plot's shrink
 * hides the page's scrollbar, which widens the plot again, settles after one round instead of flickering. `null` (or
 * no width yet): the fixed height times the row zoom.
 */
function height(s: DotPlotState, p: unknown, props: RenderProps, rowScale = 1, width = 0, room = 0): number {
  const ar = (p as { aspect_ratio?: number | null }).aspect_ratio;
  const ratio = ar === undefined ? 1 : ar;
  if (ratio === null || !(ratio > 0) || width <= 0) return Math.round(FIXED * rowScale) | 0;
  // A box's room less a margin (a tenth of it, at least two rows): a plot that fits shows whole under the pinned rows
  // at more than one scroll position.
  const cap = room > 0 ? room - Math.max(Math.round(room / 10), 2 * props.rowHeight) : MAX;
  const fitted = Math.max(MIN, Math.min(cap, Math.round(width / ratio)));
  const stretch = props.dotPlotMaxStretch > 1 ? props.dotPlotMaxStretch : 1;
  const floor = Math.ceil(width / (ratio * stretch));
  const h = Math.max(fitted, floor) | 0;
  // sized: the height given, the last three distinct targets (newest first), the target a bounce is held at (0: none).
  const z = s.sized;
  if (h !== z[1]) {
    z[3] = z[2];
    z[2] = z[1];
    z[1] = h;
  }
  const shown = z[0];
  const near = shown > h && shown - h < SETTLE && shown <= Math.max(cap, floor);
  // Held: still at a bounce's target. A bounce: the target before the one that set the height, again.
  const held = near && z[4] > 0 && Math.abs(h - z[4]) <= 1;
  const bounce = near && z[2] === shown && Math.abs(h - z[3]) <= 1;
  if (held || bounce) {
    z[4] = h;
    return shown;
  }
  z[4] = 0;
  z[0] = h;
  return h;
}

export const DotPlotRenderer = {
  ...kernelEntry<DotPlotState>('dot_plot', height, {
    content: false,
    zooms: true,
    minHeight: 24,
    miniature: false,
    prepare: prepareDotPlot,
    hit,
    draw: drawDotPlot,
    overlay: overlayDotPlot,
    label,
    // No row band (the hoverRow behavior): its rows are the other axis's positions, which the crosshair lights.
    rowHover: () => false,
    cache: { create: () => new DotPlotTiles(), dispose: (c) => (c as DotPlotTiles).dispose() },
  }),
  geometries: ['linear'] as const,
};
