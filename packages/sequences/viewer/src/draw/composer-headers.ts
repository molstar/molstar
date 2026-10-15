/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer's headers and gutters (draw/composer.ts): an entry's label (cut with an ellipsis, a group's or a
// split track's chevron), a layer's band labels, the rows' header (its tree and row labels), the data columns' cells
// and heads, each layer's `gutter` and a wrapped line's number. Per frame nothing is allocated (RN1's draw pattern).

import { Lanes } from '@molstar/seq-core/format';
import { type Axis, AxisGeometry } from '../axis.js';
import type { ViewModel } from '../graph.js';
import { CHEVRON_GAP, CHEVRON_HALF, CHEVRON_ROOM, splitChevron } from '../header-toggle.js';
import {
  EntryKind,
  fullTrack,
  positionLabel,
  type SectionLayout,
  SectionLayouts,
  type SplitFeatures,
  trackLabelPx,
  WRAPPED_AXIS_ROOM,
} from '../layout.js';
import type { RowsState } from '../rows.js';
import { fitIn, fitted } from '../text.js';
import { PathMode, PathVerb, SkinToken, TextAlign, TextBaseline } from './backend.js';
import {
  type ComposerFrame,
  depth,
  fail,
  gutterWidth,
  index,
  name,
  open,
  place,
  prime,
  release,
} from './composer-frame.js';

/** A group's or split track's chevron, CSS px around its center: open (⌄), then closed (›). */
const CHEV = new Float64Array([-3.5, -1.5, 0, 1.5, 3.5, -1.5, -1.5, -3.5, 1.5, 0, -1.5, 3.5]);

// Enum-like members read once (a read per frame allocated under the test runner's module wrapper).
const WRAPPED = AxisGeometry.wrapped;
const GROUP = EntryKind.group;
const LABELS = EntryKind.labels;
const TITLE = EntryKind.title;
const TRACK = EntryKind.track;
const STROKE = PathMode.stroke;
const LINE = PathVerb.line;
const MOVE = PathVerb.move;
const TOKEN_INK = SkinToken.ink;
const TOKEN_INK_QUIET = SkinToken.inkQuiet;
const LEFT = TextAlign.left;
const RIGHT = TextAlign.right;
const MIDDLE = TextBaseline.middle;
const BOTTOM = TextBaseline.bottom;

export function header(f: ComposerFrame, l: SectionLayout, e: number) {
  const b = f.box;
  const kind = l.entryKind[e];
  const hw = Math.min(l.headerWidth, l.dataLeft - gutterWidth(f, l));
  if (hw <= 0 || b[1] < 7 || kind === LABELS || !open(f, 0, hw, 0)) return;
  const m = f.model as ViewModel;
  const rv = m.resolved;
  const t = l.entryTrack[e];
  const r = l.rows;
  if (kind === TRACK && r.layer >= 0 && t === r.track && e >= r.entry && !r.lazy) {
    f.rg[0] = b[0] - l.entryRow[e] * b[1];
    f.rg[1] = b[1];
    rowsHeader(f, l, l.entryRow[e], l.entryRow[e]);
    return f.v.backend.popClip();
  }
  let text = '';
  let indent = 0;
  // A split row: its value, indented under the track's label, small.
  const splitAt = kind === TRACK ? SectionLayouts.splitRow(l, e) : -1;
  const split = t >= 0 ? rv.tracks.split[t] : -1;
  if (kind === TITLE) text = name(f, rv.blocks.label[rv.sections.blockStart[index(f)] + l.entryBlock[e]]);
  else if (kind === GROUP) {
    text = name(f, rv.groups.label[l.entryGroup[e]]);
    indent = depth(f, rv.groups.parent[l.entryGroup[e]]);
  } else {
    // An empty label stays empty (a ruler's); none at all: the track's ref.
    text = name(f, rv.tracks.label[t] >= 0 ? rv.tracks.label[t] : rv.tracks.ref[t]);
    indent = depth(f, l.entryGroup[e]) + (splitAt >= 0 ? 1 : 0);
    if (splitAt >= 0) text = (m.prepared[split].state as SplitFeatures).split?.labels[splitAt] ?? '';
  }
  // One size for every label (a split row's a px smaller), smaller only in a band too short for it.
  const px = Math.max(
    splitAt >= 0 ? 7 : 8,
    Math.min(f.v.props.labelFontPx - (splitAt >= 0 ? 1 : 0), Math.floor(b[1] * 0.8)),
  );
  // Regular, as track labels are; only a block's title keeps the heavier weight.
  const font = f.ctx.text.font('ui', px, kind === TITLE ? 600 : 400);
  // A group's chevron after its label; a split track's only on a scrolled entry whose split has values.
  const toggle = kind === GROUP || (split >= 0 && splitChevron(m, l, e));
  const x = 10 + indent * 14;
  const chevroned = toggle && l.entryRow[e] < 0;
  const room = l.headerWidth - x - 6 - (chevroned ? CHEVRON_ROOM : 0);
  if (kind !== TRACK || !bandsHeader(f, l, e, x)) {
    style(f, font);
    at(f, x, b[0] + b[4] / 2, room);
    label(f, text, indent > 0 && kind === TRACK);
    flush(f);
    if (chevroned) {
      const expanded =
        kind === GROUP ? f.ctx.ui.expanded.value.has(name(f, rv.groups.ref[l.entryGroup[e]])) : l.entryRow[e] === -2;
      const w = f.ctx.text.label(text, font).width;
      at(f, x + Math.min(w, room) + CHEVRON_GAP + CHEVRON_HALF, b[0] + b[4] / 2, 0);
      chevron(f, expanded ? 0 : 6);
    }
  }
  f.v.backend.popClip();
}
/** A layer drawing equal bands (a matrix's states): each band's label beside it; false: none, or under 9 px. */
function bandsHeader(f: ComposerFrame, l: SectionLayout, e: number, x: number): boolean {
  const m = f.model as ViewModel;
  const rv = m.resolved;
  const t = l.entryTrack[e];
  let k = l.entryLayer[e];
  while (k >= 0 && k < rv.layers.count && rv.layers.track[k] === t && !f.renderers[k]?.bands) k++;
  const bands = k >= 0 && k < rv.layers.count && rv.layers.track[k] === t ? f.renderers[k]?.bands : undefined;
  const s = bands ? m.prepared[k].state : null;
  const n = bands ? bands.count(s) : 0;
  const b = f.box;
  const bh = n > 0 ? (b[4] - 2) / n : 0;
  if (!bands || bh < 9) return false;
  style(f, f.ctx.text.font('ui', Math.max(8, Math.min(f.v.props.labelFontPx, Math.floor(bh * 0.7)))));
  for (let q = 0; q < n; q++) {
    at(f, x, b[0] + 1 + (q + 0.5) * bh, l.headerWidth - x - 6);
    label(f, bands.label(s, q), false);
  }
  flush(f);
  return true;
}
/**
 * The rows' header over display rows [d0, d1] (`rg`: row 0's top, the row height), in the header clip: the tree
 * beside them (`header.draw`), each row's label (`rowModel`) indented by its depth.
 */
export function rowsHeader(f: ComposerFrame, l: SectionLayout, d0: number, d1: number) {
  const m = f.model as ViewModel;
  const rv = m.resolved;
  const r = l.rows;
  const k = r.layer;
  const rend = f.renderers[k];
  const v = f.v;
  const h = f.rg[1];
  if (!rend || f.failed.has(k)) return;
  const s = m.prepared[k].state;
  const p = rv.layers.params[k];
  const hw = Math.min(l.headerWidth, l.dataLeft - gutterWidth(f, l));
  const inset = rend.header ? Math.min(rend.header.inset(s, p), hw * 0.5) : 0;
  if (rend.header?.draw && inset > 0) {
    prime(f, k, l, r.entry);
    v.top = f.rg[0];
    v.height = h;
    v.left = 0;
    try {
      rend.header.draw(v, s, p, d0, d1);
    } catch (err) {
      fail(f, k, err);
    } finally {
      release(f);
      v.left = l.dataLeft;
    }
  }
  const rm = rend.rowModel;
  if (!rm || h < 7) return;
  style(f, f.ctx.text.font('ui', Math.max(8, Math.min(v.props.labelFontPx, Math.floor(h * 0.8)))));
  for (let d = d0; d <= d1; d++) {
    const lane = SectionLayouts.laneOf(r, d);
    // A collapsed clade's row: its wedge's label, quiet.
    const node = lane < 0 && r.node ? r.node[d] : -1;
    if (lane < 0 && (node < 0 || !rm.wedge)) continue;
    const x = 10 + inset + (lane >= 0 ? Lanes.depthOf(rv.lanes, lane) * 14 : 0);
    const text = lane >= 0 ? rm.label(s, Lanes.rowOf(rv.lanes, lane)) : (rm.wedge?.(s, node) ?? '');
    at(f, x, f.rg[0] + (d + 0.5) * h, l.headerWidth - x - 6);
    label(f, text, lane < 0);
  }
  flush(f);
}
/**
 * The data columns' cells over display rows [d0, d1] (`rg`: row 0's top, the row height), each column clipped to its
 * box: its renderer's `draw` per row (RN3's; each call hands its own batch).
 */
export function columns(f: ComposerFrame, l: SectionLayout, d0: number, d1: number) {
  const r = l.rows;
  const h = f.rg[1];
  const m = f.model as ViewModel;
  const cols = (m.prepared[r.layer]?.state as RowsState | undefined)?.dataColumns;
  const v = f.v;
  for (let c = 0; cols && h >= 7 && c < l.columnCount && c < cols.length; c++) {
    const cr = f.ctx.registry.get('columns', cols[c].renderer);
    if (!cr?.draw || f.failed.has(r.layer) || !open(f, l.columnX[c], l.columnX[c] + l.columnWidth[c], 0)) continue;
    v.columnX = l.columnX[c];
    v.columnWidth = l.columnWidth[c];
    v.height = h;
    try {
      for (let d = d0; d <= d1; d++) {
        const lane = SectionLayouts.laneOf(r, d);
        v.top = f.rg[0] + d * h;
        if (lane >= 0) cr.draw(v, cols[c].state, Lanes.rowOf(m.resolved.lanes, lane));
      }
    } catch (err) {
      fail(f, r.layer, err);
    } finally {
      v.backend.popClip();
    }
  }
}
/**
 * Each data column's label over its cells, centered in the heads band the layout places directly above the first
 * scrolled entry: the last pinned entry's line, or a row of their own. Never the middle of the whole pinned band,
 * which with several pinned tracks lands on another track's line.
 */
export function heads(f: ComposerFrame, l: SectionLayout) {
  const k = l.rows.layer;
  if (l.columnCount === 0 || k < 0 || l.scrolledTop - l.navHeight < 7) return;
  const cols = (f.model as ViewModel).resolved.layers.params[k].data_columns as
    | readonly { readonly label?: string | null; readonly field: string }[]
    | null;
  const cy = l.headsTop + l.headsHeight / 2;
  style(f, f.ctx.text.font('ui', 11, 600));
  for (let c = 0; cols && c < l.columnCount && c < cols.length; c++) {
    at(f, l.columnX[c] + 6, cy, l.columnWidth[c] - 10);
    label(f, cols[c].label ?? cols[c].field, true);
  }
  flush(f);
}
/** A thin chevron centered at `lab`'s point, open (0) or closed (6), in the quiet ink. */
function chevron(f: ComposerFrame, at: number) {
  const p = f.tri;
  const dpr = f.v.dpr;
  const cx = f.lab[0];
  const cy = f.lab[1];
  p.reset();
  p.begin(TOKEN_INK_QUIET, 255, STROKE, 3);
  for (let j = 0; j < 3; j++) {
    const o = p.verb(j === 0 ? MOVE : LINE);
    p.coords[o] = (cx + CHEV[at + 2 * j]) * dpr;
    p.coords[o + 1] = (cy + CHEV[at + 2 * j + 1]) * dpr;
  }
  f.v.backend.paths(p);
}
function style(f: ComposerFrame, font: number) {
  f.text.reset();
  f.text.style(font, LEFT, MIDDLE);
}
/** Where the next label or chevron goes: its x, center y and width (CSS px), through a typed array (no boxing). */
function at(f: ComposerFrame, x: number, cy: number, w: number) {
  f.lab[0] = x;
  f.lab[1] = cy;
  f.lab[2] = w;
}
/** A label cut with an ellipsis to `lab`'s width at its x and center y, on a device pixel, in the batch's font. */
function label(f: ComposerFrame, text: string, quiet: boolean) {
  const tm = f.ctx.text;
  const t = f.text;
  const fit = f.fit;
  const x = f.lab[0];
  const cy = f.lab[1];
  const p = tm.label(text, t.font);
  fit[0] = f.lab[2];
  fit[1] = tm.label('…', t.font).width;
  const n = fitIn(p, fit);
  if (n === 0) return;
  t.begin(quiet ? TOKEN_INK_QUIET : TOKEN_INK);
  const i = t.add(fitted(p, n));
  t.x[i] = Math.round(x * f.v.dpr);
  t.y[i] = Math.round((cy + 0.5) * f.v.dpr);
}
function flush(f: ComposerFrame) {
  if (f.text.n > 0) f.v.backend.text(f.text);
}
/** Each layer's `gutter`; wrapped, the line's first position's label (`none` group labels only). */
export function gutter(f: ComposerFrame, l: SectionLayout, e: number) {
  const v = f.v;
  const b = f.box;
  const m = f.model as ViewModel;
  const rv = m.resolved;
  const t = l.entryTrack[e];
  // Wrapped: the line's number; a track drawn in full, its layers' gutters (an axis).
  if (l.geometry === WRAPPED) {
    lineNumber(f, l, e);
    if (l.entryKind[e] !== TRACK || !fullTrack(v.props, rv, m.set, t) || v.runs[l.lines.block[l.entryRow[e]]]?.chip)
      return;
  }
  let opened = false;
  for (let k = l.entryLayer[e]; k < rv.layers.count && rv.layers.track[k] === t; k++) {
    const r = f.renderers[k];
    if (!r?.gutter || f.failed.has(k)) continue;
    // Clipped to the entry's visible band: a scrolled row's numbers never reach the pinned band.
    if (!opened && !open(f, 0, l.width, 0)) return;
    opened = true;
    v.top = b[0];
    v.height = b[4];
    place(f, l, e);
    prime(f, k, l, e);
    try {
      r.gutter(v, m.prepared[k].state, rv.layers.params[k]);
    } catch (err) {
      fail(f, k, err);
    } finally {
      release(f);
    }
  }
  if (opened) v.backend.popClip();
}
function lineNumber(f: ComposerFrame, l: SectionLayout, e: number) {
  const v = f.v;
  const line = l.entryRow[e];
  if (v.props.wrappedGroupLabels !== 'none' || (e > l.pinnedCount && l.entryRow[e - 1] === line)) return;
  if (f.lineText.length === 0) f.lineText = new Array<string | null>(l.lines.count).fill(null);
  const m = f.model as ViewModel;
  const b0 = m.resolved.sections.blockStart[index(f)];
  const ax = v.axis as Axis;
  const text = f.lineText[line] ?? positionLabel(m.resolved, m.set, ax, b0, l.lines.block[line], l.lines.start[line]);
  f.lineText[line] = text;
  if (!text || !open(f, 0, l.dataLeft, 0)) return;
  const t = f.text;
  t.reset();
  t.style(f.ctx.text.font('mono', Math.round(v.props.wrappedFontPx * 0.8)), RIGHT, MIDDLE);
  t.begin(TOKEN_INK_QUIET);
  const i = t.add(text);
  t.x[i] = Math.round((l.dataLeft - 12 - (v.props.wrappedTrackAxes ? WRAPPED_AXIS_ROOM : 0)) * v.dpr);
  t.y[i] = Math.round((f.box[0] + f.box[4] / 2 + 0.5) * v.dpr);
  v.backend.text(t);
  v.backend.popClip();
}

/**
 * Wrapped (`wrappedTrackLabels`): a track drawn in full labeled in its label's row (`EntryKind.trackLabel`), above its
 * band on the first line it has one, at the line's start, in the quiet ink: never over its data
 * (seqspec:viewer#the-text-view). In the entry's clip (`box`, opened by the caller).
 */
export function trackLabel(f: ComposerFrame, l: SectionLayout, e: number) {
  const v = f.v;
  const m = f.model as ViewModel;
  const id = m.resolved.tracks.label[l.entryTrack[e]];
  const text = id >= 0 ? m.set.strings.get(id) : '';
  if (!text) return;
  const t = f.text;
  t.reset();
  t.style(f.ctx.text.font('ui', trackLabelPx(v.props), 600), LEFT, BOTTOM);
  t.begin(TOKEN_INK_QUIET);
  const i = t.add(text);
  t.x[i] = Math.round((l.dataLeft + l.lines.left[l.entryRow[e]]) * v.dpr);
  t.y[i] = Math.round((f.box[0] + f.box[1] - 1) * v.dpr);
  v.backend.text(t);
}
