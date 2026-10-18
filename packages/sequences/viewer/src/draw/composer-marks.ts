/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer's marks on the overlay (draw/composer.ts): each channel's runs and carets over its block's lane or a
// wrapped line, a split row's own unit's, the hover on rings, the selection being dragged, and the channels' cells
// on an alignment's rows. Per frame nothing is allocated (RN1's draw pattern).

import { Lanes } from '@molstar/seq-core/format';
import { AxisGeometry } from '../axis.js';
import type { ViewModel } from '../graph.js';
import { HoverKind } from '../hit.js';
import { EntryKind, type SectionLayout, SectionLayouts, type SplitFeatures } from '../layout.js';
import { MarkBit, type MarkSet, type SectionMarks } from '../marks.js';
import type { ViewportState } from '../viewport.js';
import { PathMode, PathVerb, SkinToken } from './backend.js';
import { boxOf, type ComposerFrame, index, laneOf, next, paintMark, place, shownRun, tint } from './composer-frame.js';
import { boundaryEdge } from './place.js';

// Enum-like members read once (a read per frame allocated under the test runner's module wrapper).
const WRAPPED = AxisGeometry.wrapped;
const TRACK = EntryKind.track;
const ALIGNED = HoverKind.aligned;
const BIT_BAND = MarkBit.band;
const BIT_HIGHLIGHT = MarkBit.highlight;
const BIT_HOVER = MarkBit.hover;
const BIT_ROW = MarkBit.row;
const BIT_SELECTION = MarkBit.selection;
const STROKE = PathMode.stroke;
const LINE = PathVerb.line;
const MOVE = PathVerb.move;
const TOKEN_HOVER = SkinToken.hover;
const TOKEN_SELECTION_BORDER = SkinToken.selectionBorder;

/** Each channel's runs (selection, highlight, accent, hover) over its block's lane (linear) or per wrapped line. */
export function positionMarks(f: ComposerFrame, l: SectionLayout, s: ViewportState, m: SectionMarks) {
  const v = f.v;
  const ax = v.axis;
  if (!ax) return;
  if (l.geometry !== WRAPPED) {
    for (let b = 0; b < ax.count; b++) if (laneOf(f, l, b)) litRows(f, l, s, m, b, v.visLo[b], v.visHi[b]);
    return;
  }
  for (let e = next(l, s, -1); e >= 0; e = next(l, s, e)) {
    if (l.entryKind[e] !== TRACK || !boxOf(f, l, e)) continue;
    place(f, l, e);
    f.ext[0] = Math.round(f.box[2] * v.dpr) | 0;
    f.ext[1] = Math.round(f.box[3] * v.dpr) | 0;
    lit(f, m, v.onlyBlock, v.lineStart, v.lineStop);
  }
}
/**
 * `lit` over `ext`'s rows, where a split row on a unit of its own draws that unit's marks over its band
 * (`SectionMarks.ownOf`: a run on another such unit leaves it unlit).
 */
function litRows(
  f: ComposerFrame,
  l: SectionLayout,
  s: ViewportState,
  m: SectionMarks,
  b: number,
  lo: number,
  hi: number,
) {
  if (m.ownUnits.length === 0) return lit(f, m, b, lo, hi);
  const dpr = f.v.dpr;
  const y1 = f.ext[1];
  let at = f.ext[0];
  for (let e = next(l, s, -1); e >= 0; e = next(l, s, e)) {
    const own = m.ownOf(rowUnit(f, l, e));
    if (own === null || !boxOf(f, l, e)) continue;
    const top = Math.max(at, Math.round(f.box[2] * dpr) | 0);
    const bottom = Math.min(y1, Math.round(f.box[3] * dpr) | 0);
    if (bottom <= top) continue;
    f.ext[0] = at;
    f.ext[1] = top;
    if (top > at) lit(f, m, b, lo, hi);
    f.ext[0] = top;
    f.ext[1] = bottom;
    lit(f, own, b, lo, hi);
    at = bottom;
  }
  f.ext[0] = at;
  f.ext[1] = y1;
  if (y1 > at) lit(f, m, b, lo, hi);
}
/** Entry e's split row's one unit (`SplitRows.unit`); -1: none. */
function rowUnit(f: ComposerFrame, l: SectionLayout, e: number): number {
  const at = SectionLayouts.splitRow(l, e);
  const m = f.model as ViewModel;
  const k = at >= 0 ? m.resolved.tracks.split[l.entryTrack[e]] : -1;
  const unit = k >= 0 ? (m.prepared[k]?.state as SplitFeatures | undefined)?.split?.unit : undefined;
  return unit !== undefined && at < unit.length ? unit[at] : -1;
}
/**
 * Block b's runs of each channel within shown [lo, hi), over `ext`'s rows. Wrapped, where no band is drawn in the main
 * pass, a `band` set's runs are painted as a highlight's, so a location set shows in every layout (LS2).
 */
function lit(f: ComposerFrame, m: SectionMarks, b: number, lo: number, hi: number) {
  if (f.v.geometry === WRAPPED) litRuns(f, m.runs(BIT_BAND), b, lo, hi, BIT_HIGHLIGHT);
  for (let bit: number = BIT_SELECTION; bit <= BIT_HOVER; bit <<= 1) {
    litRuns(f, m.runs(bit), b, lo, hi, bit);
    litCarets(f, m.runs(bit), b, lo, hi, bit);
  }
}
/**
 * Block b's carets in channel `bit` (seqspec:interaction#channels): a boundary lit, an insertion's place, as a bar on
 * the edge between two shown positions within [lo, hi], over `ext`'s rows and within the data area, in the channel's
 * tint; a selection's with its border's line through it. Only the carets of the indices shown there: per index run
 * of [lo, hi) (`shownToRuns`), a binary search into the block's sorted carets.
 */
function litCarets(f: ComposerFrame, set: MarkSet, b: number, lo: number, hi: number, bit: number) {
  const c0 = set.caretFrom(b);
  const c1 = set.caretTo(b);
  if (c1 <= c0 || hi <= lo || f.ext[1] <= f.ext[0]) return;
  const v = f.v;
  let n = v.shownToRuns(b, lo, hi, f.caretRuns);
  if (2 * n > f.caretRuns.length) {
    f.caretRuns = new Int32Array(4 * n);
    n = v.shownToRuns(b, lo, hi, f.caretRuns);
  }
  const at = set.caretAt;
  const L = f.ext[2];
  const R = f.ext[3];
  const w = Math.max(3, Math.round(3 * v.dpr)) | 0;
  for (let r = 0; r < n; r++) {
    // The first caret at or after the run's start; a caret at its stop is the edge after its last index.
    const s = f.caretRuns[2 * r];
    const t = f.caretRuns[2 * r + 1];
    let k = c0;
    for (let h = c1; k < h; ) {
      const mid = (k + h) >> 1;
      if (at[mid] < s) k = mid + 1;
      else h = mid;
    }
    for (; k < c1 && at[k] <= t; k++) {
      const e = boundaryEdge(v, b, at[k]);
      if (e < lo || e > hi) continue;
      const x = v.xd(b, e);
      const x0 = Math.max(L, x - (w >> 1));
      const x1 = Math.min(R, x - (w >> 1) + w);
      if (x1 <= x0) continue;
      tint(f, bit, 255);
      f.fills.push(x0, f.ext[0], x1 - x0, f.ext[1] - f.ext[0]);
      if (bit !== BIT_SELECTION || x < L || x >= R) continue;
      f.edges.begin(TOKEN_SELECTION_BORDER);
      f.edges.push(x, f.ext[0], Math.max(1, Math.round(v.dpr)) | 0, f.ext[1] - f.ext[0]);
    }
  }
}
function litRuns(f: ComposerFrame, set: MarkSet, b: number, lo: number, hi: number, bit: number) {
  for (let r = set.from(b); r < set.to(b); r++)
    for (let q = 0; shownRun(f, b, set.start[r], set.stop[r], lo, hi, q); q++) {
      if (f.span[1] <= f.span[0]) continue;
      if (bit === BIT_HOVER && f.v.ring) radial(f, b, f.span[0], f.span[1]);
      else paintMark(f, b, f.span[0], f.span[1], f.ext[0], f.ext[1], bit, false);
    }
}
/**
 * The hover on rings: shown [k0, k1) of block b as its sector, or, narrower than 3 device px along the outer ring, a
 * radial line through its middle as wide at the hub as at the rim (the linear hover line's form).
 */
function radial(f: ComposerFrame, b: number, k0: number, k1: number) {
  const v = f.v;
  const x0 = v.xd(b, k0);
  const x1 = v.xd(b, k1);
  if (x1 - x0 >= 3) return paintMark(f, b, k0, k1, f.ext[0], f.ext[1], BIT_HOVER, false);
  const p = f.tri;
  p.reset();
  p.begin(TOKEN_HOVER, 255, STROKE, Math.max(2, Math.round(4 * v.dpr)));
  let o = p.verb(MOVE);
  p.coords[o] = (x0 + x1) / 2;
  p.coords[o + 1] = f.ext[0];
  o = p.verb(LINE);
  p.coords[o] = (x0 + x1) / 2;
  p.coords[o + 1] = f.ext[1];
  v.backend.paths(p);
  p.reset();
}
/**
 * The selection being dragged, from its numbers: a block's shown range, or display rows of the rows (linear);
 * wrapped, the range on each line of its block, over the line's track bands as the selection's marks are.
 */
export function pending(f: ComposerFrame, l: SectionLayout, s: ViewportState) {
  const p = f.ctx.selection.pending.current;
  const v = f.v;
  const r = l.rows;
  if (!p.active || p.section !== index(f) || p.block < 0) return;
  // A drag that removes its range is its outline alone, over the selection it takes from.
  const fill = !p.remove;
  if (l.geometry === WRAPPED) return pendingLines(f, l, s, p.block, p.startShown, p.endShown, fill);
  const k0 = Math.max(v.visLo[p.block], Math.min(p.startShown, p.endShown));
  const k1 = Math.min(v.visHi[p.block], Math.max(p.startShown, p.endShown) + 1);
  if (p.kind === ALIGNED && p.rowStart >= 0 && r.layer >= 0) {
    const top0 = l.scrolledTop + l.entryTop[r.entry] - s.scrollY;
    const y0 = Math.max(l.scrolledTop, top0 + Math.min(p.rowStart, p.rowEnd) * r.height);
    const y1 = Math.min(l.height, top0 + (Math.max(p.rowStart, p.rowEnd) + 1) * r.height);
    const d0 = Math.round(y0 * v.dpr) | 0;
    paintMark(f, p.block, k0, k1, d0, Math.round(y1 * v.dpr) | 0, BIT_SELECTION, true, fill);
  } else if (laneOf(f, l, p.block)) paintMark(f, p.block, k0, k1, f.ext[0], f.ext[1], BIT_SELECTION, false, fill);
  // A dot plot's rectangle: its y range is positions of the axis too, lit across the tracks as the x range is.
  const yb = p.yBlock;
  if (yb < 0 || !laneOf(f, l, yb)) return;
  const y0 = Math.max(v.visLo[yb], Math.min(p.yStartShown, p.yEndShown));
  const y1 = Math.min(v.visHi[yb], Math.max(p.yStartShown, p.yEndShown) + 1);
  paintMark(f, yb, y0, y1, f.ext[0], f.ext[1], BIT_SELECTION, false, fill);
}
/** Wrapped: shown [min(a, z), max(a, z)] of block b on each of its lines in view, over each track band. */
function pendingLines(
  f: ComposerFrame,
  l: SectionLayout,
  s: ViewportState,
  b: number,
  a: number,
  z: number,
  fill: boolean,
) {
  const v = f.v;
  const k0 = Math.min(a, z);
  const k1 = Math.max(a, z) + 1;
  for (let e = next(l, s, -1); e >= 0; e = next(l, s, e)) {
    if (l.entryKind[e] !== TRACK || l.lines.block[l.entryRow[e]] !== b || !boxOf(f, l, e)) continue;
    place(f, l, e);
    const lo = Math.max(k0, v.lineStart);
    const hi = Math.min(k1, v.lineStop);
    if (hi <= lo) continue;
    const y0 = Math.round(f.box[2] * v.dpr) | 0;
    paintMark(f, b, lo, hi, y0, Math.round(f.box[3] * v.dpr) | 0, BIT_SELECTION, false, fill);
  }
}
/** The rows layer's data row → display row (-1: not shown), rebuilt per layout version (edge rate). */
function displayRows(f: ComposerFrame, l: SectionLayout, version: number): Int32Array {
  if (f.displayAt === version) return f.display;
  f.displayAt = version;
  const r = l.rows;
  const lanes = (f.model as ViewModel).resolved.lanes;
  let n = 0;
  for (let d = 0; d < r.count; d++) {
    const lane = SectionLayouts.laneOf(r, d);
    if (lane >= 0) n = Math.max(n, Lanes.rowOf(lanes, lane) + 1);
  }
  if (f.display.length < n) f.display = new Int32Array(n);
  f.display.fill(-1);
  for (let d = 0; d < r.count; d++) {
    const lane = SectionLayouts.laneOf(r, d);
    if (lane >= 0) f.display[Lanes.rowOf(lanes, lane)] = d;
  }
  return f.display;
}
/**
 * The channels' rectangles in a section of an alignment's rows, under the position marks: cells over their rows, a
 * selection outlined; a row lit whole (`MarkBit.row`) across header and data columns too.
 */
export function rowMarks(f: ComposerFrame, l: SectionLayout, m: SectionMarks, version: number) {
  const r = l.rows;
  const rm = m.rows;
  if (r.layer < 0 || rm.count === 0) return;
  const v = f.v;
  const map = displayRows(f, l, version);
  const top0 = l.scrolledTop + l.entryTop[r.entry] - v.scrollY;
  const T = Math.round(l.scrolledTop * v.dpr) | 0;
  const B = Math.round(l.height * v.dpr) | 0;
  for (let q = 0; q < rm.count; q++) {
    const bit = rm.bit[q] & ~BIT_ROW;
    const b = rm.block[q];
    for (let row = rm.rowStart[q]; row < rm.rowStop[q]; ) {
      const d = row < map.length ? map[row] : -1;
      // Rows on consecutive display rows: one rectangle.
      let n = 1;
      while (d >= 0 && row + n < rm.rowStop[q] && row + n < map.length && map[row + n] === d + n) n++;
      row += n;
      const y0 = Math.max(T, Math.round((top0 + d * r.height) * v.dpr)) | 0;
      const y1 = Math.min(B, Math.round((top0 + (d + n) * r.height) * v.dpr)) | 0;
      if (d < 0 || y1 <= y0) continue;
      if (rm.bit[q] & BIT_ROW) {
        tint(f, bit, 255);
        f.fills.push(0, y0, f.ext[2], y1 - y0);
      }
      for (let k = 0; shownRun(f, b, rm.start[q], rm.stop[q], v.visLo[b], v.visHi[b], k); k++)
        if (f.span[1] > f.span[0]) paintMark(f, b, f.span[0], f.span[1], y0, y1, bit, true);
    }
  }
}
