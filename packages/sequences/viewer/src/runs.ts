/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Wrapped runs (seqspec:viewer#text-measuring): each block of a wrapped section is prepared once into groups of shown
// positions (pretext's segments; the gap between groups its break opportunity), then laid out to a width by arithmetic.
// A
// resize reruns `layoutRuns` only; a font change or a filter reruns `prepareRun`.
//
// Adapted from pretext (https://github.com/chenglou/pretext, commit e699e27e7e4b48a50024fca545c6a12bc6ac975d):
// `countPreparedLines`' greedy loop over segment widths, a segment wider than the line cut by its per-grapheme
// advances (here per position), and the fit tolerance (`lineFitEpsilon`, 0.005 px). Each block starts a new line, but
// chips (seqspec:viewer#the-text-view: a short unit's one position, as wide as its label) flow after each other on
// shared rows, as words do.
//
// pretext is MIT licensed: Copyright (c) 2026 Pretext contributors. Permission is hereby granted, free of charge,
// to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions: The above copyright notice and this permission notice
// shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
// WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import type { Axis, AxisView, ShownIndex } from './axis.js';
import { segWindow } from './pyramids.js';

export interface PreparedRun {
  readonly block: number;
  /** Shown positions. */
  readonly count: number;
  /**
   * Groups of `wrappedGroup` shown positions, also cut at a break (seqspec:model#pieces-and-layout-maps: a filter's, a
   * cyclic origin's).
   */
  readonly segments: number;
  /** Per group: its inline label (`inline`), then its cells × the font's cell, ligand spans included. */
  readonly widths: Float32Array;
  /**
   * Cumulative advance per position, inline labels included: to cut a group wider than a line, and to place. Null
   * without inline labels, where it is `(k + 1) · advance` (doubles, never a float per position).
   */
  readonly positionEnds: Float64Array | null;
  /** Per group, its first shown index. */
  readonly startShown: Int32Array;
  /** A position's advance: its cells × the font's cell. */
  readonly advance: number;
  /**
   * The group labels' placement (`wrappedGroupLabels`): true, inline before each group (its first position's label,
   * its width part of the group's); false, above the line over its last letter (its last position's), taking no room.
   */
  readonly inline: boolean;
  /** Per group, its label in the block's numbering ("": none). */
  readonly labels: readonly string[];
  /** Per group, its label's width; inline, with the space before the letter (0: none). */
  readonly labelWidths: Float32Array;
  /** A chip: a short unit's one position, `advance` its label's room; consecutive chips share rows. */
  readonly chip: boolean;
}

/** A chip's label font px (`ui`, 600): half the letters' band, as a ligand's pill sizes its label in its row. */
export const chipFontPx = (fontPx: number) => Math.max(7, Math.min(12, Math.floor(Math.ceil(fontPx * 1.6) * 0.5)));
/** Px of a chip either side of its label (its pill's border and padding). */
export const CHIP_PAD = 5;

/**
 * The group labels of a wrapped section: which position's label a group carries, its text and width. Asked once per
 * group when a run is prepared, never per frame.
 */
export interface GroupLabeler {
  /** Inline before the group's first letter, else above its last. */
  readonly inline: boolean;
  /** Shown position k's label in the block's numbering; "": none. */
  text(block: number, k: number): string;
  /** Its width; inline, with `GROUP_LABEL_PAD`. */
  width(text: string): number;
}

/** A group label's font px: about 62% of the letters'. */
export const groupLabelPx = (fontPx: number) => Math.max(6, Math.round(fontPx * 0.62));
/** The band the labels above the lines take: their line height and 2 px. */
export const groupLabelBand = (fontPx: number) => Math.ceil(groupLabelPx(fontPx) * 1.25) + 2;
/** Px between an inline label and its group's first letter. */
export const GROUP_LABEL_PAD = 1;

/**
 * Reused: grown only when a layout has more lines than ever before. A line is one block's; `left` is its x from the
 * data area's left (0 but for a chip after another on its row), `follows` 1 where it shares the previous line's row
 * (its `top`). `top` and `bottom` (its row's) are the layout's.
 */
export interface LineTable {
  count: number;
  block: Int32Array;
  start: Int32Array;
  stop: Int32Array;
  left: Float64Array;
  follows: Uint8Array;
  top: Float64Array;
  bottom: Float64Array;
}

export function lineTable(capacity = 16): LineTable {
  return {
    count: 0,
    block: new Int32Array(capacity),
    start: new Int32Array(capacity),
    stop: new Int32Array(capacity),
    left: new Float64Array(capacity),
    follows: new Uint8Array(capacity),
    top: new Float64Array(capacity),
    bottom: new Float64Array(capacity),
  };
}

/**
 * The line of `lines` holding x (`xs[0]`, from the data area's left) on the row of line k: the last starting at or
 * before x. The x in an array: a fractional double passed to a call V8 does not inline is boxed, 16 bytes a hit on a
 * wrapped track (2026-10-15).
 */
export function lineAtX(lines: LineTable, k: number, xs: Float64Array): number {
  const x = xs[0];
  let first = k;
  while (first > 0 && lines.follows[first] === 1) first--;
  let at = first;
  for (let j = first + 1; j < lines.count && lines.follows[j] === 1; j++) if (lines.left[j] <= x) at = j;
  return at;
}

/**
 * One block of a wrapped section. Edge-rate: on a font change, a filter, a new axis, numbering or label placement.
 * With a labeler, each group's label is formatted and measured here; an inline one is part of the group's width.
 */
export function prepareRun(
  block: number,
  axis: Axis,
  cell: number,
  group: number,
  labeler: GroupLabeler | null = null,
  chip = 0,
): PreparedRun {
  const n = axis.shown[block];
  // A chip's one position is its label's room, with no group label.
  if (chip > 0) labeler = null;
  const w = chip > 0 ? chip : axis.span[block] * cell;
  const starts: number[] = [];
  for (let k = 0; k < n; k++)
    if (k === 0 || axis.breakAt(block, k as ShownIndex) || k - starts[starts.length - 1] >= group) starts.push(k);
  const startShown = Int32Array.from(starts);
  const inline = labeler?.inline === true;
  const positionEnds = inline ? new Float64Array(n) : null;
  const widths = new Float32Array(starts.length);
  const labels: string[] = [];
  const labelWidths = new Float32Array(starts.length);
  let x = 0;
  for (let g = 0; g < starts.length; g++) {
    const stop = g + 1 < starts.length ? starts[g + 1] : n;
    const text = labeler ? labeler.text(block, inline ? starts[g] : stop - 1) : '';
    labels.push(text);
    labelWidths[g] = text && labeler ? labeler.width(text) : 0;
    const own = inline ? labelWidths[g] : 0;
    widths[g] = own + (stop - starts[g]) * w;
    x += own;
    for (let k = starts[g]; positionEnds && k < stop; k++) positionEnds[k] = x += w;
  }
  return {
    block,
    count: n,
    segments: starts.length,
    widths,
    positionEnds,
    startShown,
    advance: w,
    inline,
    labels,
    labelWidths,
    chip: chip > 0,
  };
}

/** The advance at the end of shown position k (k = -1: 0). */
const endOf = (run: PreparedRun, k: number) =>
  k < 0 ? 0 : run.positionEnds ? run.positionEnds[k] : (k + 1) * run.advance;

function grow(out: LineTable) {
  const c = out.block.length * 2 || 16;
  const next = lineTable(c);
  next.block.set(out.block);
  next.start.set(out.start);
  next.stop.set(out.stop);
  next.left.set(out.left);
  next.follows.set(out.follows);
  Object.assign(out, { ...next, count: out.count });
}

function push(out: LineTable, n: number, block: number, start: number, stop: number, left = 0, follows = 0): number {
  if (n >= out.block.length) grow(out);
  out.block[n] = block;
  out.start[n] = start;
  out.stop[n] = stop;
  out.left[n] = left;
  out.follows[n] = follows;
  return n + 1;
}

/**
 * Fills lines greedily, each block on a new line; a short last group takes a full group's room; a group wider than a
 * line is cut at a position. Chips flow: one after another on a row, `gap` apart, a new row where the next does not
 * fit. The line count.
 */
export function layoutRuns(runs: readonly PreparedRun[], maxWidth: number, gap: number, out: LineTable): number {
  const fit = Math.max(0, maxWidth) + 0.005;
  let n = 0;
  // The chips' row so far: its width, -1 with none open.
  let row = -1;
  for (let r = 0; r < runs.length; r++) {
    const run = runs[r];
    const { widths, startShown } = run;
    if (run.chip && run.count > 0) {
      const w = run.advance;
      const on = row >= 0 && row + gap + w <= fit;
      n = push(out, n, run.block, 0, run.count, on ? row + gap : 0, on ? 1 : 0);
      row = on ? row + gap + w : w;
      continue;
    }
    row = -1;
    let lineW = 0;
    let has = false;
    let lineStart = 0;
    for (let g = 0; g < run.segments; g++) {
      const w = widths[g];
      if (has) {
        // A short last group joins a line only where a full one would fit: every line holds whole groups alike.
        const room = g === run.segments - 1 && g > 0 && widths[g - 1] > w ? widths[g - 1] : w;
        if (lineW + gap + room <= fit) {
          lineW += gap + w;
          continue;
        }
        n = push(out, n, run.block, lineStart, startShown[g]);
        has = false;
      }
      lineStart = startShown[g];
      if (w <= fit) {
        lineW = w;
        has = true;
        continue;
      }
      // Wider than a line: cut at positions, at least one per line (pretext's overflowing breakable segment).
      const end = g + 1 < run.segments ? startShown[g + 1] : run.count;
      while (lineStart < end) {
        const base = endOf(run, lineStart - 1);
        let lo = lineStart + 1;
        let hi = end;
        while (lo < hi) {
          const mid = (lo + hi + 1) >> 1;
          if (endOf(run, mid - 1) - base <= fit) lo = mid;
          else hi = mid - 1;
        }
        if (lo >= end) {
          lineW = endOf(run, end - 1) - base;
          has = true;
          break;
        }
        n = push(out, n, run.block, lineStart, lo);
        lineStart = lo;
      }
    }
    if (has) n = push(out, n, run.block, lineStart, run.count);
  }
  out.count = n;
  return n;
}

/** The x of shown position k in its line (from the line's start), gaps between groups and inline labels included. */
export function runX(run: PreparedRun, lineStart: number, k: number, gap: number): number {
  runXInto(run, lineStart, k, gap, xs, 0);
  return xs[0];
}
const xs = new Float64Array(1);
/** `runX` into `out[i]`: a per-frame caller takes no double from a call V8 does not inline (it boxes it). */
export function runXInto(run: PreparedRun, lineStart: number, k: number, gap: number, out: Float64Array, i: number) {
  const g = groupOf(run, k);
  // `endOf` inline: a fractional end returned from a call V8 does not inline (a big caller's inlining budget spent, the
  // hit's and the overlay's) is a boxed double, two a call (2026-10-15: 4.7 bytes a step over a select drag's scroll).
  const ends = run.positionEnds;
  const a = k < 1 ? 0 : ends ? ends[k - 1] : k * run.advance;
  const b = lineStart < 1 ? 0 : ends ? ends[lineStart - 1] : lineStart * run.advance;
  let at = a - b;
  if (run.inline && k < run.count && run.startShown[g] === k) at += run.labelWidths[g];
  out[i] = at + gap * (g - groupOf(run, lineStart));
}

/** The width of the line [start, stop) of a run: its last position's right edge from the line's start. */
export function lineWidth(run: PreparedRun, start: number, stop: number, gap: number): number {
  runXInto(run, start, stop - 1, gap, xs, 0);
  return xs[0] + run.advance;
}

/** The group holding shown position k: a binary search over `startShown`. */
export function groupOf(run: PreparedRun, k: number): number {
  const s = run.startShown;
  let lo = 0;
  let hi = run.segments - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (s[mid] <= k) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/**
 * The shown position drawn under the pointer in line [start, stop): x `at[0]` from the line's start, which starts at
 * view x `at[1]`, at DPR `at[2]`. Cell k covers the device px [round((left + x_k)·dpr), round((left + x_k +
 * advance)·dpr)), as the draw places it (`DrawView.xd`, `xe`), so the hit names the cell drawn under the pointer's
 * device px (GG-01). -1: in a gap, on an inline label or past the end. Doubles are read from the record (a fractional
 * double passed to a call V8 does not inline is boxed).
 */
export function linePositionAt(run: PreparedRun, start: number, stop: number, gap: number, at: Float64Array): number {
  const near = linePositionNear(run, start, stop, gap, at);
  const pd = Math.floor((at[1] + at[0]) * at[2]);
  for (let k = near > start ? near - 1 : start; k <= near + 1 && k < stop; k++) {
    runXInto(run, start, k, gap, xs, 0);
    if (pd >= Math.round((at[1] + xs[0]) * at[2]) && pd < Math.round((at[1] + xs[0] + run.advance) * at[2])) return k;
  }
  return -1;
}

/**
 * `linePositionNear` by the drawn cells (`at[1]` the line's view x, `at[2]` the DPR, as `linePositionAt`): the last
 * position whose device-snapped left edge is at or before the pointer's device px, so a select drag takes the cell
 * drawn under it, as a hover does (GG-01, LC-27); left of the line its first position.
 */
export function linePositionDrawnNear(
  run: PreparedRun,
  start: number,
  stop: number,
  gap: number,
  at: Float64Array,
): number {
  const near = linePositionNear(run, start, stop, gap, at);
  const pd = Math.floor((at[1] + at[0]) * at[2]);
  let best = start;
  for (let k = near > start ? near - 1 : start; k <= near + 1 && k < stop; k++) {
    runXInto(run, start, k, gap, xs, 0);
    if (Math.round((at[1] + xs[0]) * at[2]) <= pd) best = k;
  }
  return best;
}

/**
 * The shown position nearest x (`at[0]`, from the line's start) as a text selection takes it: the last position
 * starting at or before x, so a gap or an inline label belongs to the position before it, left of the line its first
 * position and past its end its last.
 */
export function linePositionNear(run: PreparedRun, start: number, stop: number, gap: number, at: Float64Array): number {
  let lo = start;
  let hi = stop - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    runXInto(run, start, mid, gap, xs, 0);
    if (xs[0] <= at[0]) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/**
 * Features prepared over block coordinate system indices in packed lanes (LY2's `Features` has these): per (block,
 * lane), `seg` its first feature, sorted by start, with `reach` (pyramids.ts, `prefixReach`).
 */
export interface LaneSegments {
  readonly lanes: number;
  readonly seg: Int32Array;
  readonly start: Int32Array;
  readonly stop: Int32Array;
  readonly reach: Int32Array;
}

let lineRunsOut = new Int32Array(64);
const win = new Int32Array(2);
/**
 * A wrapped line's lanes: per prepared lane of its block, its lane on the line, compacted to the lanes with an item on
 * it (in lane order), into `out`; -1: none there. Their count. The line's shown bounds become block coordinate system
 * runs through the axis, each searched per lane with `segWindow`, so a line costs
 * its own items. Allocation-free once its scratch has room for the line's runs.
 */
export function lineLanes(f: LaneSegments, axis: AxisView, lines: LineTable, line: number, out: Int32Array): number {
  const b = lines.block[line];
  let n = axis.shownToRuns(b, lines.start[line], lines.stop[line], lineRunsOut);
  if (2 * n > lineRunsOut.length) {
    lineRunsOut = new Int32Array(4 * n);
    n = axis.shownToRuns(b, lines.start[line], lines.stop[line], lineRunsOut);
  }
  let used = 0;
  for (let lane = 0; lane < f.lanes; lane++) {
    const lo = f.seg[b * f.lanes + lane];
    const hi = f.seg[b * f.lanes + lane + 1];
    let on = false;
    for (let r = 0; r < n && !on; r++) {
      const from = lineRunsOut[2 * r];
      segWindow(f.start, f.reach, lo, hi, from, lineRunsOut[2 * r + 1], win);
      for (let k = win[0]; k < win[1] && !on; k++) on = Math.max(f.stop[k], f.start[k] + 1) > from;
    }
    out[lane] = on ? used++ : -1;
  }
  return used;
}

let lanesOut = new Int32Array(32);
/**
 * `lineHeights` of the types drawn as lanes of underlines (`ranges`, `points`, `secondary`): the gap under the letters
 * (`wrappedTrackGap`) and a lane's height (`wrappedLaneHeight`) per lane on the line.
 */
export function laneLineHeights(s: unknown, _p: unknown, axis: AxisView, lines: LineTable, out: Float32Array) {
  const f = s as LaneSegments;
  if (lanesOut.length < f.lanes) lanesOut = new Int32Array(2 * f.lanes);
  for (let i = 0; i < lines.count; i++)
    out[i] = lineLanes(f, axis, lines, i, lanesOut) * axis.laneHeight + axis.trackGap;
}

/** What a strip lane reads of its layer's prepared state (`PerPosition`'s): per block its values or rows, a stream. */
export interface StripSource {
  readonly value: readonly (Float64Array | null)[];
  readonly row: readonly (Int32Array | null)[];
  readonly stream?: unknown;
}
/**
 * `lineHeights` of the types drawn as a strip lane (`values`, `color`; seqspec:viewer#the-text-view): the gap under
 * the letters (`wrappedTrackGap`) and the strip's height (`wrappedStripHeight`) on a line of a block the layer has data
 * on (values, rows or a stream), none on another.
 */
export function stripLineHeights(s: unknown, _p: unknown, axis: AxisView, lines: LineTable, out: Float32Array) {
  const st = s as StripSource | undefined;
  for (let i = 0; i < lines.count; i++) {
    const b = lines.block[i];
    const on = !!st && (st.stream != null || (b < st.value.length && (st.value[b] !== null || st.row[b] !== null)));
    out[i] = on ? axis.stripHeight + axis.trackGap : 0;
  }
}
