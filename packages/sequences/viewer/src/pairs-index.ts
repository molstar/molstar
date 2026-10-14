/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The pairs layer's indexes in a view (seqspec:viewer#extension-points `pairs`;
// seqspec:sources#levels-and-the-zoom-gate), part of the renderer API and shared by the draw, the overlay and the hit:
// the pairs whose key interval meets the view's keys or holds a position (the interval tree prepare built), and the
// ends within a block run (the end index, by key). A frame costs the pairs and ends it draws, never a scan of every
// pair. Headless and allocation-free: results into the caller's `Int32Array` (room for every pair), keys through module
// scratch so no double crosses a call.

import type { AxisView } from "./axis.js";
import { IndexVisits, overlapping } from "./intervals.js";
import type { OverlayMarks } from "./marks.js";
import type { PairsState } from "./prepare.js";

/** A block's keys start here (`keyOf`: block × 2^32 + index, exact in a double). */
const BLOCK = 4294967296;
const q = new Float64Array(2);

/** The keys of every block drawn into `q`: the lowest and highest over each one's shown hull; false: none drawn. */
function viewKeys(v: AxisView): boolean {
  const ax = v.axis;
  q[0] = Number.POSITIVE_INFINITY;
  q[1] = Number.NEGATIVE_INFINITY;
  for (let b = 0; ax && b < ax.count; b++) {
    const k0 = v.visLo[b];
    const k1 = v.visHi[b];
    if ((v.onlyBlock >= 0 && b !== v.onlyBlock) || k1 <= k0) continue;
    const w = ax.wrapAt(b);
    const i = ax.indexOf(b, k0 as never);
    const j = ax.indexOf(b, (k1 - 1) as never);
    const across = w > k0 && w < k1;
    const lo = b * BLOCK + (across ? ax.base[b] : Math.min(i, j));
    const hi = b * BLOCK + (across ? ax.base[b] + ax.length[b] - 1 : Math.max(i, j));
    if (lo < q[0]) q[0] = lo;
    if (hi > q[1]) q[1] = hi;
  }
  return q[0] <= q[1];
}

/** The pairs to draw here (an end in view, or the view between their ends) into `out`, by first key: their count. */
export function pairsInView(v: AxisView, s: PairsState, out: Int32Array): number {
  if (s.count === 0 || !viewKeys(v)) return 0;
  const n = overlapping(s.lo, s.tree, q, out);
  for (let i = 0; i < n; i++) out[i] = s.order[out[i]];
  return n;
}

/** The pairs whose key interval holds block b's index i (a stab: the arcs over it) into `out`: their count. */
export function pairsAt(s: PairsState, b: number, i: number, out: Int32Array): number {
  q[0] = q[1] = b * BLOCK + i;
  const n = s.count > 0 ? overlapping(s.lo, s.tree, q, out) : 0;
  for (let j = 0; j < n; j++) out[j] = s.order[out[j]];
  return n;
}

/** The first end whose key is at least `q[0]`. */
function lowerEnd(s: PairsState): number {
  let a = 0;
  let b = s.endKey.length;
  while (a < b) {
    const m = (a + b) >> 1;
    if (s.endKey[m] < q[0]) a = m + 1;
    else b = m;
  }
  return a;
}
/**
 * The ends that may meet block b's run [start, stop), as positions [out[0], out[1]) of `s.endOf` (pair × 2 + side):
 * a candidate meets it when it is on block b and its stop passes `start` (an end may start `maxEnd` before the run).
 */
export function endsIn(s: PairsState, b: number, start: number, stop: number, out: Int32Array): void {
  q[0] = b * BLOCK + start - s.maxEnd + 1;
  out[0] = lowerEnd(s);
  q[0] = b * BLOCK + stop;
  out[1] = lowerEnd(s);
}

const ends = new Int32Array(2);
const seen = new WeakMap<PairsState, Int32Array>();
/** Edge-rate, apart: per pair the call that listed it last (the state's call count at the end), so a pair is listed once. */
function seenOf(s: PairsState): Int32Array {
  const x = new Int32Array(s.count + 1);
  seen.set(s, x);
  return x;
}
/**
 * The pairs with a position of either end lit by one of `bits` (`MarkBit` flags: the hover, a highlight, the
 * selection), each once, into `out`: their count. Through each channel's runs and the end index, never every pair.
 */
export function litPairs(v: AxisView, s: PairsState, marks: OverlayMarks, bits: number, out: Int32Array): number {
  if (s.count === 0 || !v.axis) return 0;
  const mark = seen.get(s) ?? seenOf(s);
  if (mark[s.count] === 0x3fffffff) mark.fill(0);
  const stamp = ++mark[s.count];
  let n = 0;
  for (let bit = 1; bit <= bits; bit <<= 1) {
    if ((bits & bit) === 0) continue;
    const set = marks.runs(bit);
    for (let b = 0; b < v.axis.count; b++)
      for (let r = set.from(b); r < set.to(b); r++) {
        endsIn(s, b, set.start[r], set.stop[r], ends);
        for (let j = ends[0]; j < ends[1]; j++) {
          IndexVisits.ends++;
          const e = s.endOf[j];
          const k = e >> 1;
          if (s.block[e] !== b || s.stop[e] <= set.start[r] || mark[k] === stamp) continue;
          mark[k] = stamp;
          out[n++] = k;
        }
      }
  }
  return n;
}

/** The narrowest a glyph of a split cell is (CSS px), and a character of its "+n" marker. */
export const CELL_MIN_GLYPH = 4;
export const CELL_MARKER_CHAR = 5;
/** `io` slots of `cellLayout`: the room and the widest pitch in; the pitch, the kinds hidden and the marker's width out (`at` and `center`: the renderer's scratch). */
export const CellIO = { width: 0, widest: 1, pitch: 2, hidden: 3, marker: 4, at: 5, center: 6, SIZE: 7 } as const;
/** The characters of a "+n" marker, n capped at 99 ("+99"). */
export const markerChars = (hidden: number) => 1 + (hidden > 9 ? 2 : 1);

/**
 * How cell `c` of `s` splits into one glyph per kind present (its kinds are `s.cellKinds[c]`): the number of glyphs
 * shown, side by side, each `io[pitch]` wide, then the "+n" marker `io[marker]` wide when `io[hidden]` kinds do not
 * fit; 0: the cell is not split (one kind, or no room for two glyphs and a marker: it draws as one glyph). The cell's
 * width and the widest pitch (the lone glyph's diameter) are read from `io`; the answers go to `io`, so no double
 * crosses the call.
 */
export function cellLayout(s: PairsState, c: number, io: Float64Array): number {
  const kinds = s.cellKinds[c];
  const room = io[CellIO.width];
  io[CellIO.hidden] = 0;
  io[CellIO.marker] = 0;
  if (kinds < 2 || room < 2 * CELL_MIN_GLYPH) return 0;
  if (kinds * CELL_MIN_GLYPH <= room) {
    io[CellIO.pitch] = Math.min(io[CellIO.widest], room / kinds);
    return kinds;
  }
  // Not all fit: the first ones, and a marker for the rest (its width for the most it could say).
  const marker = CELL_MARKER_CHAR * markerChars(kinds - 1) + 2;
  const shown = Math.min(kinds - 1, Math.floor((room - marker) / CELL_MIN_GLYPH));
  if (shown < 1) return 0;
  io[CellIO.pitch] = Math.min(io[CellIO.widest], (room - marker) / shown);
  io[CellIO.hidden] = kinds - shown;
  io[CellIO.marker] = CELL_MARKER_CHAR * markerChars(kinds - shown) + 2;
  return shown;
}
