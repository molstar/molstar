/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An axis's offsets and positions (_spec/README.md "Resolved view", Axis): the positions shown are pieces `[start,
// end, …]` (start below end) in display order, laid end to end as display units; on a reversed axis the last piece is
// first and a piece `[s, e)` at offset `o` shows the continuous position `e − (u − o)` at offset `u` (EC-32), forward
// `s + (u − o)`. A break is where two pieces leave positions hidden between them: zero units wide. Lookups are a binary
// search over the pieces' first offsets, kept per pieces array (made once, nothing allocated after).

import { hash2, hashFnv32a } from '@molstar/core/data/util/hash-functions';

/** The axis's shape these read: its pieces in display order, its direction, its system's length and cyclic flag. */
export interface AxisShape {
  readonly ranges?: readonly number[];
  readonly dir: 1 | -1;
  readonly length?: number;
  readonly cyclic: boolean;
}

const Offsets = new WeakMap<readonly number[], Float64Array>();

/** Each piece's first offset, then the total: kept per pieces array. */
export function pieceOffsets(pieces: readonly number[]): Float64Array {
  let o = Offsets.get(pieces);
  if (o) return o;
  const n = pieces.length >> 1;
  o = new Float64Array(n + 1);
  let t = 0;
  for (let k = 0; k < n; k++) {
    o[k] = t;
    t += pieces[2 * k + 1] - pieces[2 * k];
  }
  o[n] = t;
  Offsets.set(pieces, o);
  return o;
}

/** The piece holding offset `u` (the last whose first offset is at most `u`; clamped to the pieces). */
export function pieceAtOffset(offsets: Float64Array, u: number): number {
  const n = offsets.length - 1;
  let lo = 0;
  let hi = n - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (offsets[mid] <= u) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/** The positions shown, or undefined while the axis waits for its system's length. */
export function axisTotal(axis: AxisShape): number | undefined {
  const r = axis.ranges;
  if (!r) return undefined;
  return pieceOffsets(r)[r.length >> 1];
}

/**
 * The continuous position at an offset (clamped to the axis), or undefined while the axis is pending. An offset on a
 * break is the next piece's first edge as a window's `start` and the previous piece's last edge as its `end` (the
 * default), so a window ending on a break shows nothing hidden.
 */
export function positionAt(axis: AxisShape, offset: number, side: 'start' | 'end' = 'end'): number | undefined {
  const r = axis.ranges;
  if (!r || r.length < 2) return undefined;
  const offsets = pieceOffsets(r);
  const n = r.length >> 1;
  const total = offsets[n];
  const u = offset < 0 ? 0 : offset > total ? total : offset;
  let k = pieceAtOffset(offsets, u);
  // on a boundary: the piece before as an end, the one after as a start (a piece of no length is skipped)
  if (side === 'end') while (k > 0 && u === offsets[k]) k--;
  else while (k + 1 < n && u === offsets[k + 1]) k++;
  const d = u - offsets[k];
  return axis.dir === -1 ? r[2 * k + 1] - d : r[2 * k] + d;
}

/**
 * The offset of a continuous position, or undefined when the axis does not show it (or is pending). A position on
 * the edge of two pieces' (contiguous) is the first one's.
 */
export function offsetOf(axis: AxisShape, position: number): number | undefined {
  const r = axis.ranges;
  if (!r) return undefined;
  const k = pieceOfPosition(r, axis.dir, position, true);
  if (k < 0) return undefined;
  const o = pieceOffsets(r)[k];
  return axis.dir === -1 ? o + r[2 * k + 1] - position : o + position - r[2 * k];
}

/** The offset of the left edge of element `p`'s cell in display order, or undefined when it is not shown. */
export function cellOffset(axis: AxisShape, p: number): number | undefined {
  const r = axis.ranges;
  if (!r) return undefined;
  const k = pieceOfPosition(r, axis.dir, p, false);
  if (k < 0) return undefined;
  const o = pieceOffsets(r)[k];
  return axis.dir === -1 ? o + r[2 * k + 1] - 1 - p : o + p - r[2 * k];
}

/** The element at an offset (its cell holds it), or undefined while the axis is pending. */
export function elementAt(axis: AxisShape, offset: number): number | undefined {
  const r = axis.ranges;
  if (!r || r.length < 2) return undefined;
  const offsets = pieceOffsets(r);
  const n = r.length >> 1;
  const total = offsets[n];
  if (total <= 0) return undefined;
  const u = offset < 0 ? 0 : offset >= total ? total - 1e-9 : offset;
  let k = pieceAtOffset(offsets, u);
  while (k + 1 < n && u >= offsets[k + 1]) k++;
  const d = Math.floor(u - offsets[k]);
  return axis.dir === -1 ? r[2 * k + 1] - 1 - d : r[2 * k] + d;
}

/**
 * The piece holding a position (`[s, e)`; with `edge`, `e` too), by binary search when the pieces run one way in
 * display order (always, but over a cyclic origin), else by a scan; −1 when none.
 */
export function pieceOfPosition(r: readonly number[], dir: 1 | -1, p: number, edge: boolean): number {
  const n = r.length >> 1;
  if (monotonic(r, dir)) {
    // the pieces in genomic order: display order forward, reversed backward
    let lo = 0;
    let hi = n - 1;
    while (lo <= hi) {
      const mid = (lo + hi) >> 1;
      const k = dir === -1 ? n - 1 - mid : mid;
      const s = r[2 * k];
      const e = r[2 * k + 1];
      if (p < s) hi = mid - 1;
      else if (p > e || (p === e && !edge)) lo = mid + 1;
      else if (p === e && edge) {
        // contiguous with the next piece in genomic order: the first one's edge
        return k;
      } else return k;
    }
    return -1;
  }
  for (let k = 0; k < n; k++) {
    const s = r[2 * k];
    const e = r[2 * k + 1];
    if (p >= s && (p < e || (edge && p === e))) return k;
  }
  return -1;
}

/** By pieces array, for each direction: index 0 forward, 1 reversed. */
const Monotonic = new WeakMap<readonly number[], [boolean | undefined, boolean | undefined]>();

/**
 * Whether the pieces run one way in display order (genomic order forward, its reverse backward): kept per array and
 * direction.
 */
export function monotonic(r: readonly number[], dir: 1 | -1): boolean {
  let kept = Monotonic.get(r);
  if (!kept) Monotonic.set(r, (kept = [undefined, undefined]));
  const k = dir === 1 ? 0 : 1;
  let m = kept[k];
  if (m === undefined) {
    m = true;
    for (let i = 2; i + 1 < r.length; i += 2) if (dir === 1 ? r[i] < r[i - 1] : r[i + 1] > r[i - 2]) m = false;
    kept[k] = m;
  }
  return m;
}

/** The offset of a position, the nearest piece's edge when the axis does not show it. */
export function nearestOffsetOf(axis: AxisShape, position: number): number | undefined {
  const exact = offsetOf(axis, position);
  if (exact !== undefined) return exact;
  const r = axis.ranges;
  if (!r || r.length < 2) return undefined;
  const offsets = pieceOffsets(r);
  let best = 0;
  let distance = Infinity;
  for (let k = 0; k + 1 < offsets.length; k++) {
    const s = r[2 * k];
    const e = r[2 * k + 1];
    const ds = Math.abs(position - s);
    const de = Math.abs(position - e);
    // a piece's first offset is its start forward, its end reversed
    const atS = axis.dir === -1 ? offsets[k + 1] : offsets[k];
    const atE = axis.dir === -1 ? offsets[k] : offsets[k + 1];
    if (ds < distance) [best, distance] = [atS, ds];
    if (de < distance) [best, distance] = [atE, de];
  }
  return best;
}

/**
 * The breaks of an axis, `[offset, hiddenStart, hiddenEnd, …]` triples in display order: where two pieces leave
 * positions hidden between them (in genomic order, `start > end` across a cyclic origin). Two pieces contiguous, or
 * contiguous across a cyclic origin, have none. Kept per pieces array.
 */
export function axisBreaks(axis: AxisShape): readonly number[] {
  const r = axis.ranges;
  if (!r) return NoBreaks;
  const key = axis.dir === -1 ? ReversedBreaks : ForwardBreaks;
  let b = key.get(r);
  if (b) return b;
  b = [];
  const offsets = pieceOffsets(r);
  const L = axis.length;
  for (let k = 0; k + 3 < r.length; k += 2) {
    // in genomic order, the hidden stretch runs from the lower piece's end to the higher one's start
    const lowEnd = axis.dir === -1 ? r[k + 3] : r[k + 1];
    const highStart = axis.dir === -1 ? r[k] : r[k + 2];
    if (lowEnd === highStart) continue;
    if (axis.cyclic && L !== undefined && lowEnd === L && highStart === 0) continue;
    if (lowEnd > highStart && !(axis.cyclic && L !== undefined)) continue;
    b.push(offsets[(k >> 1) + 1], lowEnd, highStart);
  }
  key.set(r, b);
  return b;
}

const NoBreaks: readonly number[] = [];
const ForwardBreaks = new WeakMap<readonly number[], number[]>();
const ReversedBreaks = new WeakMap<readonly number[], number[]>();

/** How many positions a break hides (across a cyclic origin, both sides of it). */
export function hiddenLength(start: number, end: number, length?: number): number {
  return end >= start ? end - start : (length ?? start) - start + end;
}

/**
 * The break whose hidden stretch holds position `p` (`[start, end)`, across a cyclic origin too): its index in
 * `axisBreaks`' triples (a multiple of 3), or −1.
 */
export function breakHolding(axis: AxisShape, p: number): number {
  const b = axisBreaks(axis);
  for (let i = 0; i + 2 < b.length; i += 3) {
    const s = b[i + 1];
    const e = b[i + 2];
    if (e >= s ? p >= s && p < e : p >= s || p < e) return i;
  }
  return -1;
}

/** Whether two axes show the same positions of the same system, the same way. */
export function axesEqual(a: AxisShape & { readonly system: string }, b: AxisShape & { readonly system: string }) {
  if (a.system !== b.system || a.length !== b.length || a.cyclic !== b.cyclic || a.dir !== b.dir) return false;
  const ra = a.ranges;
  const rb = b.ranges;
  if (!ra || !rb) return ra === rb;
  if (ra.length !== rb.length) return false;
  for (let i = 0; i < ra.length; i++) if (ra[i] !== rb[i]) return false;
  return true;
}

/** A key equal for equal axes: the scroll group of sections on them (EC-113): system, length, pieces and direction. */
export function axisKey(axis: AxisShape & { readonly system: string }): string {
  return `${axis.system}|${axis.length ?? '?'}|${axis.dir}|${axis.ranges ? piecesKey(axis.ranges) : '?'}`;
}

const PiecesKeys = new WeakMap<readonly number[], string>();

/**
 * The pieces' key: their count, ends and two 32-bit hashes, short whatever their number (a history step keeps one;
 * 150,000 pieces joined would be megabytes a step). Kept per array.
 */
function piecesKey(r: readonly number[]): string {
  let k = PiecesKeys.get(r);
  if (k === undefined) {
    let h = 0;
    for (let i = 0; i < r.length; i++) h = hash2(h, r[i]);
    k = `${r.length}:${r[0] ?? ''}:${r[r.length - 1] ?? ''}:${hashFnv32a(r).toString(36)}:${(h >>> 0).toString(36)}`;
    PiecesKeys.set(r, k);
  }
  return k;
}

/**
 * The pieces of a whole system (or its `range`) once its length is known, in display order: one, or two when a cyclic
 * range wraps the origin (reversed: the piece after the origin first).
 */
export function wholePieces(length: number, range: readonly number[] | undefined, dir: 1 | -1, cyclic: boolean) {
  if (range && range.length === 2) {
    const [a, b] = range;
    if (a <= b || !cyclic) return [a, b];
    return dir === -1 ? [0, b, a, length] : [a, length, 0, b];
  }
  return [0, length];
}
