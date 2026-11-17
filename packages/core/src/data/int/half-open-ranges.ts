/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Sets of integers as sorted, merged, half-open `[start, end)` ranges in one `Int32Array`
 * (`[start0, end0, start1, end1, …]`): disjoint, not adjacent, none empty. Unlike `SortedRanges` (inclusive
 * `[min, max]` pairs of a `SortedArray`), the ends are exclusive and the values are int32 (an `Int32Array` would wrap
 * others, so they are refused).
 */
type HalfOpenRanges = Int32Array;

namespace HalfOpenRanges {
  export const Empty: HalfOpenRanges = new Int32Array(0);

  /** Whether a number is an integer within int32. */
  export function isInt32(n: number): boolean {
    return Number.isInteger(n) && n >= -0x80000000 && n <= 0x7fffffff;
  }

  /**
   * Sorts and merges `[start, end, …]` pairs (adjacent ones too), dropping empty ones; a new array. A value that is
   * not an int32 integer is a programmer error and throws, never wraps.
   */
  export function normalize(ranges: ArrayLike<number>): HalfOpenRanges {
    const n = ranges.length >> 1;
    if (ranges.length % 2 !== 0) throw new Error('Ranges are [start, end] pairs');
    const order: number[] = [];
    for (let i = 0; i < n; i++) {
      const s = ranges[2 * i];
      const e = ranges[2 * i + 1];
      if (!isInt32(s) || !isInt32(e)) throw new Error('A range end is not an int32 integer');
      if (e > s) order.push(i);
    }
    order.sort((x, y) => ranges[2 * x] - ranges[2 * y]);
    const out: number[] = [];
    for (const i of order) {
      const s = ranges[2 * i];
      const e = ranges[2 * i + 1];
      const last = out.length - 1;
      if (last > 0 && s <= out[last]) out[last] = Math.max(out[last], e);
      else out.push(s, e);
    }
    return Int32Array.from(out);
  }

  /**
   * Combines two normalized range sets by a predicate on (in a, in b). `keep(false, false)` must be false: positions in
   * neither set are never kept (no complement), so only union-, intersection- and difference-like predicates apply.
   */
  export function combine(
    a: HalfOpenRanges,
    b: HalfOpenRanges,
    keep: (inA: boolean, inB: boolean) => boolean,
  ): HalfOpenRanges {
    const cuts: number[] = [];
    for (let i = 0; i < a.length; i++) cuts.push(a[i]);
    for (let i = 0; i < b.length; i++) cuts.push(b[i]);
    cuts.sort((x, y) => x - y);
    const out: number[] = [];
    let ia = 0;
    let ib = 0;
    for (let k = 0; k + 1 < cuts.length; k++) {
      const s = cuts[k];
      const e = cuts[k + 1];
      if (e === s) continue;
      while (ia < a.length && a[ia + 1] <= s) ia += 2;
      while (ib < b.length && b[ib + 1] <= s) ib += 2;
      const inA = ia < a.length && a[ia] <= s;
      const inB = ib < b.length && b[ib] <= s;
      if (!keep(inA, inB)) continue;
      const last = out.length - 1;
      if (last > 0 && out[last] === s) out[last] = e;
      else out.push(s, e);
    }
    return Int32Array.from(out);
  }

  export function union(a: HalfOpenRanges, b: HalfOpenRanges): HalfOpenRanges {
    return combine(a, b, (x, y) => x || y);
  }

  export function intersect(a: HalfOpenRanges, b: HalfOpenRanges): HalfOpenRanges {
    return combine(a, b, (x, y) => x && y);
  }

  /** The integers of `a` not in `b`. */
  export function subtract(a: HalfOpenRanges, b: HalfOpenRanges): HalfOpenRanges {
    return combine(a, b, (x, y) => x && !y);
  }

  /** How many integers the ranges hold. */
  export function size(ranges: HalfOpenRanges): number {
    let n = 0;
    for (let i = 0; i < ranges.length; i += 2) n += ranges[i + 1] - ranges[i];
    return n;
  }

  /**
   * Whether a half-open item `[start, end)` overlaps the ranges. A zero-length item (a boundary) at `p` overlaps
   * `[s, e)` when `s ≤ p < e`; one at `length` (when given: the end of the domain) overlaps a range ending there.
   * O(log n): a binary search over the ends (the ranges are disjoint and sorted, so their ends are too). Written here,
   * not with `SortedArray`, whose searches read a plain sorted array, not every other entry of interleaved pairs.
   */
  export function overlaps(ranges: HalfOpenRanges, start: number, end: number, length?: number): boolean {
    const n = ranges.length >> 1;
    if (!n) return false;
    if (start === end && start === length) return ranges[2 * n - 1] === length;
    // the first range whose end is past the item's start
    let lo = 0;
    let hi = n;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (ranges[2 * mid + 1] > start) hi = mid;
      else lo = mid + 1;
    }
    if (lo === n) return false;
    const s = ranges[2 * lo];
    return start === end ? s <= start : s < end;
  }

  /** A total order: pairwise by value, then by length. */
  export function compare(a: HalfOpenRanges, b: HalfOpenRanges): number {
    const n = Math.min(a.length, b.length);
    for (let i = 0; i < n; i++) if (a[i] !== b[i]) return a[i] - b[i];
    return a.length - b.length;
  }
}

export { HalfOpenRanges };
