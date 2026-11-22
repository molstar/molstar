/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Intervals found by a region (_spec/README.md "Data", Intervals): an index by start of the short ones and a list of
// the long ones, so one chromosome-long feature never makes every read scan the sequence. Shared by the GFF and VCF
// runtimes and the tooltips' tables (a store entry's rows looked up under a hover). Internal: not in the `data` entry.

/**
 * Intervals found by a region: the short ones sorted by start with the longest of them (one meeting `start` starts at
 * `start − maxShort` or later: a binary search and a scan of what is near), the long ones apart (checked each read:
 * few). Why apart: one chromosome-long feature would make every read scan from the sequence's start (a running
 * greatest end reaches every region past it).
 */
export interface IntervalIndex {
  /** The short intervals' numbers, by start. */
  readonly order: Int32Array;
  readonly starts: Int32Array;
  readonly maxShort: number;
  /** The long intervals' numbers. */
  readonly long: Int32Array;
  /** Every `Block`th start (`starts[k * Block]`): a search reads this small array first, then one block of starts. */
  readonly blocks: Int32Array;
}

/** The starts a block of the index's directory spans (`IntervalIndex.blocks`). */
const Block = 64;

/** A span past which an interval is long: 64 times the median, at least 4,096 positions. */
function longSpan(s: ArrayLike<number>, e: ArrayLike<number>): number {
  if (!s.length) return 4096;
  const spans = Float64Array.from({ length: s.length }, (_, i) => e[i] - s[i]).sort();
  return Math.max(4096, 64 * spans[spans.length >> 1]);
}

export function intervalIndex(s: ArrayLike<number>, e: ArrayLike<number>): IntervalIndex {
  const limit = longSpan(s, e);
  const short: number[] = [];
  const long: number[] = [];
  let maxShort = 1;
  for (let i = 0; i < s.length; i++) {
    if (e[i] - s[i] > limit) long.push(i);
    else {
      short.push(i);
      maxShort = Math.max(maxShort, e[i] - s[i]);
    }
  }
  const order = Int32Array.from(short).sort((a, b) => s[a] - s[b] || a - b);
  const starts = Int32Array.from(order, (i) => s[i]);
  const blocks = Int32Array.from({ length: Math.ceil(starts.length / Block) }, (_, k) => starts[k * Block]);
  return { order, starts, maxShort, long: Int32Array.from(long), blocks };
}

/** Calls `fn` for each interval meeting `[start, end)` (a short one in order of start, then the long ones). */
export function eachMeeting(
  ix: IntervalIndex,
  s: ArrayLike<number>,
  e: ArrayLike<number>,
  start: number,
  end: number,
  fn: (i: number) => void,
) {
  const { order, starts } = ix;
  for (let k = lowerBound(starts.length, (k) => starts[k], start - ix.maxShort); k < starts.length; k++) {
    if (starts[k] >= end) break;
    const i = order[k];
    if (meets(s[i], e[i], start, end)) fn(i);
  }
  for (let k = 0; k < ix.long.length; k++) {
    const i = ix.long[k];
    if (meets(s[i], e[i], start, end)) fn(i);
  }
}

/**
 * The intervals meeting `[start, end)` into `out` (the first `out.length` of them: a short one in order of start, then
 * the long ones), and how many meet in all: no closure, no allocation (a hover's lookup, O(log n + k)).
 */
export function meetingInto(
  ix: IntervalIndex,
  s: ArrayLike<number>,
  e: ArrayLike<number>,
  start: number,
  end: number,
  out: Int32Array,
): number {
  const { order, starts, blocks } = ix;
  let n = 0;
  const v = start - ix.maxShort;
  // the block holding the first start at least `v` (the small array, cached), then within it
  let lo = 0;
  let hi = blocks.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (blocks[mid] < v) lo = mid + 1;
    else hi = mid;
  }
  hi = Math.min(starts.length, lo * Block);
  lo = Math.max(0, (lo - 1) * Block);
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (starts[mid] < v) lo = mid + 1;
    else hi = mid;
  }
  for (let k = lo; k < starts.length; k++) {
    if (starts[k] >= end) break;
    const i = order[k];
    if (meets(s[i], e[i], start, end)) {
      if (n < out.length) out[n] = i;
      n++;
    }
  }
  for (let k = 0; k < ix.long.length; k++) {
    const i = ix.long[k];
    if (meets(s[i], e[i], start, end)) {
      if (n < out.length) out[n] = i;
      n++;
    }
  }
  return n;
}

/** The first index of `sorted` (by `at`) whose value is at least `v`. */
export function lowerBound(count: number, at: (k: number) => number, v: number): number {
  let lo = 0;
  let hi = count;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (at(mid) < v) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

const meets = (s: number, e: number, start: number, end: number) => s < end && Math.max(e, s + 1) > start;
