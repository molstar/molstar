/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Ordered sets of integers: `Interval`, `SortedArray`, `OrderedSet`, `SortedRanges`.
 *
 * A subset port of Mol*'s `mol-data/int` (Mol*'s src/mol-data/int: tuple.ts, impl/interval.ts,
 * impl/sorted-array.ts, impl/ordered-set.ts), Copyright (c) 2017-2026 mol* contributors, MIT
 * (Mol*'s LICENSE); authors David Sehnal, Alexander Rose, Adam Midlik.
 *
 * Changes from Mol*: half-open everywhere (`ofRange(min, stop)`; Mol*'s `ofRange` is inclusive, its
 * `ofBounds` is ours), `lowerBound` for `findPredecessorIndex`, sorted arrays are always `Int32Array`,
 * no hashing, no `subtract`/`isSubset`/segments, `SortedRanges` as half-open pairs (Mol*'s are
 * inclusive) with a binary-search `has`.
 */

import type { Brand } from "./common.js";

/** min and max+1 packed in one double (two 32-bit ints). Compare with `Interval.areEqual`, never `===` (NaN). */
export type Interval<T extends number = number> = Brand<number, "interval"> & { readonly "@elem"?: T };
/** Strictly increasing. */
export type SortedArray<T extends number = number> = Int32Array & { readonly "@elem"?: T };
export type OrderedSet<T extends number = number> = Interval<T> | SortedArray<T>;
/** Sorted, non-overlapping half-open ranges as pairs: [s0, t0, s1, t1, ...]. */
export type SortedRanges<T extends number = number> = Int32Array & { readonly "@elem"?: T };

// IntTuple: two int32 in the bits of one float64, through a shared scratch buffer (allocation-free).
const buf0 = new ArrayBuffer(8);
const buf1 = new ArrayBuffer(8);
const i32 = new Int32Array(buf0);
const f64 = new Float64Array(buf0);
const i32b = new Int32Array(buf1);
const f64b = new Float64Array(buf1);

function pack(a: number, b: number): Interval {
  i32[0] = a;
  i32[1] = b;
  return f64[0] as Interval;
}
function fst(t: number) {
  f64[0] = t;
  return i32[0];
}
function snd(t: number) {
  f64[0] = t;
  return i32[1];
}

// ---------------------------------------------------------------- Interval

// Empty is one value, however an empty range was made, so equality is transitive.
const EMPTY = pack(0, 0);
function iOfRange(min: number, stop: number): Interval {
  return stop <= min ? EMPTY : pack(min, stop);
}
function iSize(t: number) {
  f64[0] = t;
  return i32[1] - i32[0];
}
function iMax(t: number) {
  return snd(t) - 1;
}
function iHas(t: number, x: number) {
  f64[0] = t;
  return i32[0] <= x && x < i32[1];
}
function iIndexOf(t: number, x: number) {
  f64[0] = t;
  return x >= i32[0] && x < i32[1] ? x - i32[0] : -1;
}
function iLowerBound(t: number, x: number) {
  f64[0] = t;
  const s = i32[0];
  const e = i32[1];
  return x <= s ? 0 : x >= e ? e - s : x - s;
}
function iAreEqual(a: number, b: number) {
  f64[0] = a;
  f64b[0] = b;
  return i32[0] === i32b[0] && i32[1] === i32b[1];
}

export const Interval = {
  Empty: EMPTY,
  ofRange: iOfRange as <T extends number>(min: T, stop: T) => Interval<T>,
  is: (x: unknown): x is Interval => typeof x === "number",
  start: fst,
  stop: snd,
  min: fst,
  max: iMax,
  size: iSize,
  has: iHas,
  indexOf: iIndexOf,
  getAt: (t: number, i: number) => fst(t) + i,
  lowerBound: iLowerBound,
  areEqual: iAreEqual,
};

// ---------------------------------------------------------------- SortedArray

const EmptyArray = new Int32Array(0);

function binarySearchRange(xs: Int32Array, value: number, start: number, end: number) {
  let min = start;
  let max = end - 1;
  while (min <= max) {
    if (min + 11 > max) {
      for (let i = min; i <= max; i++) if (value === xs[i]) return i;
      return -1;
    }
    const mid = (min + max) >> 1;
    const v = xs[mid];
    if (value < v) max = mid - 1;
    else if (value > v) min = mid + 1;
    else return mid;
  }
  return -1;
}

/** The first index in [start, end) whose value is >= query; `end` when none. */
function lowerBoundRange(xs: Int32Array, query: number, start: number, end: number): number {
  if (start === end || xs[start] >= query) return start;
  if (xs[end - 1] < query) return end;
  let min = start;
  let max = end;
  while (max - min > 4) {
    const mid = (min + max) >> 1;
    if (xs[mid] >= query) max = mid;
    else min = mid + 1;
  }
  for (let i = min; i < max; i++) if (xs[i] >= query) return i;
  return max;
}

function sIndexOf(xs: Int32Array, v: number) {
  const l = xs.length;
  return l === 0 ? -1 : xs[0] <= v && v <= xs[l - 1] ? binarySearchRange(xs, v, 0, l) : -1;
}
function sLowerBound(xs: Int32Array, x: number) {
  return lowerBoundRange(xs, x, 0, xs.length);
}
function sIsRange(xs: Int32Array) {
  return xs[xs.length - 1] - xs[0] + 1 === xs.length;
}
function sAreEqual(a: Int32Array, b: Int32Array) {
  if (a === b) return true;
  const n = a.length;
  if (n !== b.length || a[0] !== b[0] || a[n - 1] !== b[n - 1]) return false;
  if (sIsRange(a)) return true;
  for (let i = 1; i < n - 1; i++) if (a[i] !== b[i]) return false;
  return true;
}

const range = { startI: 0, startJ: 0, endI: 0, endJ: 0 };
/** For small sets the whole range; for large ones a few binary searches. */
function suitableIntersectionRange(a: Int32Array, b: Int32Array) {
  const la = a.length;
  const lb = b.length;
  const ratio = la / lb;
  if (la >= 128 || lb >= 128 || ratio <= 0.34 || ratio >= 2.99) {
    range.startI = sLowerBound(a, b[0]);
    range.startJ = sLowerBound(b, a[0]);
    range.endI = sLowerBound(a, b[lb - 1] + 1);
    range.endJ = sLowerBound(b, a[la - 1] + 1);
  } else {
    range.startI = 0;
    range.startJ = 0;
    range.endI = la;
    range.endJ = lb;
  }
  return range;
}

function commonCount(a: Int32Array, b: Int32Array, startI: number, startJ: number, endI: number, endJ: number) {
  let i = startI;
  let j = startJ;
  let count = 0;
  while (i < endI && j < endJ) {
    const x = a[i];
    const y = b[j];
    if (x < y) i++;
    else if (x > y) j++;
    else {
      i++;
      j++;
      count++;
    }
  }
  return count;
}

function sUnion(a: Int32Array, b: Int32Array): Int32Array {
  if (a === b) return a;
  const lenA = a.length;
  const lenB = b.length;
  if (lenA === 0) return b;
  if (lenB === 0) return a;
  if (a[0] > b[0]) return sUnion(b, a);
  const { startI, startJ, endI, endJ } = suitableIntersectionRange(a, b);
  const common = commonCount(a, b, startI, startJ, endI, endJ);
  if (common === lenB) return a;
  if (common === lenA) return b;
  const out = new Int32Array(lenA + lenB - common);
  let i = 0;
  let j = 0;
  let o = 0;
  for (i = 0; i < startI; i++) out[o++] = a[i];
  while (j < endJ && a[startI] > b[j]) out[o++] = b[j++];
  while (i < endI && j < endJ) {
    const x = a[i];
    const y = b[j];
    if (x < y) {
      out[o++] = x;
      i++;
    } else if (x > y) {
      out[o++] = y;
      j++;
    } else {
      out[o++] = x;
      i++;
      j++;
    }
  }
  for (; i < endI; i++) out[o++] = a[i];
  for (; j < endJ; j++) out[o++] = b[j];
  for (; i < lenA; i++) out[o++] = a[i];
  for (; j < lenB; j++) out[o++] = b[j];
  return out;
}

function sIntersect(a: Int32Array, b: Int32Array): Int32Array {
  if (a === b) return a;
  if (a.length === 0 || b.length === 0) return EmptyArray;
  const { startI, startJ, endI, endJ } = suitableIntersectionRange(a, b);
  const common = commonCount(a, b, startI, startJ, endI, endJ);
  if (!common) return EmptyArray;
  if (common === b.length) return b;
  if (common === a.length) return a;
  const out = new Int32Array(common);
  let o = 0;
  let i = startI;
  let j = startJ;
  while (i < endI && j < endJ) {
    const x = a[i];
    const y = b[j];
    if (x < y) i++;
    else if (x > y) j++;
    else {
      out[o++] = x;
      i++;
      j++;
    }
  }
  return out;
}

/** Sorts and deduplicates in place (edges); returns a strictly increasing array. */
function sOfUnsorted(xs: Int32Array): SortedArray {
  xs.sort();
  let n = 0;
  for (let i = 0; i < xs.length; i++) if (i === 0 || xs[i] !== xs[n - 1]) xs[n++] = xs[i];
  return n === xs.length ? xs : xs.slice(0, n);
}

export const SortedArray = {
  Empty: EmptyArray as SortedArray<never>,
  ofUnsortedArray: sOfUnsorted as <T extends number>(xs: Int32Array) => SortedArray<T>,
  indexOf: sIndexOf,
  has: (xs: Int32Array, v: number) => sIndexOf(xs, v) >= 0,
  lowerBound: sLowerBound,
  areEqual: sAreEqual,
  union: sUnion as <T extends number>(a: SortedArray<T>, b: SortedArray<T>) => SortedArray<T>,
  intersect: sIntersect as <T extends number>(a: SortedArray<T>, b: SortedArray<T>) => SortedArray<T>,
};

// ---------------------------------------------------------------- OrderedSet

type OS = number | Int32Array;

function ofSortedArray(xs: Int32Array): OS {
  if (!xs.length) return Interval.Empty;
  if (sIsRange(xs)) return iOfRange(xs[0], xs[xs.length - 1] + 1);
  return xs;
}
function size(s: OS) {
  return typeof s === "number" ? iSize(s) : s.length;
}
function min(s: OS) {
  return typeof s === "number" ? fst(s) : s[0];
}
function max(s: OS) {
  return typeof s === "number" ? iMax(s) : s[s.length - 1];
}

function areEqual(a: OS, b: OS) {
  if (typeof a === "number") {
    if (typeof b === "number") return iAreEqual(a, b);
    return areEqualIS(a, b);
  }
  if (typeof b === "number") return areEqualIS(b, a);
  return sAreEqual(a, b);
}
function areEqualIS(a: number, b: Int32Array) {
  const n = iSize(a);
  return n === b.length && (n === 0 || (fst(a) === b[0] && snd(a) === b[n - 1] + 1));
}

function areRangesIntersecting(a: OS, b: OS) {
  const sa = size(a);
  const sb = size(b);
  if (sa === 0 && sb === 0) return true;
  return sa > 0 && sb > 0 && max(a) >= min(b) && min(a) <= max(b);
}
function isRangeSubset(a: OS, b: OS) {
  if (!size(a)) return size(b) === 0;
  if (!size(b)) return true;
  return min(a) <= min(b) && max(a) >= max(b);
}

function union(a: OS, b: OS): OS {
  if (typeof a === "number") {
    if (typeof b === "number") return unionII(a, b);
    return unionSI(b, a);
  }
  if (typeof b === "number") return unionSI(a, b);
  return ofSortedArray(sUnion(a, b));
}
function unionII(a: number, b: number): OS {
  if (iAreEqual(a, b)) return a;
  const sizeA = iSize(a);
  const sizeB = iSize(b);
  if (!sizeB) return a;
  if (!sizeA) return b;
  const minA = fst(a);
  const minB = fst(b);
  if (areRangesIntersecting(a, b) || snd(a) === minB || snd(b) === minA) {
    return iOfRange(Math.min(minA, minB), Math.max(snd(a), snd(b)));
  }
  const [lMin, lSize, rMin, rSize] = minA < minB ? [minA, sizeA, minB, sizeB] : [minB, sizeB, minA, sizeA];
  const out = new Int32Array(sizeA + sizeB);
  for (let i = 0; i < lSize; i++) out[i] = i + lMin;
  for (let i = 0; i < rSize; i++) out[i + lSize] = i + rMin;
  return out;
}
function unionSI(a: Int32Array, b: number): OS {
  const bSize = iSize(b);
  if (!bSize) return ofSortedArray(a);
  if (isRangeSubset(b, a)) return b;
  const lo = fst(b);
  const hi = snd(b);
  const start = sLowerBound(a, lo);
  const end = sLowerBound(a, hi);
  const out = new Int32Array(start + (a.length - end) + bSize);
  let o = 0;
  for (let i = 0; i < start; i++) out[o++] = a[i];
  for (let i = lo; i < hi; i++) out[o++] = i;
  for (let i = end; i < a.length; i++) out[o++] = a[i];
  return ofSortedArray(out);
}

function intersect(a: OS, b: OS): OS {
  if (typeof a === "number") {
    if (typeof b === "number") {
      if (!areRangesIntersecting(a, b) || !iSize(a)) return Interval.Empty;
      return iOfRange(Math.max(fst(a), fst(b)), Math.min(snd(a), snd(b)));
    }
    return intersectSI(b, a);
  }
  if (typeof b === "number") return intersectSI(a, b);
  return ofSortedArray(sIntersect(a, b));
}
function intersectSI(a: Int32Array, b: number): OS {
  if (!iSize(b)) return Interval.Empty;
  const start = sLowerBound(a, fst(b));
  const end = sLowerBound(a, snd(b));
  if (end <= start) return Interval.Empty;
  return ofSortedArray(end - start === a.length ? a : a.subarray(start, end));
}

function forEach<T extends number, C>(s: OrderedSet<T>, f: (x: T, i: number, ctx: C) => void, ctx: C) {
  if (typeof s === "number") {
    const start = fst(s);
    const end = snd(s);
    for (let i = start; i < end; i++) f(i as T, i - start, ctx);
  } else {
    for (let i = 0, n = s.length; i < n; i++) f(s[i] as T, i, ctx);
  }
}

export const OrderedSet = {
  Empty: Interval.Empty as OrderedSet<never>,
  ofRange: iOfRange as <T extends number>(min: T, stop: T) => OrderedSet<T>,
  ofSingleton: <T extends number>(x: T) => iOfRange(x, x + 1) as OrderedSet<T>,
  /** An Interval when contiguous. The caller ensures `xs` is strictly increasing. */
  ofSortedArray: ofSortedArray as <T extends number>(xs: Int32Array) => OrderedSet<T>,
  size,
  isEmpty: (s: OS) => size(s) === 0,
  min,
  max,
  has: (s: OS, x: number) => (typeof s === "number" ? iHas(s, x) : sIndexOf(s, x) >= 0),
  getAt: <T extends number>(s: OrderedSet<T>, i: number) => (typeof s === "number" ? fst(s) + i : s[i]) as T,
  indexOf: (s: OS, x: number) => (typeof s === "number" ? iIndexOf(s, x) : sIndexOf(s, x)),
  /** The first i with getAt(s, i) >= x; size(s) when none. */
  lowerBound: (s: OS, x: number) => (typeof s === "number" ? iLowerBound(s, x) : sLowerBound(s, x)),
  forEach,
  areEqual,
  /** Allocates: edges. */
  union: union as <T extends number>(a: OrderedSet<T>, b: OrderedSet<T>) => OrderedSet<T>,
  /** Allocates: edges. */
  intersect: intersect as <T extends number>(a: OrderedSet<T>, b: OrderedSet<T>) => OrderedSet<T>,
};

// ---------------------------------------------------------------- SortedRanges

export const SortedRanges = {
  count: (r: SortedRanges) => r.length >> 1,
  size(r: SortedRanges) {
    let n = 0;
    for (let i = 0; i < r.length; i += 2) n += r[i + 1] - r[i];
    return n;
  },
  has(r: SortedRanges, x: number) {
    let lo = 0;
    let hi = r.length >> 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (r[2 * mid] <= x) lo = mid + 1;
      else hi = mid;
    }
    return lo > 0 && x < r[2 * lo - 1];
  },
};
