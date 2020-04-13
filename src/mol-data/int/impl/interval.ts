/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Tuple from '../tuple';

export const Empty = Tuple.Zero;
export function ofRange(min: number, max: number) { return max < min ? Tuple.create(min, min) : Tuple.create(min, max + 1); }
export function ofBounds(start: number, end: number) { return end <= start ? Tuple.create(start, start) : Tuple.create(start, end); }
export function ofLength(length: number) { return length < 0 ? Tuple.create(0, 0) : Tuple.create(0, length); }
export const is = Tuple.is;

export const start = Tuple.fst;
export const end = Tuple.snd;
export const min = Tuple.fst;
export function max(i: Tuple) { return Tuple.snd(i) - 1; }
export function size(i: Tuple) { return Tuple.snd(i) - Tuple.fst(i); }
export const hashCode = Tuple.hashCode;
export const toString = Tuple.toString;

export function has(int: Tuple, v: number) { return Tuple.fst(int) <= v && v < Tuple.snd(int); }
/** Returns the index of `x` in `set` or -1 if not found. */
export function indexOf(int: Tuple, x: number) { const m = start(int); return x >= m && x < end(int) ? x - m : -1; }
export function getAt(int: Tuple, i: number) { return Tuple.fst(int) + i; }

export const areEqual = Tuple.areEqual;
export function areIntersecting(a: Tuple, b: Tuple) {
    const sa = size(a), sb = size(b);
    if (sa === 0 && sb === 0) return true;
    return sa > 0 && sb > 0 && max(a) >= min(b) && min(a) <= max(b);
}
export function isSubInterval(a: Tuple, b: Tuple) {
    if (!size(a)) return size(b) === 0;
    if (!size(b)) return true;
    return start(a) <= start(b) && end(a) >= end(b);
}

export function findPredecessorIndex(int: Tuple, v: number) {
    const s = start(int);
    if (v <= s) return 0;
    const e = end(int);
    if (v >= e) return e - s;
    return v - s;
}

export function findPredecessorIndexInInterval(int: Tuple, v: number, bounds: Tuple) {
    const bS = start(bounds);
    const s = start(int);
    if (v <= bS + s) return bS;
    const bE = end(bounds);
    if (v >= bE + s) return bE;
    return v - s;
}

export function findRange(int: Tuple, min: number, max: number) {
    return ofBounds(findPredecessorIndex(int, min), findPredecessorIndex(int, max + 1));
}

export function intersect(a: Tuple, b: Tuple) {
    if (!areIntersecting(a, b)) return Empty;
    return ofBounds(Math.max(start(a), start(b)), Math.min(end(a), end(b)));
}

export function intersectionSize(a: Tuple, b: Tuple) {
    return size(findRange(a, min(b), max(b)));
}