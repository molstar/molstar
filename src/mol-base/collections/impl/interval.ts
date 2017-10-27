/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import IntTuple from '../int-tuple'

export const Empty = IntTuple.Zero;
export function ofRange(min: number, max: number) { return max < min ? IntTuple.create(min, min) : IntTuple.create(min, max + 1); }
export function ofBounds(min: number, max: number) { return max <= min ? IntTuple.create(min, min) : IntTuple.create(min, max); }
export const is = IntTuple.is;

export const start = IntTuple.fst;
export const end = IntTuple.snd;
export const min = IntTuple.fst;
export function max(i: IntTuple) { return IntTuple.snd(i) - 1; }
export function size(i: IntTuple) { return IntTuple.snd(i) - IntTuple.fst(i); }
export const hashCode = IntTuple.hashCode;

export function has(int: IntTuple, v: number) { return IntTuple.fst(int) <= v && v < IntTuple.snd(int); }
export function indexOf(int: IntTuple, x: number) { const m = start(int); return x >= m && x < end(int) ? x - m : -1; }
export function getAt(int: IntTuple, i: number) { return IntTuple.fst(int) + i; }

export const areEqual = IntTuple.areEqual;
export function areIntersecting(a: IntTuple, b: IntTuple) {
    const sa = size(a), sb = size(b);
    if (sa === 0 && sb === 0) return true;
    return sa > 0 && sb > 0 && max(a) >= min(b) && min(a) <= max(b);
}
export function isSubInterval(a: IntTuple, b: IntTuple) {
    if (!size(a)) return size(b) === 0;
    if (!size(b)) return true;
    return start(a) <= start(b) && end(a) >= end(b);
}

export function findPredecessorIndex(int: IntTuple, v: number) {
    const s = start(int);
    if (v <= s) return 0;
    const e = end(int);
    if (v >= e) return e - s;
    return v - s;
}

export function findPredecessorIndexInInterval(int: IntTuple, x: number, bounds: IntTuple) {
    const ret = findPredecessorIndex(int, x);
    const s = start(bounds), e = end(bounds);
    return ret <= s ? s : ret >= e ? e : ret;
}

export function findRange(int: IntTuple, min: number, max: number) {
    return ofBounds(findPredecessorIndex(int, min), findPredecessorIndex(int, max + 1));
}

export function intersect(a: IntTuple, b: IntTuple) {
    if (!areIntersecting(a, b)) return Empty;
    return ofBounds(Math.max(start(a), start(b)), Math.min(end(a), end(b)));
}