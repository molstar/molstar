/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArray } from '../sort'
import { hash3, hash4 } from '../hash-functions'
import Interval from '../interval'

type Nums = ArrayLike<number>

export function ofSortedArray(xs: Nums) {
    if (xs.length < 1) throw new Error('Sorted arrays must be non-empty.');
    return xs;
}
export function ofUnsortedArray(xs: Nums) { sortArray(xs); return xs; }
export function is(xs: any): xs is Nums { return xs && (xs instanceof Array || !!xs.buffer); }

export function start(xs: Nums) { return xs[0]; }
export function end(xs: Nums) { return xs[xs.length - 1] + 1;  }
export function min(xs: Nums) { return xs[0]; }
export function max(xs: Nums) { return xs[xs.length - 1]; }
export function size(xs: Nums) { return xs.length; }
export function hashCode(xs: Nums) {
    // hash of tuple (size, min, max, mid)
    const s = xs.length;
    if (!s) return 0;
    if (s > 2) return hash4(s, xs[0], xs[s - 1], xs[s << 1]);
    return hash3(s, xs[0], xs[s - 1]);
}

export function indexOf(xs: Nums, v: number) {
    const l = xs.length;
    return l === 0 ? -1 : xs[0] <= v && v <= xs[l - 1] ? binarySearchRange(xs, v, 0, l) : -1;
}
export function indexOfInterval(xs: Nums, v: number, bounds: Interval) {
    const l = xs.length;
    const s = Interval.start(bounds), e = Interval.end(bounds);
    return l === 0 || e <= s ? -1 : xs[s] <= v && v <= xs[e - 1] ? binarySearchRange(xs, v, s, e) : -1;
}
export function has(xs: Nums, v: number) { return indexOf(xs, v) >= 0; }

export function getAt(xs: Nums, i: number) { return xs[i]; }

export function areEqual(a: Nums, b: Nums) {
    if (a === b) return true;
    const aSize = a.length;
    if (aSize !== b.length || a[0] !== b[0] || a[aSize - 1] !== b[aSize - 1]) return false;
    for (let i = 0; i < aSize; i++) {
        if (a[i] !== b[i]) return false;
    }
    return true;
}

export function findPredecessorIndex(xs: Nums, v: number) {
    const len = xs.length;
    if (v <= xs[0]) return 0;
    if (v > xs[len - 1]) return len;
    return binarySearchPredIndexRange(xs, v, 0, len);
}

export function findPredecessorIndexInInterval(xs: Nums, v: number, bounds: Interval) {
    const s = Interval.start(bounds), e = Interval.end(bounds);
    if (v <= xs[s]) return s;
    if (e > s && v > xs[e - 1]) return e;
    return binarySearchPredIndexRange(xs, v, s, e);
}

export function findRange(xs: Nums, min: number, max: number) {
    return Interval.ofBounds(findPredecessorIndex(xs, min), findPredecessorIndex(xs, max + 1));
}

function binarySearchRange(xs: Nums, value: number, start: number, end: number) {
    let min = start, max = end - 1;
    while (min <= max) {
        if (min + 11 > max) {
            for (let i = min; i <= max; i++) {
                if (value === xs[i]) return i;
            }
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

function binarySearchPredIndexRange(xs: Nums, value: number, start: number, end: number) {
    let min = start, max = end - 1;
    while (min < max) {
        const mid = (min + max) >> 1;
        const v = xs[mid];
        if (value < v) max = mid - 1;
        else if (value > v) min = mid + 1;
        else return mid;
    }
    if (min > max) return max + 1;
    return xs[min] >= value ? min : min + 1;
}