/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from './iterator'
import { hash3, hash4 } from './hash-functions'

/** An immutable ordered set. */
interface OrderedSet {
    readonly size: number,
    has(x: number): boolean,
    indexOf(x: number): number,
    elementAt(i: number): number,
    elements(): Iterator<number>
}

interface Impl extends OrderedSet {
    readonly min: number,
    readonly max: number
}

class RangeImpl implements Impl {
    size: number;
    has(x: number) { return x >= this.min && x <= this.max; }
    indexOf(x: number) { return x >= this.min && x <= this.max ? x - this.min : -1; }
    elementAt(i: number) { return this.min + i; }
    elements() { return Iterator.Range(this.min, this.max); }

    constructor(public min: number, public max: number) {
        this.size = max - min + 1;
    }
}

class ArrayImpl implements Impl {
    size: number;
    public min: number;
    public max: number;
    has(x: number) { return x >= this.min && x <= this.max && binarySearch(this.values, x) >= 0; }
    indexOf(x: number) { return x >= this.min && x <= this.max ? binarySearch(this.values, x) : -1; }
    elementAt(i: number) { return this.values[i]; }
    elements() { return Iterator.Array(this.values); }

    constructor(public values: ArrayLike<number>) {
        this.min = values[0];
        this.max = values[values.length - 1];
        this.size = values.length;
    }
}

namespace OrderedSet {
    export function ofSingleton(value: number): OrderedSet { return new RangeImpl(value, value); }
    export function ofRange(min: number, max: number): OrderedSet { return max < min ? Empty : new RangeImpl(min, max); }
    /** It is the responsibility of the caller to ensure the array is sorted and contains unique values. */
    export function ofSortedArray(xs: ArrayLike<number>): OrderedSet {
        if (!xs.length) return Empty;
        // check if the array is just a range
        if (xs[xs.length - 1] - xs[0] + 1 === xs.length) return ofRange(xs[0], xs[xs.length - 1]);
        return new ArrayImpl(xs);
    }
    export const Empty: OrderedSet = new RangeImpl(0, -1);

    export function isEmpty(a: OrderedSet) { return a.size === 0; }
    export function min(a: OrderedSet) { return (a as Impl).min; }
    export function max(a: OrderedSet) { return (a as Impl).max; }

    export function hashCode(a: OrderedSet) {
        // hash of tuple (size, min, max, mid)
        const { size } = a;
        if (!size) return 0;
        if (size > 2) return hash4(size, a.elementAt(0), a.elementAt(size - 1), a.elementAt(size >> 0));
        return hash3(size, a.elementAt(0), a.elementAt(size - 1));
    }
    // TODO: possibly add more hash functions to allow for multilevel hashing.

    export function areEqual(a: OrderedSet, b: OrderedSet) {
        if (a === b) return true;
        if (a instanceof RangeImpl) {
            if (b instanceof RangeImpl) return a.min === b.min && a.max === b.max;
            return equalAR(b as ArrayImpl, a);
        } else if (b instanceof RangeImpl) {
            return equalAR(a as ArrayImpl, b);
        }
        return equalAA(a as ArrayImpl, b as ArrayImpl);
    }

    export function areIntersecting(a: OrderedSet, b: OrderedSet) {
        if (a === b) return true;
        if (!areRangesIntersecting(a, b)) return false;
        // if at least one is "range", they must now intersect
        if (a instanceof RangeImpl || b instanceof RangeImpl) return true;
        return areIntersectingAA((a as ArrayImpl).values, (b as ArrayImpl).values);
    }

    /** Check if the 2nd argument is a subset of the 1st */
    export function isSubset(a: OrderedSet, toTest: OrderedSet) {
        if (a === toTest) return true;
        if (!isRangeSubset(a, toTest)) return false;
        if (!toTest.size || a instanceof RangeImpl) return true;
        if (toTest instanceof RangeImpl) return a.indexOf(max(toTest)) - a.indexOf(min(toTest)) + 1 === toTest.size;
        return isSubsetAA((a as ArrayImpl).values, (toTest as ArrayImpl).values);
    }

    export function getIntervalRange(a: OrderedSet, min: number, max: number) {
        const { start, end } = a instanceof RangeImpl
            ? getStartEndR(a, min, max)
            : getStartEndA((a as ArrayImpl).values, min, max);
        return { start, end };
    }

    export function union(a: OrderedSet, b: OrderedSet) {
        if (a === b) return a;
        if (a instanceof RangeImpl) {
            if (b instanceof RangeImpl) return unionRR(a, b);
            return unionAR(b as ArrayImpl, a);
        } else if (b instanceof RangeImpl) {
            return unionAR(a as ArrayImpl, b);
        } else return unionAA(a as ArrayImpl, b as ArrayImpl);
    }

    export function intersect(a: OrderedSet, b: OrderedSet) {
        if (a === b) return a;
        if (a instanceof RangeImpl) {
            if (b instanceof RangeImpl) return intersectRR(a, b);
            return intersectAR(b as ArrayImpl, a);
        } else if (b instanceof RangeImpl) {
            return intersectAR(a as ArrayImpl, b);
        } else {
            if (!areRangesIntersecting(a, b)) return Empty;
            return intersectAA((a as ArrayImpl).values, (b as ArrayImpl).values);
        }
    }

    export function subtract(a: OrderedSet, b: OrderedSet) {
        if (a === b) return Empty;
        if (!areRangesIntersecting(a, b)) return a;

        if (a instanceof RangeImpl) {
            if (b instanceof RangeImpl) return substractRR(a, b);
            return subtractRA(a, (b as ArrayImpl).values);
        } else if (b instanceof RangeImpl) {
            return subtractAR(a as ArrayImpl, b);
        } else {
            return subtractAA(a as ArrayImpl, b as ArrayImpl);
        }
    }
}

function min(a: OrderedSet) { return (a as Impl).min; }
function max(a: OrderedSet) { return (a as Impl).max; }

function binarySearch(xs: ArrayLike<number>, value: number) {
    let min = 0, max = xs.length - 1;
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

function binarySearchIndex(xs: ArrayLike<number>, value: number) {
    let min = 0, max = xs.length - 1;
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

const _maxIntRangeRet = { i: 0, j: 0, endA: 0, endB: 0 };
function getMaxIntersectionRange(xs: ArrayLike<number>, ys: ArrayLike<number>) {
    const la = xs.length - 1, lb = ys.length - 1;
    _maxIntRangeRet.i = binarySearchIndex(xs, ys[0]);
    _maxIntRangeRet.j = binarySearchIndex(ys, xs[0]);
    _maxIntRangeRet.endA = Math.min(binarySearchIndex(xs, ys[lb]), la);
    _maxIntRangeRet.endB = Math.min(binarySearchIndex(ys, xs[la]), lb);
    return _maxIntRangeRet;
}

const _startEndRet = { start: 0, end: 0 };

function getStartEndR(r: RangeImpl, min: number, max: number) {
    if (max < min) {
        _startEndRet.start = 0;
        _startEndRet.end = 0;
        return _startEndRet;
    }
    _startEndRet.start = min <= r.max ? Math.max(r.min, min) - r.min : r.size;
    _startEndRet.end = max >= r.min ? Math.min(r.max, max) - r.min + 1 : r.size;
    return _startEndRet;
}

function getStartEndA(xs: ArrayLike<number>, min: number, max: number) {
    _startEndRet.start = binarySearchIndex(xs, min);
    let end = binarySearchIndex(xs, max);
    if (xs[end] === max) end++;
    _startEndRet.end = end;
    return _startEndRet;
}

function equalAR(a: ArrayImpl, b: RangeImpl) {
    return a.size === b.size && a.min === b.min && a.max === b.max;
}

function equalAA(a: ArrayImpl, b: ArrayImpl) {
    if (a.size !== b.size || a.min !== b.min || a.max !== b.max) return false;
    const { size, values: xs } = a;
    const { values: ys } = b;
    for (let i = 0; i < size; i++) {
        if (xs[i] !== ys[i]) return false;
    }
    return true;
}

function areIntersectingAA(xs: ArrayLike<number>, ys: ArrayLike<number>) {
    let { i, j, endA, endB } = getMaxIntersectionRange(xs, ys);
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else return true;
    }
    return false;
}

function isSubsetAA(xs: ArrayLike<number>, ys: ArrayLike<number>) {
    const lenB = ys.length;
    let { i, j, endA, endB } = getMaxIntersectionRange(xs, ys);
    // the 2nd array must be able to advance by at least lenB elements
    if (endB - j + 1 < lenB || endA - j + 1 < lenB) return false;

    let equal = 0;
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; equal++; }
    }
    return equal === lenB;
}

function areRangesIntersecting(a: OrderedSet, b: OrderedSet) {
    return a.size > 0 && b.size > 0 && max(a) >= min(b) && min(a) <= max(b);
}

function isRangeSubset(a: OrderedSet, b: OrderedSet) {
    if (!a.size) return b.size === 0;
    if (!b.size) return true;
    return min(a) <= min(b) && max(a) >= max(b);
}

function unionRR(a: RangeImpl, b: RangeImpl) {
    if (!a.size) return b;
    if (!b.size) return a;
    if (areRangesIntersecting(a, b)) return OrderedSet.ofRange(Math.min(a.min, b.min), Math.max(a.max, b.max));
    let l, r;
    if (a.min < b.min) { l = a; r = b; }
    else { l = b; r = a; }
    const arr = new Int32Array(a.size + b.size);
    for (let i = 0; i < l.size; i++) arr[i] = i + l.min;
    for (let i = 0; i < r.size; i++) arr[i + l.size] = i + r.min;
    return OrderedSet.ofSortedArray(arr);
}

function unionAR(a: ArrayImpl, b: RangeImpl) {
    if (!b.size) return a;
    // is the array fully contained in the range?
    if (a.min >= b.min && a.max <= b.max) return b;

    const xs = a.values;
    const { min, max } = b;
    const { start, end } = getStartEndA(xs, min, max);

    const size = start + (xs.length - end) + b.size;
    const indices = new Int32Array(size);
    let offset = 0;
    for (let i = 0; i < start; i++) indices[offset++] = xs[i];
    for (let i = min; i <= max; i++) indices[offset++] = i;
    for (let i = end, _i = xs.length; i < _i; i++) indices[offset] = xs[i];

    return OrderedSet.ofSortedArray(indices);
}

function unionAA(a: ArrayImpl, b: ArrayImpl) {
    const xs = a.values, ys = b.values;
    const lenA = xs.length, lenB = ys.length;

    let { i: sI, j: sJ, endA, endB } = getMaxIntersectionRange(xs, ys);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; commonCount++; }
    }

    if (!commonCount) return a;
    if (commonCount >= lenA) return OrderedSet.Empty

    const resultSize = lenA + lenB - commonCount;
    const l = Math.min(min(a), min(b)), r = Math.max(max(a), max(b));
    // is this just a range?
    if (resultSize === r - l + 1) {
        return OrderedSet.ofRange(l, r);
    }

    const indices = new Int32Array(lenA + lenB - commonCount);
    let offset = 0;

    // insert the "prefixes"
    for (let k = 0; k < sI; k++) indices[offset++] = xs[k];
    for (let k = 0; k < sJ; k++) indices[offset++] = xs[k];

    // insert the common part
    i = sI;
    j = sJ;
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { indices[offset++] = x; i++; }
        else if (x > y) { indices[offset++] = y; j++; }
        else { indices[offset++] = x; i++; j++; }
    }

    // insert the "tail"
    for (; i < lenA; i++) indices[offset++] = xs[i];
    for (; j < lenB; j++) indices[offset++] = ys[j];

    return OrderedSet.ofSortedArray(indices);
}

function intersectRR(a: RangeImpl, b: RangeImpl) {
    if (!areRangesIntersecting(a, b)) return OrderedSet.Empty;
    return OrderedSet.ofRange(Math.max(a.min, b.min), Math.min(a.max, b.max));
}

function intersectAR(a: ArrayImpl, r: RangeImpl) {
    if (!r.size) return OrderedSet.Empty;

    const xs = a.values;
    const { start, end } = getStartEndA(xs, r.min, r.max);
    const resultSize = end - start;
    if (!resultSize) return OrderedSet.Empty;

    const indices = new Int32Array(resultSize);
    let offset = 0;
    for (let i = start; i < end; i++) {
        indices[offset++] = xs[i];
    }
    return OrderedSet.ofSortedArray(indices);
}

function intersectAA(xs: ArrayLike<number>, ys: ArrayLike<number>) {
    let { i: sI, j: sJ, endA, endB } = getMaxIntersectionRange(xs, ys);
    let i = sI, j = sJ;
    let resultSize = 0;
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; resultSize++; }
    }

    if (!resultSize) return OrderedSet.Empty;

    const indices = new Int32Array(resultSize);
    let offset = 0;
    i = sI;
    j = sJ;
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { indices[offset++] = x; i++; j++; }
    }

    return OrderedSet.ofSortedArray(indices);
}

function substractRR(a: RangeImpl, b: RangeImpl) {
    if (a.size === 0 || b.size === 0) return a;
    // is A subset of B? ==> Empty
    if (isRangeSubset(b, a)) return OrderedSet.Empty;
    if (isRangeSubset(a, b)) {
        // this splits the interval into two, gotta represent it as a set.
        const l = b.min - a.min, r = a.max - b.max;
        const ret = new Int32Array(l + r);
        let offset = 0;
        for (let i = 0; i < l; i++) ret[offset++] = a.min + i;
        for (let i = 1; i <= r; i++) ret[offset++] = b.max + i;
        return OrderedSet.ofSortedArray(ret);
    }
    // non intersecting ranges are handled by top-level substract.
    // at this point, b either contains a.min or a.max, but not both.
    if (a.min < b.min) return OrderedSet.ofRange(a.min, b.min - 1);
    return OrderedSet.ofRange(b.max + 1, a.max);
}

function subtractAR(a: ArrayImpl, r: RangeImpl) {
    if (!r.size) return a;

    const xs = a.values;
    const { min, max } = r;
    const { start, end } = getStartEndA(xs, min, max);
    const size = xs.length - (end - start);
    if (size <= 0) return OrderedSet.Empty;
    const ret = new Int32Array(size);
    let offset = 0;
    for (let i = 0; i < start; i++) ret[offset++] = xs[i];
    for (let i = end, _i = xs.length; i < _i; i++) ret[offset++] = xs[i];
    return OrderedSet.ofSortedArray(ret);
}

function subtractRA(r: RangeImpl, ys: ArrayLike<number>) {
    if (!r.size) return r;

    const { min, max } = r;
    const { start, end } = getStartEndA(ys, min, max);
    const commonCount = end - start;
    const resultSize = r.size - commonCount;
    if (resultSize <= 0) return OrderedSet.Empty;
    const ret = new Int32Array(resultSize);
    const li = ys.length - 1;
    const fst = ys[Math.min(start, li)], last = ys[Math.min(end, li)];
    let offset = 0;
    for (let i = min; i < fst; i++) ret[offset++] = i;
    for (let i = fst; i <= last; i++) {
        if (binarySearch(ys, i) < 0) ret[offset++] = i;
    }
    for (let i = last + 1; i <= max; i++) ret[offset++] = i;
    return OrderedSet.ofSortedArray(ret);
}

function subtractAA(a: ArrayImpl, b: ArrayImpl) {
    const xs = a.values, ys = b.values;
    const lenA = xs.length;

    let { i: sI, j: sJ, endA, endB } = getMaxIntersectionRange(xs, ys);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; commonCount++; }
    }

    if (!commonCount) return a;
    if (commonCount >= lenA) return OrderedSet.Empty;

    const indices = new Int32Array(lenA - commonCount);
    let offset = 0;

    // insert the "prefix"
    for (let k = 0; k < sI; k++) indices[offset++] = xs[k];

    i = sI;
    j = sJ;
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { indices[offset++] = x; i++; }
        else if (x > y) { j++; }
        else { i++; j++; }
    }

    // insert the "tail"
    for (; i < lenA; i++) indices[offset++] = xs[i];

    return OrderedSet.ofSortedArray(indices);
}

export default OrderedSet