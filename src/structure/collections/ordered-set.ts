/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from './iterator'

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
    toArray() {
        const ret = new Array(this.size);
        for (let i = 0; i < this.size; i++) ret[i] = i + this.min;
        return ret;
    }
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
    export function isEmpty(a: OrderedSet) { return a.size === 0; }

    export function hashCode(a: OrderedSet) {
        // hash of tuple (size, min, max, mid)
        const { size } = a;
        let hash = 23;
        if (!size) return hash;
        hash = 31 * hash + size;
        hash = 31 * hash + a.elementAt(0);
        hash = 31 * hash + a.elementAt(size - 1);
        if (size > 2) hash = 31 * hash + a.elementAt(size >> 1);
        return hash;
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

    export function union(a: OrderedSet, b: OrderedSet) {
        if (a instanceof RangeImpl) {
            if (b instanceof RangeImpl) return unionRR(a, b);
            return unionAR(b as ArrayImpl, a);
        } else if (b instanceof RangeImpl) {
            return unionAR(a as ArrayImpl, b);
        } else return unionAA((a as ArrayImpl).values, (b as ArrayImpl).values);
    }

    export function intersect(a: OrderedSet, b: OrderedSet) {
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

    export function ofSingleton(value: number): OrderedSet { return new RangeImpl(value, value); }
    export function ofRange(min: number, max: number): OrderedSet { return max < min ? Empty : new RangeImpl(min, max); }
    /** It is the responsibility of the caller to ensure the array is sorted and contains unique values. */
    export function ofSortedArray(xs: ArrayLike<number>): OrderedSet {
        if (!xs.length) return Empty;
        // check if the array is just a range
        if (xs[xs.length - 1] - xs[0] + 1 === xs.length) return ofRange(xs[0], xs[xs.length - 1]);
        return new ArrayImpl(xs);
    }
    export const Empty = new RangeImpl(0, -1);
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
    const la = xs.length, lb = ys.length;
    let i = 0, j = 0;
    while (i < la && j < lb) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else return true;
    }
    return false;
}

function isSubsetAA(xs: ArrayLike<number>, ys: ArrayLike<number>) {
    const la = xs.length, lb = ys.length;
    let i = 0, j = 0, equal = 0;
    while (i < la && j < lb) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; equal++; }
    }
    return equal === lb;
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

    let start = 0, end = xs.length - 1;
    while (xs[start] < min) { start++; }
    while (xs[end] > max) { end--; }
    end++;

    const size = start + (xs.length - end) + b.size;
    const indices = new Int32Array(size);
    let offset = 0;
    for (let i = 0; i < start; i++) indices[offset++] = xs[i];
    for (let i = min; i <= max; i++) indices[offset++] = i;
    for (let i = end, _i = xs.length; i < _i; i++) indices[offset] = xs[i];

    return OrderedSet.ofSortedArray(indices);
}

function unionAA(xs: ArrayLike<number>, ys: ArrayLike<number>) {
    const la = xs.length, lb = ys.length;

    // sorted list merge.

    let i = 0, j = 0, resultSize = 0;
    while (i < la && j < lb) {
        const x = xs[i], y = ys[j];
        resultSize++;
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; }
    }
    resultSize += Math.max(la - i, lb - j);

    const indices = new Int32Array(resultSize);
    let offset = 0;
    i = 0;
    j = 0;
    while (i < la && j < lb) {
        const x = xs[i], y = ys[j];
        if (x < y) { indices[offset++] = x; i++; }
        else if (x > y) { indices[offset++] = y; j++; }
        else { indices[offset++] = x; i++; j++; }
    }
    for (; i < la; i++) { indices[offset++] = xs[i]; }
    for (; j < lb; j++) { indices[offset++] = ys[j]; }

    return OrderedSet.ofSortedArray(indices);
}

function intersectRR(a: RangeImpl, b: RangeImpl) {
    if (!areRangesIntersecting(a, b)) return OrderedSet.Empty;
    return OrderedSet.ofRange(Math.max(a.min, b.min), Math.min(a.max, b.max));
}

function intersectAR(a: ArrayImpl, r: RangeImpl) {
    const xs = a.values;
    let resultSize = 0;
    for (let i = 0, _i = xs.length; i < _i; i++) {
        if (r.has(xs[i])) resultSize++;
    }

    if (!resultSize) return OrderedSet.Empty;

    const indices = new Int32Array(resultSize);
    let offset = 0;

    for (let i = 0, _i = xs.length; i < _i; i++) {
        if (r.has(xs[i])) indices[offset++] = xs[i];
    }

    return OrderedSet.ofSortedArray(indices);
}

function intersectAA(xs: ArrayLike<number>, ys: ArrayLike<number>) {
    const la = xs.length, lb = ys.length;

    // a variation on sorted list merge.

    let i = 0, j = 0, resultSize = 0;
    while (i < la && j < lb) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; resultSize++; }
    }

    if (!resultSize) return OrderedSet.Empty;

    const indices = new Int32Array(resultSize);
    let offset = 0;
    i = 0;
    j = 0;
    while (i < la && j < lb) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { indices[offset++] = x; i++; j++; }
    }

    return OrderedSet.ofSortedArray(indices);
}

export default OrderedSet