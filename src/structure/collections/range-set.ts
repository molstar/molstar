/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from './iterator'

interface RangeSet {
    readonly size: number,
    has(x: number): boolean,
    indexOf(x: number): number,
    elementAt(i: number): number,
    elements(): Iterator<number>
}

namespace RangeSet {
    interface Impl extends RangeSet {
        readonly min: number,
        readonly max: number,
        toArray(): ArrayLike<number>
    }

    export function hashCode(a: RangeSet) {
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

    export function areEqual(a: RangeSet, b: RangeSet) {
        if (a === b) return true;
        if (a instanceof RangeImpl) {
            if (b instanceof RangeImpl) return a.min === b.min && a.max === b.max;
            return equalAR(b as ArrayImpl, a);
        } else if (b instanceof RangeImpl) {
            return equalAR(a as ArrayImpl, b);
        }
        return equalAA(a as ArrayImpl, b as ArrayImpl);
    }

    export function union(a: RangeSet, b: RangeSet) {
        if (a instanceof RangeImpl) {
            if (b instanceof RangeImpl) return unionRR(a, b);
            return unionAR(b as ArrayImpl, a);
        } else if (b instanceof RangeImpl) {
            return unionAR(a as ArrayImpl, b);
        } else return unionAA((a as Impl).toArray(), (b as Impl).toArray());
    }

    export function intersect(a: RangeSet, b: RangeSet) {
        if (a instanceof RangeImpl) {
            if (b instanceof RangeImpl) return intersectRR(a, b);
            return intersectAR(b as ArrayImpl, a);
        } else if (b instanceof RangeImpl) {
            return intersectAR(a as ArrayImpl, b);
        } else {
            const ai = a as Impl, bi = b as Impl;
            if (!areRangesIntersecting(ai, bi)) return Empty;
            return intersectAA(ai.toArray(), bi.toArray());
        }
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
        toArray() { return this.values; }
        elementAt(i: number) { return this.values[i]; }
        elements() { return Iterator.Array(this.values); }

        constructor(public values: ArrayLike<number>) {
            this.min = values[0];
            this.max = values[values.length - 1];
            this.size = values.length;
        }
    }

    export function ofSingleton(value: number): RangeSet { return new RangeImpl(value, value); }
    export function ofRange(min: number, max: number): RangeSet { return new RangeImpl(min, max); }
    /** It is the responsibility of the caller to ensure the array is sorted and contains unique values. */
    export function ofSortedArray(xs: ArrayLike<number>): RangeSet {
        if (!xs.length) return Empty;
        // check if the array is just a range
        if (xs[xs.length - 1] - xs[0] + 1 === xs.length) return ofRange(xs[0], xs[xs.length - 1]);
        return new ArrayImpl(xs);
    }
    export const Empty = ofRange(0, -1);

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

    function areRangesIntersecting(a: Impl, b: Impl) {
        return a.size > 0 && b.size > 0 && a.max >= b.min && a.min <= b.max;
    }

    function unionRR(a: RangeImpl, b: RangeImpl) {
        if (!a.size) return b;
        if (!b.size) return a;
        if (areRangesIntersecting(a, b)) return ofRange(Math.min(a.min, b.min), Math.max(a.max, b.max));
        let l, r;
        if (a.min < b.min) { l = a; r = b; }
        else { l = b; r = a; }
        const arr = new Int32Array(a.size + b.size);
        for (let i = 0; i < l.size; i++) arr[i] = i + l.min;
        for (let i = 0; i < r.size; i++) arr[i + l.size] = i + r.min;
        return ofSortedArray(arr);
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

        return ofSortedArray(indices);
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

        return ofSortedArray(indices);
    }

    function intersectRR(a: RangeImpl, b: RangeImpl) {
        if (!areRangesIntersecting(a, b)) return Empty;
        return ofRange(Math.max(a.min, b.min), Math.min(a.max, b.max));
    }

    function intersectAR(a: ArrayImpl, r: RangeImpl) {
        const xs = a.values;
        let resultSize = 0;
        for (let i = 0, _i = xs.length; i < _i; i++) {
            if (r.has(xs[i])) resultSize++;
        }

        if (!resultSize) return Empty;

        const indices = new Int32Array(resultSize);
        let offset = 0;

        for (let i = 0, _i = xs.length; i < _i; i++) {
            if (r.has(xs[i])) indices[offset++] = xs[i];
        }

        return ofSortedArray(indices);
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

        if (!resultSize) return Empty;

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

        return ofSortedArray(indices);
    }
}

export default RangeSet