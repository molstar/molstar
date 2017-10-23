/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import IntPair from './int-pair'
import { hash3, hash4 } from './hash-functions'

type Nums = ArrayLike<number>
type OrderedSet = number | Nums

/** An immutable ordered set. */
namespace OrderedSet {
    export function ofSingleton(value: number): OrderedSet { return IntPair.pack1(value, value); }
    export function ofRange(min: number, max: number): OrderedSet { return max < min ? Empty : IntPair.pack1(min, max); }
    /** It is the responsibility of the caller to ensure the array is sorted and contains unique values. */
    export function ofSortedArray(xs: Nums): OrderedSet {
        if (!xs.length) return Empty;
        // check if the array is just a range
        if (xs[xs.length - 1] - xs[0] + 1 === xs.length) return ofRange(xs[0], xs[xs.length - 1]);
        return xs;
    }
    export const Empty: OrderedSet = IntPair.pack1(0, -1);

    export function size(set: OrderedSet) { return typeof set === 'number' ? sizeR(set) : set.length; }
    export function has(set: OrderedSet, x: number) { return typeof set === 'number' ? hasR(set, x) : hasA(set, x); }
    export function indexOf(set: OrderedSet, x: number) { return typeof set === 'number' ? indexOfR(set, x) : indexOfA(set, x); }
    export function elementAt(set: OrderedSet, i: number) { return typeof set === 'number' ? elementAtR(set, i) : set[i]; }
    export function min(set: OrderedSet) { return typeof set === 'number' ? minR(set) : set[0]; }
    export function max(set: OrderedSet) { return typeof set === 'number' ? maxR(set) : set[set.length - 1]; }

    export function hashCode(set: OrderedSet) {
        // hash of tuple (size, min, max, mid)
        const s = size(set);
        if (!s) return 0;
        if (s > 2) return hash4(s, elementAt(set, 0), elementAt(set, s - 1), elementAt(set, s >> 1));
        return hash3(s, elementAt(set, 0), elementAt(set, s - 1));
    }
    // TODO: possibly add more hash functions to allow for multilevel hashing.

    export function areEqual(a: OrderedSet, b: OrderedSet) {
        if (typeof a === 'number') {
            if (typeof b === 'number') return equalRR(a, b);
            return false;
        } else if (typeof b === 'number') return false;
        else if (a === b) return true;
        return equalAA(a, b);
    }

    export function areIntersecting(a: OrderedSet, b: OrderedSet) {
        // if at least one is "range", they must now intersect
        if (typeof a === 'number') {
            if (typeof b === 'number') return equalRR(a, b) || areRangesIntersecting(a, b);
            return areRangesIntersecting(a, b);
        }
        if (!areRangesIntersecting(a, b)) return false;
        else if (typeof b === 'number') return false;
        if (a === b) return true;
        return areIntersectingAA(a, b);
    }

    /** Check if the 2nd argument is a subset of the 1st */
    export function isSubset(set: OrderedSet, toTest: OrderedSet) {
        if (set === toTest) return true;
        if (!isRangeSubset(set, toTest)) return false;
        const testSize = size(toTest);
        if (typeof set === 'number' || !testSize) return true;
        if (typeof toTest === 'number') return indexOf(set, maxR(toTest)) - indexOf(set, minR(toTest)) + 1 === testSize;
        return isSubsetAA(set, toTest);
    }

    export function getInsertionIndex(set: OrderedSet, x: number) {
        return typeof set === 'number' ? rangeSearchIndex(set, x) : binarySearchIndex(set, x);
    }

    export function getIntervalRange(set: OrderedSet, min: number, max: number) {
        const { start, end } = getStartEnd(set, min, max);
        return { start, end };
    }

    export function union(a: OrderedSet, b: OrderedSet) {
        if (a === b) return a;
        if (typeof a === 'number') {
            if (typeof b === 'number') return unionRR(a, b);
            return unionAR(b, a);
        } else if (typeof b === 'number') {
            return unionAR(a, b);
        } else return unionAA(a, b);
    }

    export function intersect(a: OrderedSet, b: OrderedSet) {
        if (a === b) return a;
        if (typeof a === 'number') {
            if (typeof b === 'number') return intersectRR(a, b);
            return intersectAR(b, a);
        } else if (typeof b === 'number') {
            return intersectAR(a, b);
        } else {
            if (!areRangesIntersecting(a, b)) return Empty;
            return intersectAA(a, b);
        }
    }

    export function subtract(a: OrderedSet, b: OrderedSet) {
        if (a === b) return Empty;
        if (!areRangesIntersecting(a, b)) return a;

        if (typeof a === 'number') {
            if (typeof b === 'number') return substractRR(a, b);
            return subtractRA(a, b);
        } else if (typeof b === 'number') {
            return subtractAR(a, b);
        } else {
            return subtractAA(a, b);
        }
    }
}

import S = OrderedSet

const minR = IntPair.fst
const maxR = IntPair.snd
const equalRR = IntPair.areEqual

const _eR = IntPair.zero();
function sizeR(set: number) { IntPair.unpack(set, _eR); return _eR.snd - _eR.fst + 1; }
function hasR(set: number, x: number) { IntPair.unpack(set, _eR); return x >= _eR.fst && x <= _eR.snd; }
function indexOfR(set: number, x: number) { IntPair.unpack(set, _eR); return x >= _eR.fst && x <= _eR.snd ? x - _eR.fst : -1; }
function elementAtR(set: number, i: number) { return IntPair.fst(set) + i; }

function hasA(set: Nums, x: number) { return x >= set[0] && x <= set[set.length - 1] && binarySearch(set, x) >= 0; }
function indexOfA(set: Nums, x: number) { return x >= set[0] && x <= set[set.length - 1] ? binarySearch(set, x) : -1; }

function binarySearch(xs: Nums, value: number) {
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

function binarySearchIndex(xs: Nums, value: number) {
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

const _rsiR = IntPair.zero();
function rangeSearchIndex(r: number, value: number) {
    IntPair.unpack(r, _rsiR);
    if (value < _rsiR.fst) return 0;
    if (value > _rsiR.snd) return _rsiR.snd - _rsiR.fst + 1;
    return value - _rsiR.fst;
}

const _maxIntRangeRet = { i: 0, j: 0, endA: 0, endB: 0 };
function getMaxIntersectionRange(xs: Nums, ys: Nums) {
    const la = xs.length - 1, lb = ys.length - 1;
    _maxIntRangeRet.i = binarySearchIndex(xs, ys[0]);
    _maxIntRangeRet.j = binarySearchIndex(ys, xs[0]);
    _maxIntRangeRet.endA = Math.min(binarySearchIndex(xs, ys[lb]), la);
    _maxIntRangeRet.endB = Math.min(binarySearchIndex(ys, xs[la]), lb);
    return _maxIntRangeRet;
}

const _startEndRet = { start: 0, end: 0 };

function getStartEnd(set: OrderedSet, min: number, max: number) {
    _startEndRet.start = S.getInsertionIndex(set, min);
    let end = S.getInsertionIndex(set, max);
    if (end < S.size(set) && S.elementAt(set, end) === max) end++;
    _startEndRet.end = end;
    return _startEndRet;
}

function equalAA(a: Nums, b: Nums) {
    let size = a.length;
    if (a.length !== b.length || a[0] !== b[0] || a[size - 1] !== b[size - 1]) return false;
    for (let i = 0; i < size; i++) {
        if (a[i] !== b[i]) return false;
    }
    return true;
}

function areIntersectingAA(xs: Nums, ys: Nums) {
    let { i, j, endA, endB } = getMaxIntersectionRange(xs, ys);
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else return true;
    }
    return false;
}

function isSubsetAA(xs: Nums, ys: Nums) {
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
    return S.size(a) > 0 && S.size(b) > 0 && S.max(a) >= S.min(b) && S.min(a) <= S.max(b);
}

function isRangeSubset(a: OrderedSet, b: OrderedSet) {
    if (!S.size(a)) return S.size(b) === 0;
    if (!S.size(b)) return true;
    return S.min(a) <= S.min(b) && S.max(a) >= S.max(b);
}

function unionRR(a: number, b: number) {
    const sizeA = S.size(a), sizeB = S.size(b);
    if (!sizeA) return b;
    if (!sizeB) return a;
    const minA = minR(a), minB = minR(b);
    if (areRangesIntersecting(a, b)) return S.ofRange(Math.min(minA, minB), Math.max(maxR(a), maxR(b)));
    let lSize, lMin, rSize, rMin;
    if (minR(a) < minR(b)) { lSize = sizeA; lMin = minA; rSize = sizeB; rMin = minB; }
    else { lSize = sizeB; lMin = minB; rSize = sizeA; rMin = minA; }
    const arr = new Int32Array(sizeA + sizeB);
    for (let i = 0; i < lSize; i++) arr[i] = i + lMin;
    for (let i = 0; i < rSize; i++) arr[i + lSize] = i + rMin;
    return S.ofSortedArray(arr);
}

const _uAR = IntPair.zero();
function unionAR(a: Nums, b: number) {
    const bSize = S.size(b);
    if (!bSize) return a;
    // is the array fully contained in the range?
    if (isRangeSubset(b, a)) return b;

    IntPair.unpack(b, _uAR);
    const min = _uAR.fst, max = _uAR.snd;
    const { start, end } = getStartEnd(a, min, max);

    const size = start + (a.length - end) + bSize;
    const indices = new Int32Array(size);
    let offset = 0;
    for (let i = 0; i < start; i++) indices[offset++] = a[i];
    for (let i = min; i <= max; i++) indices[offset++] = i;
    for (let i = end, _i = a.length; i < _i; i++) indices[offset] = a[i];

    return OrderedSet.ofSortedArray(indices);
}

function unionAA(a: Nums, b: Nums) {
    const lenA = a.length, lenB = b.length;

    let { i: sI, j: sJ, endA, endB } = getMaxIntersectionRange(a, b);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i <= endA && j <= endB) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; commonCount++; }
    }

    if (!commonCount) return a;
    if (commonCount >= lenA) return OrderedSet.Empty

    const resultSize = lenA + lenB - commonCount;
    const l = Math.min(a[0], b[0]), r = Math.max(a[lenA - 1], b[lenB - 1]);
    // is this just a range?
    if (resultSize === r - l + 1) {
        return OrderedSet.ofRange(l, r);
    }

    const indices = new Int32Array(lenA + lenB - commonCount);
    let offset = 0;

    // insert the "prefixes"
    for (let k = 0; k < sI; k++) indices[offset++] = a[k];
    for (let k = 0; k < sJ; k++) indices[offset++] = a[k];

    // insert the common part
    i = sI;
    j = sJ;
    while (i <= endA && j <= endB) {
        const x = a[i], y = b[j];
        if (x < y) { indices[offset++] = x; i++; }
        else if (x > y) { indices[offset++] = y; j++; }
        else { indices[offset++] = x; i++; j++; }
    }

    // insert the "tail"
    for (; i < lenA; i++) indices[offset++] = a[i];
    for (; j < lenB; j++) indices[offset++] = b[j];

    return OrderedSet.ofSortedArray(indices);
}

const _iRA = IntPair.zero(), _iRB = IntPair.zero();
function intersectRR(a: number, b: number) {
    if (!areRangesIntersecting(a, b)) return OrderedSet.Empty;

    IntPair.unpack(a, _iRA);
    IntPair.unpack(b, _iRB);
    return OrderedSet.ofRange(Math.max(_iRA.fst, _iRB.fst), Math.min(_iRA.snd, _iRB.snd));
}

const _iAR = IntPair.zero();
function intersectAR(a: Nums, r: number) {
    if (!S.size(r)) return OrderedSet.Empty;

    IntPair.unpack(r, _iAR);
    const { start, end } = getStartEnd(a, _iAR.fst, _iAR.snd);
    const resultSize = end - start;
    if (!resultSize) return OrderedSet.Empty;

    const indices = new Int32Array(resultSize);
    let offset = 0;
    for (let i = start; i < end; i++) {
        indices[offset++] = a[i];
    }
    return OrderedSet.ofSortedArray(indices);
}

function intersectAA(xs: Nums, ys: Nums) {
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

const _sRA = IntPair.zero(), _sRB = IntPair.zero();
function substractRR(a: number, b: number) {
    IntPair.unpack(a, _sRA);
    IntPair.unpack(b, _sRB);

    if (_sRA.snd < _sRA.fst || _sRB.snd < _sRB.fst) return a;
    // is A subset of B? ==> Empty
    if (isRangeSubset(b, a)) return OrderedSet.Empty;
    if (isRangeSubset(a, b)) {
        // this splits the interval into two, gotta represent it as a set.
        const l = _sRB.fst - _sRA.fst, r = _sRA.snd - _sRB.snd;
        const ret = new Int32Array(l + r);
        let offset = 0;
        for (let i = 0; i < l; i++) ret[offset++] = _sRA.fst + i;
        for (let i = 1; i <= r; i++) ret[offset++] = _sRB.snd + i;
        return OrderedSet.ofSortedArray(ret);
    }
    // non intersecting ranges are handled by top-level substract.
    // at this point, b either contains rA.fst or rA.snd, but not both.
    if (_sRA.fst < _sRB.fst) return OrderedSet.ofRange(_sRA.fst, _sRB.fst - 1);
    return OrderedSet.ofRange(_sRB.snd + 1, _sRA.snd);
}

const _sAR = IntPair.zero();
function subtractAR(a: Nums, r: number) {
    IntPair.unpack(r, _sAR);
    if (_sAR.snd < _sAR.fst) return a;

    const min = _sAR.fst, max = _sAR.snd;
    const { start, end } = getStartEnd(a, min, max);
    const size = a.length - (end - start);
    if (size <= 0) return OrderedSet.Empty;
    const ret = new Int32Array(size);
    let offset = 0;
    for (let i = 0; i < start; i++) ret[offset++] = a[i];
    for (let i = end, _i = a.length; i < _i; i++) ret[offset++] = a[i];
    return OrderedSet.ofSortedArray(ret);
}

const _sAR1 = IntPair.zero();
function subtractRA(r: number, b: Nums) {
    IntPair.unpack(r, _sAR1);
    if (_sAR1.snd < _sAR1.fst) return r;

    const min = _sAR1.fst, max = _sAR1.snd;
    const rSize = max - min + 1;
    const { start, end } = getStartEnd(b, min, max);
    const commonCount = end - start;
    const resultSize = rSize - commonCount;
    if (resultSize <= 0) return OrderedSet.Empty;
    const ret = new Int32Array(resultSize);
    const li = b.length - 1;
    const fst = b[Math.min(start, li)], last = b[Math.min(end, li)];
    let offset = 0;
    for (let i = min; i < fst; i++) ret[offset++] = i;
    for (let i = fst; i <= last; i++) {
        if (binarySearch(b, i) < 0) ret[offset++] = i;
    }
    for (let i = last + 1; i <= max; i++) ret[offset++] = i;
    return OrderedSet.ofSortedArray(ret);
}

function subtractAA(a: Nums, b: Nums) {
    const lenA = a.length;

    let { i: sI, j: sJ, endA, endB } = getMaxIntersectionRange(a, b);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i <= endA && j <= endB) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; commonCount++; }
    }

    if (!commonCount) return a;
    if (commonCount >= lenA) return OrderedSet.Empty;

    const indices = new Int32Array(lenA - commonCount);
    let offset = 0;

    // insert the "prefix"
    for (let k = 0; k < sI; k++) indices[offset++] = a[k];

    i = sI;
    j = sJ;
    while (i <= endA && j <= endB) {
        const x = a[i], y = b[j];
        if (x < y) { indices[offset++] = x; i++; }
        else if (x > y) { j++; }
        else { i++; j++; }
    }

    // insert the "tail"
    for (; i < lenA; i++) indices[offset++] = a[i];

    return OrderedSet.ofSortedArray(indices);
}

export default OrderedSet