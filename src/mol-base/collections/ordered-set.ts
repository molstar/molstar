/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import IntTuple from './int-tuple'
import { hash3, hash4 } from './hash-functions'

/** An immutable ordered set. */
interface OrderedSet { '@type': 'int-ordered-set' }

namespace OrderedSet {
    export const Empty: OrderedSet = IntTuple.pack(0, -1) as any;
    export function ofSingleton(value: number): OrderedSet { return IntTuple.pack(value, value) as any; }
    export function ofRange(min: number, max: number): OrderedSet { return max < min ? Empty : IntTuple.pack(min, max) as any; }
    /** It is the responsibility of the caller to ensure the array is sorted and contains unique values. */
    export function ofSortedArray(xs: SortedArray): OrderedSet {
        if (!xs.length) return Empty;
        // check if the array is just a range
        if (xs[xs.length - 1] - xs[0] + 1 === xs.length) return ofRange(xs[0], xs[xs.length - 1]);
        return xs as any;
    }

    export const has: (set: OrderedSet, x: number) => boolean = hasI as any;
    export const indexOf: (set: OrderedSet, x: number) => number = indexOfI as any;
    export const getAt: (set: OrderedSet, i: number) => number = getAtI as any;

    export const min: (set: OrderedSet) => number = minI as any;
    export const max: (set: OrderedSet) => number = maxI as any;
    export const size: (set: OrderedSet) => number = sizeI as any;
    export const hashCode: (set: OrderedSet) => number = hashCodeI as any;

    export const areEqual: (a: OrderedSet, b: OrderedSet) => boolean = areEqualI as any;
    export const areIntersecting: (a: OrderedSet, b: OrderedSet) => boolean = areIntersectingI as any;
    export const isSubset: (a: OrderedSet, b: OrderedSet) => boolean = isSubsetI as any;

    export const union: (a: OrderedSet, b: OrderedSet) => OrderedSet = unionI as any;
    export const intersect: (a: OrderedSet, b: OrderedSet) => OrderedSet = intersectI as any;
    export const subtract: (a: OrderedSet, b: OrderedSet) => OrderedSet = subtractI as any;

    export const getInsertionIndex: (set: OrderedSet, x: number) => number = getInsertionIndexI as any;
    export const getIntervalRange: (set: OrderedSet, min: number, max: number) => { start: number, end: number } = getIntervalRangeI as any;
}

export default OrderedSet

/** Long and painful implementation starts here */

type Range = IntTuple
type SortedArray = ArrayLike<number>
type OrderedSetImpl = Range | SortedArray

function sizeI(set: OrderedSetImpl) { return typeof set === 'number' ? sizeR(set) : (set as SortedArray).length; }
function hasI(set: OrderedSetImpl, x: number) { return typeof set === 'number' ? hasR(set, x) : hasA(set as SortedArray, x); }
function indexOfI(set: OrderedSetImpl, x: number) { return typeof set === 'number' ? indexOfR(set, x) : indexOfA(set as SortedArray, x); }
function getAtI(set: OrderedSetImpl, i: number) { return typeof set === 'number' ? elementAtR(set, i) : (set as SortedArray)[i]; }
function minI(set: OrderedSetImpl) { return typeof set === 'number' ? minR(set) : (set as SortedArray)[0]; }
function maxI(set: OrderedSetImpl) { return typeof set === 'number' ? maxR(set) : (set as SortedArray)[(set as SortedArray).length - 1]; }

function hashCodeI(set: OrderedSetImpl) {
    // hash of tuple (size, min, max, mid)
    const s = sizeI(set);
    if (!s) return 0;
    if (s > 2) return hash4(s, getAtI(set, 0), getAtI(set, s - 1), getAtI(set, s >> 1));
    return hash3(s, getAtI(set, 0), getAtI(set, s - 1));
}
// TODO: possibly add more hash functions to allow for multilevel hashing.

function areEqualI(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return equalRR(a, b);
        return false;
    } else if (typeof b === 'number') return false;
    return equalAA(a as SortedArray, b as SortedArray);
}

function areIntersectingI(a: OrderedSetImpl, b: OrderedSetImpl) {
    // if at least one is "range", they must now intersect
    if (typeof a === 'number') {
        if (typeof b === 'number') return equalRR(a, b) || areRangesIntersecting(a, b);
        return areRangesIntersecting(a, b);
    }
    if (!areRangesIntersecting(a, b)) return false;
    else if (typeof b === 'number') return false;
    return areIntersectingAA(a as SortedArray, b as SortedArray);
}

/** Check if the 2nd argument is a subset of the 1st */
function isSubsetI(set: OrderedSetImpl, toTest: OrderedSetImpl) {
    if (!isRangeSubset(set, toTest)) return false;
    const testSize = sizeI(toTest);
    if (typeof set === 'number' || !testSize) return true;
    if (typeof toTest === 'number') return indexOfI(set, maxR(toTest)) - indexOfI(set, minR(toTest)) + 1 === testSize;
    return isSubsetAA(set as SortedArray, toTest as SortedArray);
}

function getInsertionIndexI(set: OrderedSetImpl, x: number) {
    return typeof set === 'number' ? rangeSearchIndex(set, x) : binarySearchIndex(set as SortedArray, x);
}

function getIntervalRangeI(set: OrderedSetImpl, min: number, max: number) {
    const { start, end } = getStartEnd(set, min, max);
    return { start, end };
}

function unionI(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return unionRR(a, b);
        return unionAR(b as SortedArray, a);
    } else if (typeof b === 'number') {
        return unionAR(a as SortedArray, b);
    } else return unionAA(a as SortedArray, b as SortedArray);
}

function intersectI(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return intersectRR(a, b);
        return intersectAR(b as SortedArray, a);
    } else if (typeof b === 'number') {
        return intersectAR(a as SortedArray, b);
    } else {
        if (!areRangesIntersecting(a, b)) return OrderedSet.Empty;
        return intersectAA(a as SortedArray, b as SortedArray);
    }
}

function subtractI(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (!areRangesIntersecting(a, b)) return a;

    if (typeof a === 'number') {
        if (typeof b === 'number') return substractRR(a, b);
        return subtractRA(a, b as SortedArray);
    } else if (typeof b === 'number') {
        return subtractAR(a as SortedArray, b);
    } else {
        return subtractAA(a as SortedArray, b as SortedArray);
    }
}

const minR = IntTuple.fst
const maxR = IntTuple.snd
const equalRR = IntTuple.areEqual

const _eR = IntTuple.zero();
function sizeR(set: Range) { IntTuple.unpack(set, _eR); return _eR.snd - _eR.fst + 1; }
function hasR(set: Range, x: number) { IntTuple.unpack(set, _eR); return x >= _eR.fst && x <= _eR.snd; }
function indexOfR(set: Range, x: number) { IntTuple.unpack(set, _eR); return x >= _eR.fst && x <= _eR.snd ? x - _eR.fst : -1; }
function elementAtR(set: Range, i: number) { return IntTuple.fst(set) + i; }

function hasA(set: SortedArray, x: number) { return x >= set[0] && x <= set[set.length - 1] && binarySearch(set, x) >= 0; }
function indexOfA(set: SortedArray, x: number) { return x >= set[0] && x <= set[set.length - 1] ? binarySearch(set, x) : -1; }

function binarySearch(xs: SortedArray, value: number) {
    return binarySearchRange(xs, value, 0, xs.length);
}

function binarySearchRange(xs: SortedArray, value: number, start: number, end: number) {
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

function binarySearchIndex(xs: SortedArray, value: number) {
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

const _rsiR = IntTuple.zero();
function rangeSearchIndex(r: Range, value: number) {
    IntTuple.unpack(r, _rsiR);
    if (value < _rsiR.fst) return 0;
    if (value > _rsiR.snd) return _rsiR.snd - _rsiR.fst + 1;
    return value - _rsiR.fst;
}

const _maxIntRangeRet = { i: 0, j: 0, endA: 0, endB: 0 };
function getMaxIntersectionRange(xs: SortedArray, ys: SortedArray) {
    const la = xs.length - 1, lb = ys.length - 1;
    _maxIntRangeRet.i = binarySearchIndex(xs, ys[0]);
    _maxIntRangeRet.j = binarySearchIndex(ys, xs[0]);
    _maxIntRangeRet.endA = Math.min(binarySearchIndex(xs, ys[lb]), la);
    _maxIntRangeRet.endB = Math.min(binarySearchIndex(ys, xs[la]), lb);
    return _maxIntRangeRet;
}

const _startEndRet = { start: 0, end: 0 };

function getStartEnd(set: OrderedSetImpl, min: number, max: number) {
    _startEndRet.start = getInsertionIndexI(set, min);
    let end = getInsertionIndexI(set, max);
    if (end < sizeI(set) && getAtI(set, end) === max) end++;
    _startEndRet.end = end;
    return _startEndRet;
}

function equalAA(a: SortedArray, b: SortedArray) {
    if (a === b) return true;
    let size = a.length;
    if (a.length !== b.length || a[0] !== b[0] || a[size - 1] !== b[size - 1]) return false;
    for (let i = 0; i < size; i++) {
        if (a[i] !== b[i]) return false;
    }
    return true;
}

function areIntersectingAA(xs: SortedArray, ys: SortedArray) {
    if (xs === ys) return true;

    let { i, j, endA, endB } = getMaxIntersectionRange(xs, ys);
    while (i <= endA && j <= endB) {
        const x = xs[i], y = ys[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else return true;
    }
    return false;
}

function isSubsetAA(a: SortedArray, b: SortedArray) {
    if (a === b) return true;

    const lenB = b.length;
    let { i, j, endA, endB } = getMaxIntersectionRange(a, b);
    // must be able to advance by lenB elements
    if (endB - j + 1 < lenB || endA - i + 1 < lenB) return false;

    let equal = 0;
    while (i <= endA && j <= endB) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; equal++; }
    }
    return equal === lenB;
}

function areRangesIntersecting(a: OrderedSetImpl, b: OrderedSetImpl) {
    return sizeI(a) > 0 && sizeI(b) > 0 && maxI(a) >= minI(b) && minI(a) <= maxI(b);
}

function isRangeSubset(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (!sizeI(a)) return sizeI(b) === 0;
    if (!sizeI(b)) return true;
    return minI(a) <= minI(b) && maxI(a) >= maxI(b);
}

function unionRR(a: Range, b: Range) {
    if (IntTuple.areEqual(a, b)) return a;

    const sizeA = sizeR(a), sizeB = sizeR(b);
    if (!sizeA) return b;
    if (!sizeB) return a;
    const minA = minR(a), minB = minR(b);
    if (areRangesIntersecting(a, b)) return OrderedSet.ofRange(Math.min(minA, minB), Math.max(maxR(a), maxR(b)));
    let lSize, lMin, rSize, rMin;
    if (minR(a) < minR(b)) { lSize = sizeA; lMin = minA; rSize = sizeB; rMin = minB; }
    else { lSize = sizeB; lMin = minB; rSize = sizeA; rMin = minA; }
    const arr = new Int32Array(sizeA + sizeB);
    for (let i = 0; i < lSize; i++) arr[i] = i + lMin;
    for (let i = 0; i < rSize; i++) arr[i + lSize] = i + rMin;
    return OrderedSet.ofSortedArray(arr);
}

const _uAR = IntTuple.zero();
function unionAR(a: SortedArray, b: Range) {
    const bSize = sizeI(b);
    if (!bSize) return a;
    // is the array fully contained in the range?
    if (isRangeSubset(b, a)) return b;

    IntTuple.unpack(b, _uAR);
    const min = _uAR.fst, max = _uAR.snd;
    const { start, end } = getStartEnd(a, min, max);

    const indices = new Int32Array(start + (a.length - end) + bSize);
    let offset = 0;
    for (let i = 0; i < start; i++) indices[offset++] = a[i];
    for (let i = min; i <= max; i++) indices[offset++] = i;
    for (let i = end, _i = a.length; i < _i; i++) indices[offset] = a[i];

    return OrderedSet.ofSortedArray(indices);
}

function unionAA(a: SortedArray, b: SortedArray) {
    if (a === b) return a;

    let { i: sI, j: sJ, endA, endB } = getMaxIntersectionRange(a, b);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i <= endA && j <= endB) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; commonCount++; }
    }

    const lenA = a.length, lenB = b.length;
    // A === B || B is subset of A ==> A
    if ((commonCount === lenA && commonCount === lenB) || commonCount === lenB) return a;
    // A is subset of B ===> B
    if (commonCount === lenA) return b;

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

const _iRA = IntTuple.zero(), _iRB = IntTuple.zero();
function intersectRR(a: Range, b: Range) {
    if (!areRangesIntersecting(a, b)) return OrderedSet.Empty;
    if (IntTuple.areEqual(a, b)) return a;

    IntTuple.unpack(a, _iRA);
    IntTuple.unpack(b, _iRB);
    return OrderedSet.ofRange(Math.max(_iRA.fst, _iRB.fst), Math.min(_iRA.snd, _iRB.snd));
}

const _iAR = IntTuple.zero();
function intersectAR(a: SortedArray, r: Range) {
    if (!sizeI(r)) return OrderedSet.Empty;

    IntTuple.unpack(r, _iAR);
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

function intersectAA(a: SortedArray, b: SortedArray) {
    if (a === b) return a;

    let { i: sI, j: sJ, endA, endB } = getMaxIntersectionRange(a, b);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i <= endA && j <= endB) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; commonCount++; }
    }

    const lenA = a.length, lenB = b.length;
    // no common elements
    if (!commonCount) return OrderedSet.Empty;
    // A === B || B is subset of A ==> B
    if ((commonCount === lenA && commonCount === lenB) || commonCount === lenB) return b;
    // A is subset of B ==> A
    if (commonCount === lenA) return a;

    const indices = new Int32Array(commonCount);
    let offset = 0;
    i = sI;
    j = sJ;
    while (i <= endA && j <= endB) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { indices[offset++] = x; i++; j++; }
    }

    return OrderedSet.ofSortedArray(indices);
}

const _sRA = IntTuple.zero(), _sRB = IntTuple.zero();
function substractRR(a: Range, b: Range) {
    if (IntTuple.areEqual(a, b)) return OrderedSet.Empty;

    IntTuple.unpack(a, _sRA);
    IntTuple.unpack(b, _sRB);

    if (_sRA.snd < _sRA.fst || _sRB.snd < _sRB.fst) return a;
    // is A subset of B? ==> Empty
    if (isRangeSubset(b, a)) return OrderedSet.Empty;
    if (isRangeSubset(a, b)) {
        // this splits the interval into two, gotta represent it as a set.
        const l = _sRB.fst - _sRA.fst, r = _sRA.snd - _sRB.snd;
        if (l <= 0) return OrderedSet.ofRange(_sRB.snd + 1, _sRB.snd + r);
        if (r <= 0) return OrderedSet.ofRange(_sRA.fst, _sRA.fst + l - 1);
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

const _sAR = IntTuple.zero();
function subtractAR(a: SortedArray, b: Range) {
    IntTuple.unpack(b, _sAR);

    // is empty?
    if (_sAR.snd < _sAR.fst) return a;

    const min = _sAR.fst, max = _sAR.snd;
    const { start, end } = getStartEnd(a, min, max);
    const size = a.length - (end - start);
    // A is subset of B
    if (size <= 0) return OrderedSet.Empty;
    // No common elements
    if (size === a.length) return a;

    const ret = new Int32Array(size);
    let offset = 0;
    for (let i = 0; i < start; i++) ret[offset++] = a[i];
    for (let i = end, _i = a.length; i < _i; i++) ret[offset++] = a[i];
    return OrderedSet.ofSortedArray(ret);
}

const _sAR1 = IntTuple.zero();
function subtractRA(a: Range, b: SortedArray) {
    IntTuple.unpack(a, _sAR1);

    // is empty?
    if (_sAR1.snd < _sAR1.fst) return a;

    const min = _sAR1.fst, max = _sAR1.snd;
    const rSize = max - min + 1;
    const { start, end } = getStartEnd(b, min, max);
    const commonCount = end - start;

    // No common elements.
    if (commonCount === 0) return a;

    const resultSize = rSize - commonCount;
    // A is subset of B
    if (resultSize <= 0) return OrderedSet.Empty;

    const ret = new Int32Array(resultSize);
    const li = b.length - 1;
    const fst = b[Math.min(start, li)], last = b[Math.min(end, li)];
    let offset = 0;
    for (let i = min; i < fst; i++) ret[offset++] = i;
    for (let i = fst; i <= last; i++) {
        if (binarySearchRange(b, i, start, end) < 0) ret[offset++] = i;
    }
    for (let i = last + 1; i <= max; i++) ret[offset++] = i;
    return OrderedSet.ofSortedArray(ret);
}

function subtractAA(a: SortedArray, b: SortedArray) {
    if (a === b) return OrderedSet.Empty;

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

    // A isnt intersecting B ===> A
    if (!commonCount) return a;
    // A === B || A is subset of B ===> Empty
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