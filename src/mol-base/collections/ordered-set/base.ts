/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import IntTuple from '../int-tuple'
import { hash3, hash4 } from '../hash-functions'

type Range = IntTuple
type SortedArray = ArrayLike<number>
type OrderedSetImpl = Range | SortedArray

export const Empty: OrderedSetImpl = IntTuple.create(0, -1) as any;
export function ofSingleton(value: number): OrderedSetImpl { return IntTuple.create(value, value) as any; }
export function ofRange(min: number, max: number): OrderedSetImpl { return max < min ? Empty : IntTuple.create(min, max) as any; }
/** It is the responsibility of the caller to ensure the array is sorted and contains unique values. */
export function ofSortedArray(xs: SortedArray): OrderedSetImpl {
    if (!xs.length) return Empty;
    // check if the array is just a range
    if (xs[xs.length - 1] - xs[0] + 1 === xs.length) return ofRange(xs[0], xs[xs.length - 1]);
    return xs as any;
}

export function size(set: OrderedSetImpl) { return typeof set === 'number' ? sizeR(set) : (set as SortedArray).length; }
export function has(set: OrderedSetImpl, x: number) { return typeof set === 'number' ? hasR(set, x) : hasA(set as SortedArray, x); }
export function indexOf(set: OrderedSetImpl, x: number) { return typeof set === 'number' ? indexOfR(set, x) : indexOfA(set as SortedArray, x); }
export function getAt(set: OrderedSetImpl, i: number) { return typeof set === 'number' ? elementAtR(set, i) : (set as SortedArray)[i]; }
export function minValue(set: OrderedSetImpl) { return typeof set === 'number' ? minR(set) : (set as SortedArray)[0]; }
export function maxValue(set: OrderedSetImpl) { return typeof set === 'number' ? maxR(set) : (set as SortedArray)[(set as SortedArray).length - 1]; }

export function hashCode(set: OrderedSetImpl) {
    // hash of tuple (size, min, max, mid)
    const s = size(set);
    if (!s) return 0;
    if (s > 2) return hash4(s, getAt(set, 0), getAt(set, s - 1), getAt(set, s >> 1));
    return hash3(s, getAt(set, 0), getAt(set, s - 1));
}
// TODO: possibly add more hash functions to allow for multilevel hashing.

export function areEqual(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return equalRR(a, b);
        return false;
    } else if (typeof b === 'number') return false;
    return equalAA(a as SortedArray, b as SortedArray);
}

export function areIntersecting(a: OrderedSetImpl, b: OrderedSetImpl) {
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
export function isSubset(set: OrderedSetImpl, toTest: OrderedSetImpl) {
    if (!isRangeSubset(set, toTest)) return false;
    const testSize = size(toTest);
    if (typeof set === 'number' || !testSize) return true;
    if (typeof toTest === 'number') return indexOf(set, maxR(toTest)) - indexOf(set, minR(toTest)) + 1 === testSize;
    return isSubsetAA(set as SortedArray, toTest as SortedArray);
}

export function getPredIndex(set: OrderedSetImpl, x: number) {
    return typeof set === 'number' ? rangeSearchIndex(set, x) : binarySearchPredIndex(set as SortedArray, x);
}

export function getPredIndexInRange(set: OrderedSetImpl, x: number, { start, end }: { start: number, end: number }) {
    if (typeof set === 'number') {
        const ret = rangeSearchIndex(set, x);
        return ret <= start ? start : ret >= end ? end : ret;
     } else return binarySearchPredIndexRange(set as SortedArray, x, start, end);
}

export function getIntervalRange(set: OrderedSetImpl, min: number, max: number) {
    const { start, end } = getStartEnd(set, min, max);
    return { start, end };
}

export function union(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return unionRR(a, b);
        return unionAR(b as SortedArray, a);
    } else if (typeof b === 'number') {
        return unionAR(a as SortedArray, b);
    } else return unionAA(a as SortedArray, b as SortedArray);
}

export function intersect(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return intersectRR(a, b);
        return intersectAR(b as SortedArray, a);
    } else if (typeof b === 'number') {
        return intersectAR(a as SortedArray, b);
    } else {
        if (!areRangesIntersecting(a, b)) return Empty;
        return intersectAA(a as SortedArray, b as SortedArray);
    }
}

export function subtract(a: OrderedSetImpl, b: OrderedSetImpl) {
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

function sizeR(set: Range) { return maxR(set) - minR(set) + 1; }
function hasR(set: Range, x: number) { return x >= minR(set) && x <= maxR(set); }
function indexOfR(set: Range, x: number) { const m = minR(set); return x >= m && x <= maxR(set) ? x - m : -1; }
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

function binarySearchPredIndex(xs: SortedArray, value: number) {
    return binarySearchPredIndexRange(xs, value, 0, xs.length);
}

function binarySearchPredIndexRange(xs: SortedArray, value: number, start: number, end: number) {
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

function rangeSearchIndex(r: Range, value: number) {
    const min = minR(r);
    if (value < min) return 0;
    const max = maxR(r);
    if (value > max) return max - min + 1;
    return value - min;
}

const _maxIntRangeRet = { startI: 0, startJ: 0, endI: 0, endJ: 0 };
function getMaxIntersectionRange(xs: SortedArray, ys: SortedArray) {
    _maxIntRangeRet.startI = binarySearchPredIndex(xs, ys[0]);
    _maxIntRangeRet.startJ = binarySearchPredIndex(ys, xs[0]);
    _maxIntRangeRet.endI = binarySearchPredIndex(xs, ys[ys.length - 1] + 1);
    _maxIntRangeRet.endJ = binarySearchPredIndex(ys, xs[xs.length - 1] + 1);

    return _maxIntRangeRet;
}

const _startEndRet = { start: 0, end: 0 };
function getStartEnd(set: OrderedSetImpl, min: number, max: number) {
    _startEndRet.start = getPredIndex(set, min);
    _startEndRet.end = getPredIndex(set, max + 1);
    return _startEndRet;
}

function equalAA(a: SortedArray, b: SortedArray) {
    if (a === b) return true;
    const aSize = a.length;
    if (aSize !== b.length || a[0] !== b[0] || a[aSize - 1] !== b[aSize - 1]) return false;
    for (let i = 0; i < aSize; i++) {
        if (a[i] !== b[i]) return false;
    }
    return true;
}

function areIntersectingAA(xs: SortedArray, ys: SortedArray) {
    if (xs === ys) return true;

    let { startI: i, startJ: j, endI, endJ } = getMaxIntersectionRange(xs, ys);
    while (i < endI && j < endJ) {
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
    let { startI: i, startJ: j, endI, endJ } = getMaxIntersectionRange(a, b);
    // must be able to advance by lenB elements
    if (endJ - j < lenB || endI - i < lenB) return false;

    let equal = 0;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; equal++; }
    }
    return equal === lenB;
}

function areRangesIntersecting(a: OrderedSetImpl, b: OrderedSetImpl) {
    return size(a) > 0 && size(b) > 0 && maxValue(a) >= minValue(b) && minValue(a) <= maxValue(b);
}

function isRangeSubset(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (!size(a)) return size(b) === 0;
    if (!size(b)) return true;
    return minValue(a) <= minValue(b) && maxValue(a) >= maxValue(b);
}

function unionRR(a: Range, b: Range) {
    if (IntTuple.areEqual(a, b)) return a;

    const sizeA = sizeR(a), sizeB = sizeR(b);
    if (!sizeA) return b;
    if (!sizeB) return a;
    const minA = minR(a), minB = minR(b);
    if (areRangesIntersecting(a, b)) return ofRange(Math.min(minA, minB), Math.max(maxR(a), maxR(b)));
    let lSize, lMin, rSize, rMin;
    if (minR(a) < minR(b)) { lSize = sizeA; lMin = minA; rSize = sizeB; rMin = minB; }
    else { lSize = sizeB; lMin = minB; rSize = sizeA; rMin = minA; }
    const arr = new Int32Array(sizeA + sizeB);
    for (let i = 0; i < lSize; i++) arr[i] = i + lMin;
    for (let i = 0; i < rSize; i++) arr[i + lSize] = i + rMin;
    return ofSortedArray(arr);
}

function unionAR(a: SortedArray, b: Range) {
    const bSize = size(b);
    if (!bSize) return a;
    // is the array fully contained in the range?
    if (isRangeSubset(b, a)) return b;

    const min = minR(b), max = maxR(b);
    const { start, end } = getStartEnd(a, min, max);

    const indices = new Int32Array(start + (a.length - end) + bSize);
    let offset = 0;
    for (let i = 0; i < start; i++) indices[offset++] = a[i];
    for (let i = min; i <= max; i++) indices[offset++] = i;
    for (let i = end, _i = a.length; i < _i; i++) indices[offset] = a[i];

    return ofSortedArray(indices);
}

function unionAA(a: SortedArray, b: SortedArray) {
    if (a === b) return a;

    const { startI: sI, startJ: sJ, endI, endJ } = getMaxIntersectionRange(a, b);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i < endI && j < endJ) {
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
        return ofRange(l, r);
    }

    const indices = new Int32Array(lenA + lenB - commonCount);
    let offset = 0;

    // insert the "prefixes"
    for (let k = 0; k < sI; k++) indices[offset++] = a[k];
    for (let k = 0; k < sJ; k++) indices[offset++] = b[k];

    // insert the common part
    i = sI;
    j = sJ;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) { indices[offset++] = x; i++; }
        else if (x > y) { indices[offset++] = y; j++; }
        else { indices[offset++] = x; i++; j++; }
    }

    // insert the "tail"
    for (; i < lenA; i++) indices[offset++] = a[i];
    for (; j < lenB; j++) indices[offset++] = b[j];

    return ofSortedArray(indices);
}

function intersectRR(a: Range, b: Range) {
    if (!areRangesIntersecting(a, b)) return Empty;
    if (IntTuple.areEqual(a, b)) return a;
    return ofRange(Math.max(minR(a), minR(b)), Math.min(maxR(a), maxR(b)));
}

function intersectAR(a: SortedArray, r: Range) {
    if (!size(r)) return Empty;

    const { start, end } = getStartEnd(a, minR(r), maxR(r));
    const resultSize = end - start;
    if (!resultSize) return Empty;

    const indices = new Int32Array(resultSize);
    let offset = 0;
    for (let i = start; i < end; i++) {
        indices[offset++] = a[i];
    }
    return ofSortedArray(indices);
}

function intersectAA(a: SortedArray, b: SortedArray) {
    if (a === b) return a;

    const { startI: sI, startJ: sJ, endI, endJ } = getMaxIntersectionRange(a, b);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; commonCount++; }
    }

    const lenA = a.length, lenB = b.length;
    // no common elements
    if (!commonCount) return Empty;
    // A === B || B is subset of A ==> B
    if ((commonCount === lenA && commonCount === lenB) || commonCount === lenB) return b;
    // A is subset of B ==> A
    if (commonCount === lenA) return a;

    const indices = new Int32Array(commonCount);
    let offset = 0;
    i = sI;
    j = sJ;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { indices[offset++] = x; i++; j++; }
    }

    return ofSortedArray(indices);
}

function substractRR(a: Range, b: Range) {
    if (IntTuple.areEqual(a, b)) return Empty;

    const minA = minR(a), maxA = maxR(a);
    const minB = minR(b), maxB = maxR(b);

    if (maxA < minA || maxB < minB) return a;
    // is A subset of B? ==> Empty
    if (isRangeSubset(b, a)) return Empty;
    if (isRangeSubset(a, b)) {
        // this splits the interval into two, gotta represent it as a set.
        const l = minB - minA, r = maxA - maxB;
        if (l <= 0) return ofRange(maxB + 1, maxB + r);
        if (r <= 0) return ofRange(minA, minA + l - 1);
        const ret = new Int32Array(l + r);
        let offset = 0;
        for (let i = 0; i < l; i++) ret[offset++] = minA + i;
        for (let i = 1; i <= r; i++) ret[offset++] = maxB + i;
        return ofSortedArray(ret);
    }
    // non intersecting ranges are handled by top-level substract.
    // at this point, b either contains rA.fst or rA.snd, but not both.
    if (minA < minB) return ofRange(minA, minB - 1);
    return ofRange(maxB + 1, maxA);
}

function subtractAR(a: SortedArray, b: Range) {
    // is empty?
    const min = minR(b), max = maxR(b);
    if (max < min) return a;

    const { start, end } = getStartEnd(a, min, max);
    const resultSize = a.length - (end - start);
    // A is subset of B
    if (resultSize <= 0) return Empty;
    // No common elements
    if (resultSize === a.length) return a;

    const ret = new Int32Array(resultSize);
    let offset = 0;
    for (let i = 0; i < start; i++) ret[offset++] = a[i];
    for (let i = end, _i = a.length; i < _i; i++) ret[offset++] = a[i];
    return ofSortedArray(ret);
}

function subtractRA(a: Range, b: SortedArray) {
    const min = minR(a), max = maxR(a);

    // is empty?
    if (max < min) return a;

    const rSize = max - min + 1;
    const { start, end } = getStartEnd(b, min, max);
    const commonCount = end - start;

    // No common elements.
    if (commonCount === 0) return a;

    const resultSize = rSize - commonCount;
    // A is subset of B
    if (resultSize <= 0) return Empty;

    const ret = new Int32Array(resultSize);
    const li = b.length - 1;
    const fst = b[Math.min(start, li)], last = b[Math.min(end, li)];
    let offset = 0;
    for (let i = min; i < fst; i++) ret[offset++] = i;
    for (let i = fst; i <= last; i++) {
        if (binarySearchRange(b, i, start, end) < 0) ret[offset++] = i;
    }
    for (let i = last + 1; i <= max; i++) ret[offset++] = i;
    return ofSortedArray(ret);
}

function subtractAA(a: SortedArray, b: SortedArray) {
    if (a === b) return Empty;

    const lenA = a.length;

    const { startI: sI, startJ: sJ, endI, endJ } = getMaxIntersectionRange(a, b);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else { i++; j++; commonCount++; }
    }

    // A isnt intersecting B ===> A
    if (!commonCount) return a;
    // A === B || A is subset of B ===> Empty
    if (commonCount >= lenA) return Empty;

    const indices = new Int32Array(lenA - commonCount);
    let offset = 0;

    // insert the "prefix"
    for (let k = 0; k < sI; k++) indices[offset++] = a[k];

    i = sI;
    j = sJ;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) { indices[offset++] = x; i++; }
        else if (x > y) { j++; }
        else { i++; j++; }
    }

    // insert the "tail"
    for (; i < lenA; i++) indices[offset++] = a[i];

    return ofSortedArray(indices);
}