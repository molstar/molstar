/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import S from '../sorted-array'
import I from '../interval'

type OrderedSetImpl = I | S
type Nums = ArrayLike<number>

export const Empty: OrderedSetImpl = I.Empty;

export const ofSingleton = I.ofSingleton
export const ofRange = I.ofRange
export const ofBounds = I.ofBounds

export function ofSortedArray(xs: Nums): OrderedSetImpl {
    if (!xs.length) return Empty;
    // check if the array is just a range
    if (xs[xs.length - 1] - xs[0] + 1 === xs.length) return I.ofRange(xs[0], xs[xs.length - 1]);
    return xs as any;
}

export function size(set: OrderedSetImpl) { return I.is(set) ? I.size(set) : S.size(set); }
export function has(set: OrderedSetImpl, x: number) { return I.is(set) ? I.has(set, x) : S.has(set, x); }
export function indexOf(set: OrderedSetImpl, x: number) { return I.is(set) ? I.indexOf(set, x) : S.indexOf(set, x); }
export function getAt(set: OrderedSetImpl, i: number) { return I.is(set) ? I.getAt(set, i) : S.getAt(set, i); }
export function min(set: OrderedSetImpl) { return I.is(set) ? I.min(set) : S.min(set); }
export function max(set: OrderedSetImpl) { return I.is(set) ? I.max(set) : S.max(set); }

export function hashCode(set: OrderedSetImpl) { return I.is(set) ? I.hashCode(set) : S.hashCode(set); }
// TODO: possibly add more hash functions to allow for multilevel hashing.

export function areEqual(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (I.is(a)) {
        if (I.is(b)) return I.areEqual(a, b);
        return areEqualIS(a, b);
    } else if (I.is(b)) return areEqualIS(b, a);
    return S.areEqual(a, b);
}

export function areIntersecting(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (I.is(a)) {
        if (I.is(b)) return I.areIntersecting(a, b);
        return areIntersectingSI(b, a);
    } else if (I.is(b)) return areIntersectingSI(a, b);
    return areIntersectingSS(a, b);
}

/** Check if the 2nd argument is a subset of the 1st */
export function isSubset(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (I.is(a)) {
        if (I.is(b)) return I.isSubInterval(a, b);
        return isSubsetIS(a, b);
    } else if (I.is(b)) return isSubsetSI(a, b);
    return isSubsetSS(a, b);
}

export function findPredecessorIndex(set: OrderedSetImpl, x: number) {
    return I.is(set) ? I.findPredecessorIndex(set, x) : S.findPredecessorIndex(set, x);
}

export function findPredecessorIndexInInterval(set: OrderedSetImpl, x: number, bounds: I) {
    return I.is(set) ? I.findPredecessorIndexInInterval(set, x, bounds) : S.findPredecessorIndexInInterval(set, x, bounds);
}

export function findRange(set: OrderedSetImpl, min: number, max: number) {
    return I.is(set) ? I.findRange(set, min, max) : S.findRange(set, min, max);
}

export function union(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (I.is(a)) {
        if (I.is(b)) return unionII(a, b);
        return unionSI(b, a);
    } else if (I.is(b)) return unionSI(a, b);
    return unionSS(a, b);
}

export function intersect(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (I.is(a)) {
        if (I.is(b)) return I.intersect(a, b);
        return intersectSI(b, a);
    } else if (I.is(b)) return intersectSI(a, b);
    return intersectSS(a, b);
}

export function subtract(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (I.is(a)) {
        if (I.is(b)) return subtractII(a, b);
        return subtractIS(a, b);
    } else if (I.is(b)) return subtractSI(a, b);
    return subtractSS(a, b);
}

const _maxIntRangeRet = { startI: 0, startJ: 0, endI: 0, endJ: 0 };
// for small sets, just gets the whole range, for large sets does a bunch of binary searches
function getSuitableIntersectionRange(a: S, b: S) {
    const la = a.length, lb = b.length;
    const ratio = la / lb;
    if (ratio <= 0.5 || ratio >= 2 || (la >= 128 && lb >= 128)) {
        _maxIntRangeRet.startI = S.findPredecessorIndex(a, S.start(b));
        _maxIntRangeRet.startJ = S.findPredecessorIndex(b, S.start(a));
        _maxIntRangeRet.endI = S.findPredecessorIndex(a, S.end(b));
        _maxIntRangeRet.endJ = S.findPredecessorIndex(b, S.end(a));
    } else {
        _maxIntRangeRet.startI = 0;
        _maxIntRangeRet.startJ = 0;
        _maxIntRangeRet.endI = la;
        _maxIntRangeRet.endJ = lb;
    }
    return _maxIntRangeRet;
}

function areEqualIS(a: I, b: S) { return I.size(a) === S.size(b) && I.start(a) === S.start(b) && I.end(a) === S.end(b); }

function areIntersectingSI(a: S, b: I) {
    return areRangesIntersecting(a, b);
}

function areIntersectingSS(a: S, b: S) {
    if (a === b) return true;

    let { startI: i, startJ: j, endI, endJ } = getSuitableIntersectionRange(a, b);
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) { i++; }
        else if (x > y) { j++; }
        else return true;
    }
    return false;
}

function isSubsetSI(a: S, b: I) {
    const minB = I.min(b), maxB = I.max(b);
    if (maxB - minB + 1 === 0) return true;
    const minA = S.min(a), maxA = S.max(a);
    if (minB < minA || maxB > maxA) return false;
    const r = S.findRange(a, minB, maxB);
    return I.size(r) === I.size(b);
}

function isSubsetIS(a: I, b: S) {
    const minA = I.min(a), maxA = I.max(a);
    if (maxA - minA + 1 === 0) return false;
    const minB = S.min(b), maxB = S.max(b);
    return minB >= minA && maxA <= maxB;
}

function isSubsetSS(a: S, b: S) {
    if (a === b) return true;

    const lenB = b.length;
    let { startI: i, startJ: j, endI, endJ } = getSuitableIntersectionRange(a, b);
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
    const sa = size(a), sb = size(b);
    if (sa === 0 && sb === 0) return true;
    return sa > 0 && sb > 0 && max(a) >= min(b) && min(a) <= max(b);
}

function isRangeSubset(a: OrderedSetImpl, b: OrderedSetImpl) {
    if (!size(a)) return size(b) === 0;
    if (!size(b)) return true;
    return min(a) <= min(b) && max(a) >= max(b);
}

function unionII(a: I, b: I) {
    if (I.areEqual(a, b)) return a;

    const sizeA = I.size(a), sizeB = I.size(b);
    if (!sizeA) return b;
    if (!sizeB) return a;
    const minA = I.min(a), minB = I.min(b);
    if (areRangesIntersecting(a, b)) return I.ofRange(Math.min(minA, minB), Math.max(I.max(a), I.max(b)));
    let lSize, lMin, rSize, rMin;
    if (minA < minB) { lSize = sizeA; lMin = minA; rSize = sizeB; rMin = minB; }
    else { lSize = sizeB; lMin = minB; rSize = sizeA; rMin = minA; }
    const arr = new Int32Array(sizeA + sizeB);
    for (let i = 0; i < lSize; i++) arr[i] = i + lMin;
    for (let i = 0; i < rSize; i++) arr[i + lSize] = i + rMin;
    return ofSortedArray(arr);
}

function unionSI(a: S, b: I) {
    const bSize = I.size(b);
    if (!bSize) return a;
    // is the array fully contained in the range?
    if (isRangeSubset(b, a)) return b;

    const min = I.min(b), max = I.max(b);
    const r = S.findRange(a, min, max);
    const start = I.start(r), end = I.end(r);
    const indices = new Int32Array(start + (a.length - end) + bSize);
    let offset = 0;
    for (let i = 0; i < start; i++) indices[offset++] = a[i];
    for (let i = min; i <= max; i++) indices[offset++] = i;
    for (let i = end, _i = a.length; i < _i; i++) indices[offset] = a[i];

    return ofSortedArray(indices);
}

function unionSS(a: S, b: S) {
    if (a === b) return a;

    const { startI: sI, startJ: sJ, endI, endJ } = getSuitableIntersectionRange(a, b);
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
        return I.ofRange(l, r);
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

function intersectSI(a: S, b: I) {
    if (!I.size(b)) return Empty;

    const r = S.findRange(a, I.min(b), I.max(b));
    const start = I.start(r), end = I.end(r);
    const resultSize = end - start;
    if (!resultSize) return Empty;

    const indices = new Int32Array(resultSize);
    let offset = 0;
    for (let i = start; i < end; i++) {
        indices[offset++] = a[i];
    }
    return ofSortedArray(indices);
}

function intersectSS(a: S, b: S) {
    if (a === b) return a;

    const { startI: sI, startJ: sJ, endI, endJ } = getSuitableIntersectionRange(a, b);
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

function subtractII(a: I, b: I) {
    if (I.areEqual(a, b)) return Empty;
    if (!I.areIntersecting(a, b)) return a;

    const minA = I.min(a), maxA = I.max(a);
    const minB = I.min(b), maxB = I.max(b);

    if (maxA < minA || maxB < minB) return a;
    // is A subset of B? ==> Empty
    if (I.isSubInterval(b, a)) return Empty;
    if (I.isSubInterval(a, b)) {
        // this splits the interval into two, gotta represent it as a set.
        const l = minB - minA, r = maxA - maxB;
        if (l <= 0) return I.ofRange(maxB + 1, maxB + r);
        if (r <= 0) return I.ofRange(minA, minA + l - 1);
        const ret = new Int32Array(l + r);
        let offset = 0;
        for (let i = 0; i < l; i++) ret[offset++] = minA + i;
        for (let i = 1; i <= r; i++) ret[offset++] = maxB + i;
        return ofSortedArray(ret);
    }
    if (minA < minB) return I.ofRange(minA, minB - 1);
    return I.ofRange(maxB + 1, maxA);
}

function subtractSI(a: S, b: I) {
    const min = I.min(b), max = I.max(b);
    // is empty?
    if (max < min) return a;

    const r = S.findRange(a, min, max);
    const start = I.start(r), end = I.end(r);
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

function subtractIS(a: I, b: S) {
    const min = I.min(a), max = I.max(a);

    // is empty?
    if (max < min) return a;

    const rSize = max - min + 1;
    const interval = S.findRange(b, min, max);
    const start = I.start(interval), end = I.end(interval);
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
        if (S.indexOfInterval(b, i, interval) < 0) ret[offset++] = i;
    }
    for (let i = last + 1; i <= max; i++) ret[offset++] = i;
    return ofSortedArray(ret);
}

function subtractSS(a: S, b: S) {
    if (a === b) return Empty;

    const lenA = a.length;
    const { startI: sI, startJ: sJ, endI, endJ } = getSuitableIntersectionRange(a, b);
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