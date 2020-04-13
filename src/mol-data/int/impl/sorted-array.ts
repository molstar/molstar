/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArray, hash3, hash4, createRangeArray } from '../../util';
import Interval from '../interval';

type Nums = ArrayLike<number>


export const Empty: Nums = [];

export function ofSingleton(v: number) { return [v]; }
export function ofSortedArray(xs: Nums) { return xs; }
export function ofUnsortedArray(xs: Nums) { sortArray(xs); return xs; }
export function ofRange(min: number, max: number) {
    if (max < min) return [];
    const ret = new Int32Array(max - min + 1);
    for (let i = min; i <= max; i++) ret[i - min] = i;
    return ret;
}
export function is(xs: any): xs is Nums { return xs && (Array.isArray(xs) || !!xs.buffer); }

export function start(xs: Nums) { return xs[0]; }
export function end(xs: Nums) { return xs[xs.length - 1] + 1;  }
export function min(xs: Nums) { return xs[0]; }
export function max(xs: Nums) { return xs[xs.length - 1]; }
export function size(xs: Nums) { return xs.length; }
export function hashCode(xs: Nums) {
    // hash of tuple (size, min, max, mid)
    const s = xs.length;
    if (!s) return 0;
    if (s > 2) return hash4(s, xs[0], xs[s - 1], xs[s >> 1]);
    return hash3(s, xs[0], xs[s - 1]);
}
export function toString(xs: Nums) {
    const s = xs.length;
    if (s > 5) return `[${xs[0]}, ${xs[1]}, ..., ${xs[s - 1]}], length ${s}`;
    return `[${(xs as number[]).join(', ')}]`;
}

/** Returns the index of `x` in `set` or -1 if not found. */
export function indexOf(xs: Nums, v: number) {
    const l = xs.length;
    return l === 0 ? -1 : xs[0] <= v && v <= xs[l - 1] ? binarySearchRange(xs, v, 0, l) : -1;
}
export function indexOfInInterval(xs: Nums, v: number, bounds: Interval) {
    return indexOfInRange(xs, v, Interval.start(bounds), Interval.end(bounds));
}
export function indexOfInRange(xs: Nums, v: number, s: number, e: number) {
    const l = xs.length;
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

/**
 * Returns 0 if `v` is smaller or equal the first element of `xs`
 * Returns length of `xs` if `v` is bigger than the last element of `xs`
 * Otherwise returns the first index where the value of `xs` is equal or bigger than `v`
 */
export function findPredecessorIndex(xs: Nums, v: number) {
    const len = xs.length;
    if (v <= xs[0]) return 0;
    if (v > xs[len - 1]) return len;
    return binarySearchPredIndexRange(xs, v, 0, len);
}

export function findPredecessorIndexInInterval(xs: Nums, v: number, bounds: Interval) {
    const s = Interval.start(bounds), e = Interval.end(bounds);
    const sv = xs[s];
    if (v <= sv) return s;
    if (e > s && v > xs[e - 1]) return e;
    // do a linear search if there are only 10 or less items remaining
    if (v - sv <= 11) return linearSearchPredInRange(xs, v, s + 1, e);
    return binarySearchPredIndexRange(xs, v, s, e);
}

export function findRange(xs: Nums, min: number, max: number) {
    return Interval.ofBounds(findPredecessorIndex(xs, min), findPredecessorIndex(xs, max + 1));
}

function binarySearchRange(xs: Nums, value: number, start: number, end: number) {
    let min = start, max = end - 1;
    while (min <= max) {
        // do a linear search if there are only 10 or less items remaining
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
        // do a linear search if there are only 10 or less items remaining
        if (min + 11 > max) {
            for (let i = min; i <= max; i++) {
                if (value <= xs[i]) return i;
            }
            return max + 1;
        }
        const mid = (min + max) >> 1;
        const v = xs[mid];
        if (value < v) max = mid - 1;
        else if (value > v) min = mid + 1;
        else return mid;
    }
    if (min > max) return max + 1;
    return xs[min] >= value ? min : min + 1;
}

function linearSearchPredInRange(xs: Nums, value: number, start: number, end: number) {
    for (let i = start; i < end; i++) {
        if (value <= xs[i]) return i;
    }
    return end;
}

export function areIntersecting(a: Nums, b: Nums) {
    if (a === b) return true;

    let { startI: i, startJ: j, endI, endJ } = getSuitableIntersectionRange(a, b);
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) i++;
        else if (x > y) j++;
        else return true;
    }
    return false;
}

export function isSubset(a: Nums, b: Nums) {
    if (a === b) return true;

    const lenB = b.length;
    let { startI: i, startJ: j, endI, endJ } = getSuitableIntersectionRange(a, b);
    // must be able to advance by lenB elements
    if (endJ - j < lenB || endI - i < lenB) return false;

    let equal = 0;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) {
            i++;
        } else if (x > y) {
            j++;
        } else {
            i++; j++; equal++;
        }
    }
    return equal === lenB;
}

export function union(a: Nums, b: Nums): Nums {
    if (a === b) return a;

    const lenA = a.length, lenB = b.length;
    if (lenA === 0) return b;
    if (lenB === 0) return a;
    if (a[0] > b[0]) return union(b, a);

    const { startI, startJ, endI, endJ } = getSuitableIntersectionRange(a, b);
    const commonCount = getCommonCount(a, b, startI, startJ, endI, endJ);

    // A === B || B is subset of A ==> A
    if ((commonCount === lenA && commonCount === lenB) || commonCount === lenB) return a;
    // A is subset of B ===> B
    if (commonCount === lenA) return b;

    const indices = new Int32Array(lenA + lenB - commonCount);
    let i = 0, j = 0, offset = 0;

    // insert the "prefixes"
    for (i = 0; i < startI; i++) indices[offset++] = a[i];
    while (j < endJ && a[startI] > b[j]) indices[offset++] = b[j++];

    // insert the common part
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) {
            indices[offset++] = x; i++;
        } else if (x > y) {
            indices[offset++] = y; j++;
        } else {
            indices[offset++] = x; i++; j++;
        }
    }

    // insert the remaining common part
    for (; i < endI; i++) indices[offset++] = a[i];
    for (; j < endJ; j++) indices[offset++] = b[j];

    // insert the "tail"
    for (; i < lenA; i++) indices[offset++] = a[i];
    for (; j < lenB; j++) indices[offset++] = b[j];

    return ofSortedArray(indices);
}

export function intersectionSize(a: Nums, b: Nums) {
    if (a === b) return size(a);
    const { startI, startJ, endI, endJ } = getSuitableIntersectionRange(a, b);
    return getCommonCount(a, b, startI, startJ, endI, endJ);
}

function getCommonCount(a: Nums, b: Nums, startI: number, startJ: number, endI: number, endJ: number) {
    let i = startI, j = startJ;
    let commonCount = 0;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) {
            i++;
        } else if (x > y) {
            j++;
        } else {
            i++; j++; commonCount++;
        }
    }
    return commonCount;
}

export function intersect(a: Nums, b: Nums) {
    if (a === b) return a;

    const { startI, startJ, endI, endJ } = getSuitableIntersectionRange(a, b);
    const commonCount = getCommonCount(a, b, startI, startJ, endI, endJ);

    const lenA = a.length, lenB = b.length;
    // no common elements
    if (!commonCount) return Empty;
    // A === B || B is subset of A ==> B
    if ((commonCount === lenA && commonCount === lenB) || commonCount === lenB) return b;
    // A is subset of B ==> A
    if (commonCount === lenA) return a;

    const indices = new Int32Array(commonCount);
    let offset = 0;
    let i = startI;
    let j = startJ;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) {
            i++;
        } else if (x > y) {
            j++;
        } else {
            indices[offset++] = x; i++; j++;
        }
    }

    return ofSortedArray(indices);
}

export function subtract(a: Nums, b: Nums) {
    if (a === b) return Empty;

    const lenA = a.length;
    const { startI: sI, startJ: sJ, endI, endJ } = getSuitableIntersectionRange(a, b);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) {
            i++;
        } else if (x > y) {
            j++;
        } else {
            i++; j++; commonCount++;
        }
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
        if (x < y) {
            indices[offset++] = x; i++;
        } else if (x > y) {
            j++;
        } else {
            i++; j++;
        }
    }

    // insert the "tail"
    for (; i < lenA; i++) indices[offset++] = a[i];

    return ofSortedArray(indices);
}

export function deduplicate(xs: Nums) {
    if (xs.length < 2) return xs;
    let count = 1;
    for (let i = 0, _i = xs.length - 1; i < _i; i++) {
        if (xs[i] !== xs[i + 1]) count++;
    }
    if (count === xs.length) return xs;
    const ret = new Int32Array(count);
    let o = 0;
    for (let i = 0, _i = xs.length - 1; i < _i; i++) {
        if (xs[i] !== xs[i + 1]) ret[o++] = xs[i];
    }
    ret[o] = xs[xs.length - 1];
    return ret;
}

export function indicesOf(a: Nums, b: Nums): Nums {
    if (a === b) return ofSortedArray(createRangeArray(0, a.length - 1));

    const { startI: sI, startJ: sJ, endI, endJ } = getSuitableIntersectionRange(a, b);
    let i = sI, j = sJ;
    let commonCount = 0;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) {
            i++;
        } else if (x > y) {
            j++;
        } else {
            i++; j++; commonCount++;
        }
    }

    const lenA = a.length;
    // no common elements
    if (!commonCount) return Empty;
    // A is subset of B ==> A
    if (commonCount === lenA) return ofSortedArray(createRangeArray(0, a.length - 1));

    const indices = new Int32Array(commonCount);
    let offset = 0;
    i = sI;
    j = sJ;
    while (i < endI && j < endJ) {
        const x = a[i], y = b[j];
        if (x < y) {
            i++;
        } else if (x > y) {
            j++;
        } else {
            indices[offset++] = i; i++; j++;
        }
    }

    return ofSortedArray(indices);
}

const _maxIntRangeRet = { startI: 0, startJ: 0, endI: 0, endJ: 0 };
// for small sets, just gets the whole range, for large sets does a bunch of binary searches
function getSuitableIntersectionRange(a: Nums, b: Nums) {
    const la = a.length, lb = b.length;
    const ratio = la / lb;
    if (la >= 128 || lb >= 128 || ratio <= 0.34 || ratio >= 2.99) {
        _maxIntRangeRet.startI = findPredecessorIndex(a, start(b));
        _maxIntRangeRet.startJ = findPredecessorIndex(b, start(a));
        _maxIntRangeRet.endI = findPredecessorIndex(a, end(b));
        _maxIntRangeRet.endJ = findPredecessorIndex(b, end(a));
    } else {
        _maxIntRangeRet.startI = 0;
        _maxIntRangeRet.startJ = 0;
        _maxIntRangeRet.endI = la;
        _maxIntRangeRet.endJ = lb;
    }
    return _maxIntRangeRet;
}