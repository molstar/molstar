/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Jsonable, canonicalJsonString } from './json';
import { NumberArray } from './type-helpers';

// TODO move to mol-math as Vector???

/** Get the maximum value in an array */
export function arrayMax(array: ArrayLike<number>) {
    let max = -Infinity;
    for (let i = 0, il = array.length; i < il; ++i) {
        if (array[i] > max) max = array[i];
    }
    return max;
}

/** Get the minimum value in an array */
export function arrayMin(array: ArrayLike<number>) {
    let min = Infinity;
    for (let i = 0, il = array.length; i < il; ++i) {
        if (array[i] < min) min = array[i];
    }
    return min;
}

/** Get the minimum & maximum value in an array */
export function arrayMinMax(array: ArrayLike<number>) {
    let min = Infinity;
    let max = -Infinity;
    for (let i = 0, il = array.length; i < il; ++i) {
        if (array[i] < min) min = array[i];
        if (array[i] > max) max = array[i];
    }
    return [min, max];
}

/** Get the sum of values in an array */
export function arraySum(array: ArrayLike<number>, stride = 1, offset = 0) {
    const n = array.length;
    let sum = 0;
    for (let i = offset; i < n; i += stride) {
        sum += array[i];
    }
    return sum;
}

/** Get the mean of values in an array */
export function arrayMean(array: ArrayLike<number>, stride = 1, offset = 0) {
    return arraySum(array, stride, offset) / (array.length / stride);
}

/** Get the root mean square of values in an array */
export function arrayRms(array: ArrayLike<number>) {
    const n = array.length;
    let sumSq = 0;
    for (let i = 0; i < n; ++i) {
        const di = array[i];
        sumSq += di * di;
    }
    return Math.sqrt(sumSq / n);
}

/** Fill an array with serial numbers starting from 0 until n - 1 (defaults to array.length) */
export function fillSerial<T extends NumberArray>(array: T, n?: number) {
    for (let i = 0, il = n ? Math.min(n, array.length) : array.length; i < il; ++i) array[i] = i;
    return array;
}

export function arrayRemoveInPlace<T>(xs: T[], x: T) {
    let i = 0, found = false;
    for (const il = xs.length; i < il; i++) {
        if (xs[i] === x) {
            found = true;
            break;
        }
    }
    if (!found) return false;
    arrayRemoveAtInPlace(xs, i);
    return true;
}

export function arrayRemoveAtInPlace<T>(xs: T[], idx: number) {
    for (let i = idx, _i = xs.length - 1; i < _i; i++) {
        xs[i] = xs[i + 1];
    }
    xs.pop();
}

export function arraySetAdd<T>(xs: T[], x: T) {
    if (xs.indexOf(x) >= 0) return false;
    xs.push(x);
    return true;
}

export function arraySetRemove<T>(xs: T[], x: T) {
    const idx = xs.indexOf(x);
    if (idx < 0) return false;
    for (let i = idx, _i = xs.length - 1; i < _i; i++) {
        xs[i] = xs[i + 1];
    }
    xs.pop();
    return true;
}

/**
 * Caution, O(n^2) complexity. Only use for small input sizes.
 * For larger inputs consider using `SortedArray`.
 */
export function arrayAreIntersecting<T>(xs: T[], ys: T[]) {
    for (let i = 0, il = xs.length; i < il; ++i) {
        if (ys.includes(xs[i])) return true;
    }
    return false;
}

/**
 * Caution, O(n^2) complexity. Only use for small input sizes.
 * For larger inputs consider using `SortedArray`.
 */
export function arrayIntersectionSize<T>(xs: T[], ys: T[]) {
    let count = 0;
    for (let i = 0, il = xs.length; i < il; ++i) {
        if (ys.includes(xs[i])) count += 1;
    }
    return count;
}

export function arrayEqual<T>(xs?: ArrayLike<T>, ys?: ArrayLike<T>) {
    if (!xs || xs.length === 0) return !ys || ys.length === 0;
    if (!ys) return false;

    const lenX = xs.length;
    if (lenX !== ys.length) return false;
    for (let i = 0; i < lenX; i++) {
        if (xs[i] !== ys[i]) return false;
    }
    return true;
}

export function arrayIsIdentity(xs: ArrayLike<number>) {
    for (let i = 0, _i = xs.length; i < _i; i++) {
        if (xs[i] !== i) return false;
    }
    return true;
}

export function arrayMapUpsert<T>(xs: [string, T][], key: string, value: T) {
    for (let i = 0, il = xs.length; i < il; ++i) {
        if (xs[i][0] === key) {
            xs[i][1] = value;
            return;
        }
    }
    xs.push([key, value]);
}

/** Return an array containing integers from [start, end) if `end` is given,
 * or from [0, start) if `end` is omitted. */
export function range(start: number, end?: number): number[] {
    if (end === undefined) {
        end = start;
        start = 0;
    }
    const length = Math.max(end - start, 0);
    const result = Array(length);
    for (let i = 0; i < length; i++) {
        result[i] = start + i;
    }
    return result;
}

/** Copy all elements from `src` to the end of `dst`.
 * Equivalent to `dst.push(...src)`, but avoids storing element on call stack. Faster that `extend` from Underscore.js.
 * `extend(a, a)` will double the array.
 * Returns the modified `dst` array.
 */
export function arrayExtend<T>(dst: T[], src: ArrayLike<T>): T[] {
    const offset = dst.length;
    const nCopy = src.length;
    dst.length += nCopy;
    for (let i = 0; i < nCopy; i++) {
        dst[offset + i] = src[i];
    }
    return dst;
}

/** Check whether `array` is sorted, sort if not. */
export function sortIfNeeded<T>(array: T[], compareFn: (a: T, b: T) => number): T[] {
    return arrayIsSorted(array, compareFn) ? array : array.sort(compareFn);
}

/** Decide whether `array` is sorted. */
export function arrayIsSorted<T>(array: T[], compareFn: (a: T, b: T) => number): boolean {
    for (let i = 1, n = array.length; i < n; i++) {
        if (compareFn(array[i - 1], array[i]) > 0) {
            return false;
        }
    }
    return true;
}

/** Remove all elements from the array which do not fulfil `predicate`. Return the modified array itself. */
export function filterInPlace<T>(array: T[], predicate: (x: T) => boolean): T[] {
    const n = array.length;
    let iDest = 0;
    for (let iSrc = 0; iSrc < n; iSrc++) {
        if (predicate(array[iSrc])) {
            array[iDest++] = array[iSrc];
        }
    }
    array.length = iDest;
    return array;
}

/** Return an array of all distinct values from `values`
 * (i.e. with removed duplicates).
 * Uses deep equality for objects and arrays,
 * independent from object key order and undefined properties.
 * E.g. {a: 1, b: undefined, c: {d: [], e: null}} is equal to {c: {e: null, d: []}}, a: 1}.
 * If two or more objects in `values` are equal, only the first of them will be in the result. */
export function arrayDistinct<T extends Jsonable>(values: T[]): T[] {
    const seen = new Set<string>();
    const result: T[] = [];
    for (const value of values) {
        const key = canonicalJsonString(value);
        if (!seen.has(key)) {
            seen.add(key);
            result.push(value);
        }
    }
    return result;
}
