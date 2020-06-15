/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

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
export function fillSerial<T extends NumberArray> (array: T, n?: number) {
    for (let i = 0, il = n ? Math.min(n, array.length) : array.length; i < il; ++i) array[i] = i;
    return array;
}

export function arrayRemoveInPlace<T>(xs: T[], x: T) {
    let i = 0, l = xs.length, found = false;
    for (; i < l; i++) {
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