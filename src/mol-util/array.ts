/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from './type-helpers';

// TODO move to mol-math as Vector???

/** Get the maximum value in an array */
export function arrayMax(array: NumberArray) {
    let max = -Infinity
    for (let i = 0, il = array.length; i < il; ++i) {
        if (array[i] > max) max = array[i]
    }
    return max
}

/** Get the minimum value in an array */
export function arrayMin(array: NumberArray) {
    let min = Infinity
    for (let i = 0, il = array.length; i < il; ++i) {
        if (array[i] < min) min = array[i]
    }
    return min
}

/** Get the sum of values in an array */
export function arraySum(array: NumberArray, stride = 1, offset = 0) {
    const n = array.length
    let sum = 0
    for (let i = offset; i < n; i += stride) {
        sum += array[i]
    }
    return sum
}

/** Get the mean of values in an array */
export function arrayMean(array: NumberArray, stride = 1, offset = 0) {
    return arraySum(array, stride, offset) / (array.length / stride)
}

/** Get the root mean square of values in an array */
export function arrayRms(array: NumberArray) {
    const n = array.length
    let sumSq = 0
    for (let i = 0; i < n; ++i) {
        const di = array[i]
        sumSq += di * di
    }
    return Math.sqrt(sumSq / n)
}

/** Fill an array with serial numbers starting from 0 until n - 1 (defaults to array.length) */
export function fillSerial<T extends NumberArray> (array: T, n?: number) {
    for (let i = 0, il = n ? Math.min(n, array.length) : array.length; i < il; ++i) array[ i ] = i
    return array
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
    i++;
    for (; i < l; i++) {
        xs[i] = xs[i - 1];
    }
    xs.pop();
    return true;
}
(window as any).arrayRem = arrayRemoveInPlace