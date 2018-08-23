/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// TODO move to mol-math as Vector???

/** Get the maximum value in an array */
export function arrayMax(array: Helpers.NumberArray) {
    let max = -Infinity
    for (let i = 0, il = array.length; i < il; ++i) {
        if (array[i] > max) max = array[i]
    }
    return max
}

/** Get the minimum value in an array */
export function arrayMin(array: Helpers.NumberArray) {
    let min = Infinity
    for (let i = 0, il = array.length; i < il; ++i) {
        if (array[i] < min) min = array[i]
    }
    return min
}

/** Get the sum of values in an array */
export function arraySum(array: Helpers.NumberArray, stride = 1, offset = 0) {
    const n = array.length
    let sum = 0
    for (let i = offset; i < n; i += stride) {
        sum += array[i]
    }
    return sum
}

/** Get the mean of values in an array */
export function arrayMean(array: Helpers.NumberArray, stride = 1, offset = 0) {
    return arraySum(array, stride, offset) / (array.length / stride)
}

/** Get the root mean square of values in an array */
export function arrayRms(array: Helpers.NumberArray) {
    const n = array.length
    let sumSq = 0
    for (let i = 0; i < n; ++i) {
        const di = array[i]
        sumSq += di * di
    }
    return Math.sqrt(sumSq / n)
}

/** Fill an array with serial numbers starting from 0 */
export function fillSerial<T extends Helpers.NumberArray> (array: T) {
    const n = array.length
    for (let i = 0; i < n; ++i) array[ i ] = i
    return array
}