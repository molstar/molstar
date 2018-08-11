/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

/**
 * Determine the number of digits in a floating point number
 * Find a number M such that round(M * v) - M * v <= delta.
 * If no such M exists, return -1.
 */
export function getMantissaMultiplier(v: number, maxDigits: number, delta: number) {
    let m = 1;
    for (let i = 0; i < maxDigits; i++) {
        let mv = m * v;
        if (Math.abs(Math.round(mv) - mv) <= delta) return m;
        m *= 10;
    }
    return -1;
}

/**
 * Determine the maximum number of digits in a floating point array.
 * Find a number M such that round(M * v) - M * v <= delta.
 * If no such M exists, return -1.
 */
export function getArrayMantissaMultiplier(xs: ArrayLike<number>, maxDigits: number, delta: number) {
    let m = 1;
    for (let i = 0, _i = xs.length; i < _i; i++) {
        const t = getMantissaMultiplier(xs[i], maxDigits, delta);
        if (t < 0) return -1;
        if (t > m) m = t;
    }
    return m;
}