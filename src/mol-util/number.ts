/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Determine the number of digits in a floating point number
 * Find a number M such that round(M * v) - M * v <= delta.
 * If no such M exists, return -1.
 */
export function getMantissaMultiplier(v: number, maxDigits: number, delta: number) {
    let m = 1, i;
    for (i = 0; i < maxDigits; i++) {
        let mv = m * v;
        if (Math.abs(Math.round(mv) - mv) <= delta) return i;
        m *= 10;
    }
    return -1;
}

export function integerDigitCount(v: number, delta: number) {
    const f = Math.abs(v);
    if (f < delta) return 0;
    return Math.floor(Math.log10(Math.abs(v))) + 1;
}

/**
 * Determine the maximum number of digits in a floating point array.
 * Find a number M such that round(M * v) - M * v <= delta.
 * If no such M exists, return -1.
 */
export function getArrayDigitCount(xs: ArrayLike<number>, maxDigits: number, delta: number) {
    let mantissaDigits = 1;
    let integerDigits = 0;
    for (let i = 0, _i = xs.length; i < _i; i++) {
        if (mantissaDigits >= 0) {
            const t = getMantissaMultiplier(xs[i], maxDigits, delta);
            if (t < 0) mantissaDigits = -1;
            else if (t > mantissaDigits) mantissaDigits = t;
        }
        const abs = Math.abs(xs[i]);
        if (abs > delta) {
            const d = Math.floor(Math.log10(Math.abs(abs))) + 1;
            if (d > integerDigits) integerDigits = d;
        }
    }
    return { mantissaDigits, integerDigits };
}

export function isInteger(s: string) {
    s = s.trim();
    const n = parseInt(s, 10);
    return isNaN(n) ? false : n.toString() === s;
}

export function getPrecision(v: number) {
    if (!isFinite(v)) return 0;
    let e = 1;
    let p = 0;
    while (Math.round(v * e) / e !== v) {
        e *= 10;
        ++p;
    }
    return p;
}

export function toPrecision(v: number, precision: number) {
    return parseFloat(v.toPrecision(precision));
}