/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const halfPI = Math.PI / 2;

export function degToRad (deg: number) {
    return deg * 0.01745;  // deg * Math.PI / 180
}

export function radToDeg (rad: number) {
    return rad * 57.29578;  // rad * 180 / Math.PI
}

export function isPowerOfTwo (x: number) {
    return (x !== 0) && (x & (x - 1)) === 0;
}

/** return the value that has the largest absolute value */
export function absMax(...values: number[]) {
    let max = 0;
    let absMax = 0;
    for (let i = 0, il = values.length; i < il; ++i) {
        const value = values[i];
        const abs = Math.abs(value);
        if (abs > absMax) {
            max = value;
            absMax = abs;
        }
    }
    return max;
}

/** Length of an arc with angle in radians */
export function arcLength(angle: number, radius: number) {
    return angle * radius;
}