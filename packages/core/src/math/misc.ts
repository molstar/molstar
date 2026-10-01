/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const halfPI = Math.PI / 2;
export const PiDiv180 = Math.PI / 180;

export function degToRad(deg: number) {
    return deg * PiDiv180; // deg * Math.PI / 180
}

export function radToDeg(rad: number) {
    return rad / PiDiv180; // rad * 180 / Math.PI
}

export function isPowerOfTwo(x: number) {
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

/** Create an outward spiral of given `radius` on a 2d grid */
export function spiral2d(radius: number) {
    let x = 0;
    let y = 0;
    let deltaX = 0;
    let deltaY = -1;
    const size = radius * 2 + 1;
    const halfSize = size / 2;
    const out: [number, number][] = [];

    for (let i = Math.pow(size, 2); i > 0; --i) {
        if ((-halfSize < x && x <= halfSize) && (-halfSize < y && y <= halfSize)) {
            out.push([x, y]);
        }

        if (x === y || (x < 0 && x === -y) || (x > 0 && x === 1 - y)) {
            // change direction
            const prevDeltaX = deltaX;
            const prevDeltaY = deltaY;
            deltaX = -prevDeltaY;
            deltaY = prevDeltaX;
        }

        x += deltaX;
        y += deltaY;
    }
    return out;
}