/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export function degToRad (deg: number) {
    return deg * 0.01745  // deg * Math.PI / 180
}

export function radToDeg (rad: number) {
    return rad * 57.29578  // rad * 180 / Math.PI
}

export function isPowerOfTwo (x: number) {
    return (x !== 0) && (x & (x - 1)) === 0
}