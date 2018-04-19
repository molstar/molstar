/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export interface PositionData {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,
    radius?: ArrayLike<number>,
    // subset indices into the x/y/z/radius arrays
    indices?: ArrayLike<number>,
}
