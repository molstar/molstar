/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from 'mol-data/int'
import { Mat4, Tensor } from '../linear-algebra';

export interface PositionData {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,
    // subset indices into the x/y/z/radius arrays
    indices: OrderedSet,
    // optional element radius
    radius?: ArrayLike<number>
}

export type DensityData = { transform: Mat4, field: Tensor, idField: Tensor }