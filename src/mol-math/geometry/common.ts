/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from 'mol-data/int'
import { Mat4, Tensor, Vec3 } from '../linear-algebra';
import { Box3D } from '../geometry';
import { Texture } from 'mol-gl/webgl/texture';

export interface PositionData {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,
    /** subset of indices into the x/y/z/radius arrays */
    indices: OrderedSet,
    /** optional element radius */
    radius?: ArrayLike<number>
}

export type DensityData = {
    transform: Mat4,
    field: Tensor,
    idField: Tensor,
}

export type DensityTextureData = {
    transform: Mat4,
    texture: Texture,
    bbox: Box3D,
    gridDimension: Vec3
}