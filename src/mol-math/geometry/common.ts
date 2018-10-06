/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from 'mol-data/int'
import { Mat4, Tensor, Vec3 } from '../linear-algebra';
import { RenderTarget } from 'mol-gl/webgl/render-target';
import { Box3D } from '../geometry';

export interface PositionData {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,
    // subset indices into the x/y/z/radius arrays
    indices: OrderedSet,
    // optional element radius
    radius?: ArrayLike<number>
}

export type DensityData = {
    transform: Mat4,
    field: Tensor,
    idField: Tensor,

    renderTarget?: RenderTarget,
    bbox?: Box3D,
    gridDimension?: Vec3
}