/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../mol-data/int';
import { Mat4, Tensor, Vec3, Vec2 } from '../linear-algebra';
import { Box3D } from '../geometry';
import { Texture } from '../../mol-gl/webgl/texture';

export interface PositionData {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,
    /** subset of indices into the x/y/z/radius arrays */
    indices: OrderedSet,
    /** optional element radius */
    radius?: ArrayLike<number>,
    /** optional element id */
    id?: ArrayLike<number>,
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
    gridDim: Vec3,
    gridTexDim: Vec3
    gridTexScale: Vec2
}

export function fillGridDim(length: number, start: number, step: number) {
    const a = new Float32Array(length);
    for (let i = 0; i < a.length; i++) {
        a[i] = start + (step * i);
    }
    return a;
}