/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { OrderedSet } from '@molstar/core/data/int/ordered-set';
import type { Mat4 } from '../linear-algebra/3d/mat4.js';
import { Vec3 } from '../linear-algebra/3d/vec3.js';
import type { Tensor } from '../linear-algebra/tensor.js';
import { Box3D } from './primitives/box3d.js';

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
    resolution: number,
    maxRadius: number,
}

export interface RegularGrid3d {
    box: Box3D,
    dimensions: Vec3
}

export function getRegularGrid3dDelta({ box, dimensions }: RegularGrid3d) {
    return Vec3.div(Vec3(), Box3D.size(Vec3(), box), Vec3.subScalar(Vec3(), dimensions, 1));
}

export function fillGridDim(length: number, start: number, step: number) {
    const a = new Float32Array(length);
    for (let i = 0; i < a.length; i++) {
        a[i] = start + (step * i);
    }
    return a;
}