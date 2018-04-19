/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../linear-algebra'
import { PositionData } from '../common'

interface Box3D { min: Vec3, max: Vec3 }

namespace Box3D {
    export function computeBounding(data: PositionData): Box3D {
        const min = [Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE];
        const max = [-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE];

        const { x, y, z, indices } = data;

        if (indices) {
            for (let t = 0, _t = indices.length; t < _t; t++) {
                const i = indices[t];
                min[0] = Math.min(x[i], min[0]);
                min[1] = Math.min(y[i], min[1]);
                min[2] = Math.min(z[i], min[2]);
                max[0] = Math.max(x[i], max[0]);
                max[1] = Math.max(y[i], max[1]);
                max[2] = Math.max(z[i], max[2]);
            }
        } else {
            for (let i = 0, _i = x.length; i < _i; i++) {
                min[0] = Math.min(x[i], min[0]);
                min[1] = Math.min(y[i], min[1]);
                min[2] = Math.min(z[i], min[2]);
                max[0] = Math.max(x[i], max[0]);
                max[1] = Math.max(y[i], max[1]);
                max[2] = Math.max(z[i], max[2]);
            }
        }

        return { min: Vec3.create(min[0], min[1], min[2]), max: Vec3.create(max[0], max[1], max[2]) }
    }
}

export { Box3D }