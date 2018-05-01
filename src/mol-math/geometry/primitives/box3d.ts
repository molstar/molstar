/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../linear-algebra'
import { PositionData } from '../common'
import { OrderedSet } from 'mol-data/int';

interface Box3D { min: Vec3, max: Vec3 }

namespace Box3D {
    export function create(min: Vec3, max: Vec3): Box3D { return { min, max }; }
    
    export function computeBounding(data: PositionData): Box3D {
        const min = [Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE];
        const max = [-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE];

        const { x, y, z, indices } = data;
        for (let t = 0, _t = OrderedSet.size(indices); t < _t; t++) {
            const i = OrderedSet.getAt(indices, t);
            min[0] = Math.min(x[i], min[0]);
            min[1] = Math.min(y[i], min[1]);
            min[2] = Math.min(z[i], min[2]);
            max[0] = Math.max(x[i], max[0]);
            max[1] = Math.max(y[i], max[1]);
            max[2] = Math.max(z[i], max[2]);
        }
        return { min: Vec3.create(min[0], min[1], min[2]), max: Vec3.create(max[0], max[1], max[2]) }
    }

    export function size(box: Box3D) {
        return Vec3.sub(Vec3.zero(), box.max, box.min);
    }

    export function expand(box: Box3D, delta: Vec3): Box3D {
        return {
            min: Vec3.sub(Vec3.zero(), box.min, delta),
            max: Vec3.add(Vec3.zero(), box.max, delta)
        }
    }
}

export { Box3D }