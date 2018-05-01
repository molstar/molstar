/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../linear-algebra'
import { PositionData } from '../common'
import { OrderedSet } from 'mol-data/int';

interface Sphere3D { center: Vec3, radius: number }

namespace Sphere3D {
    export function create(center: Vec3, radius: number): Sphere3D {
        return { center, radius };
    }

    export function computeBounding(data: PositionData): Sphere3D {
        const { x, y, z, indices } = data;
        let cx = 0, cy = 0, cz = 0;
        let radiusSq = 0;

        const size = OrderedSet.size(indices);
        for (let t = 0; t < size; t++) {
            const i = OrderedSet.getAt(indices, t);
            cx += x[i];
            cy += y[i];
            cz += z[i];
        }

        if (size > 0) {
            cx /= size;
            cy /= size;
            cz /= size;
        }

        for (let t = 0; t < size; t++) {
            const i = OrderedSet.getAt(indices, t);
            const dx = x[i] - cx, dy = y[i] - cy, dz = z[i] - cz;
            const d = dx * dx + dy * dy + dz * dz;
            if (d > radiusSq) radiusSq = d;
        }

        return { center: Vec3.create(cx, cy, cz), radius: Math.sqrt(radiusSq) };
    }
}

export { Sphere3D }