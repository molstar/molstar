/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../linear-algebra'
import { PositionData } from '../common'

interface Sphere3D { center: Vec3, radius: number }

namespace Sphere3D {
    export function computeBounding(data: PositionData): Sphere3D {
        const { x, y, z, indices } = data;
        let cx = 0, cy = 0, cz = 0;
        let radiusSq = 0;

        if (indices) {
            for (let t = 0, _t = indices.length; t < _t; t++) {
                const i = indices[t];
                cx += x[i];
                cy += y[i];
                cz += z[i];
            }

            if (indices.length > 0) {
                cx /= indices.length;
                cy /= indices.length;
                cz /= indices.length;
            }

            for (let t = 0, _t = indices.length; t < _t; t++) {
                const i = indices[t];
                const dx = x[i] - cx, dy = y[i] - cy, dz = z[i] - cz;
                const d = dx * dx + dy * dy + dz * dz;
                if (d > radiusSq) radiusSq = d;
            }
        } else {
            for (let i = 0, _i = x.length; i < _i; i++) {
                cx += x[i];
                cy += y[i];
                cz += z[i];
            }

            if (x.length > 0) {
                cx /= x.length;
                cy /= x.length;
                cz /= x.length;
            }

            for (let i = 0, _i = x.length; i < _i; i++) {
                const dx = x[i] - cx, dy = y[i] - cy, dz = z[i] - cz;
                const d = dx * dx + dy * dy + dz * dz;
                if (d > radiusSq) radiusSq = d;
            }
        }

        return { center: Vec3.create(cx, cy, cz), radius: Math.sqrt(radiusSq) };
    }
}

export { Sphere3D }