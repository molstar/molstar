/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Structure from '../structure'
import ElementSet from '../element/set'
import { ElementGroup } from '../../structure'
import { Box3D, Sphere3D } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';

function computeStructureBoundary(s: Structure): { box: Box3D, sphere: Sphere3D } {
    const min = [Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE];
    const max = [-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE];

    const { units, elements } = s;

    let cx = 0, cy = 0, cz = 0;
    let radiusSq = 0;
    let size = 0;

    for (let i = 0, _i = ElementSet.unitCount(elements); i < _i; i++) {
        const group = ElementSet.unitGetByIndex(elements, i);
        const { x, y, z } = units[ElementSet.unitGetId(elements, i)];

        size += ElementGroup.size(group);
        for (let j = 0, _j = ElementGroup.size(group); j < _j; j++) {
            const e = ElementGroup.getAt(group, j);
            const xx = x(e), yy = y(e), zz = z(e);

            min[0] = Math.min(xx, min[0]);
            min[1] = Math.min(yy, min[1]);
            min[2] = Math.min(zz, min[2]);
            max[0] = Math.max(xx, max[0]);
            max[1] = Math.max(yy, max[1]);
            max[2] = Math.max(zz, max[2]);

            cx += xx;
            cy += yy;
            cz += zz;
        }
    }

    if (size > 0) {
        cx /= size;
        cy /= size;
        cz /= size;
    }

    for (let i = 0, _i = ElementSet.unitCount(elements); i < _i; i++) {
        const group = ElementSet.unitGetByIndex(elements, i);
        const { x, y, z } = units[ElementSet.unitGetId(elements, i)];

        size += ElementGroup.size(group);
        for (let j = 0, _j = ElementGroup.size(group); j < _j; j++) {
            const e = ElementGroup.getAt(group, j);
            const dx = x(e) - cx, dy = y(e) - cy, dz = z(e) - cz;
            const d = dx * dx + dy * dy + dz * dz;
            if (d > radiusSq) radiusSq = d;
        }
    }

    return {
        box: { min: Vec3.ofArray(min), max: Vec3.ofArray(max) },
        sphere: { center: Vec3.create(cx, cy, cz), radius: Math.sqrt(radiusSq) }
    };
}

export { computeStructureBoundary }