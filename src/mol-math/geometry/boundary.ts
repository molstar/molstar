/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PositionData } from './common';
import { Vec3 } from '../linear-algebra';
import { OrderedSet } from '../../mol-data/int';
import { BoundaryHelper } from './boundary-helper';
import { Box3D, Sphere3D } from '../geometry';

const boundaryHelperCoarse = new BoundaryHelper('14');
const boundaryHelperFine = new BoundaryHelper('98');
function getBoundaryHelper(count: number) {
    return count > 100_000 ? boundaryHelperCoarse : boundaryHelperFine;
}

export type Boundary = { readonly box: Box3D, readonly sphere: Sphere3D }

export function getBoundary(data: PositionData): Boundary {
    const { x, y, z, radius, indices } = data;
    const p = Vec3();

    const boundaryHelper = getBoundaryHelper(OrderedSet.size(indices));
    boundaryHelper.reset();
    for (let t = 0, _t = OrderedSet.size(indices); t < _t; t++) {
        const i = OrderedSet.getAt(indices, t);
        Vec3.set(p, x[i], y[i], z[i]);
        boundaryHelper.includePositionRadius(p, (radius && radius[i]) || 0);
    }
    boundaryHelper.finishedIncludeStep();
    for (let t = 0, _t = OrderedSet.size(indices); t < _t; t++) {
        const i = OrderedSet.getAt(indices, t);
        Vec3.set(p, x[i], y[i], z[i]);
        boundaryHelper.radiusPositionRadius(p, (radius && radius[i]) || 0);
    }

    const sphere = boundaryHelper.getSphere();

    if (!radius && OrderedSet.size(indices) <= 98) {
        const extrema: Vec3[] = [];
        for (let t = 0, _t = OrderedSet.size(indices); t < _t; t++) {
            const i = OrderedSet.getAt(indices, t);
            extrema.push(Vec3.create(x[i], y[i], z[i]));
        }
        Sphere3D.setExtrema(sphere, extrema);
    }

    return { box: boundaryHelper.getBox(), sphere };
}