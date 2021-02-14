/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from '../../../../mol-math/geometry';
import { Structure } from '../structure';
import { BoundaryHelper } from '../../../../mol-math/geometry/boundary-helper';
import { Boundary } from '../../../../mol-math/geometry/boundary';

const tmpSphere = Sphere3D();

const boundaryHelperCoarse = new BoundaryHelper('14');
const boundaryHelperFine = new BoundaryHelper('98');
function getBoundaryHelper(count: number) {
    return count > 500 ? boundaryHelperCoarse : boundaryHelperFine;
}

export function computeStructureBoundary(s: Structure): Boundary {
    const { units } = s;

    const boundaryHelper = getBoundaryHelper(units.length);
    boundaryHelper.reset();

    for (let i = 0, _i = units.length; i < _i; i++) {
        const u = units[i];
        const invariantBoundary = u.boundary;
        const o = u.conformation.operator;

        if (o.isIdentity) {
            boundaryHelper.includeSphere(invariantBoundary.sphere);
        } else {
            Sphere3D.transform(tmpSphere, invariantBoundary.sphere, o.matrix);
            boundaryHelper.includeSphere(tmpSphere);
        }
    }

    boundaryHelper.finishedIncludeStep();

    for (let i = 0, _i = units.length; i < _i; i++) {
        const u = units[i];
        const invariantBoundary = u.boundary;
        const o = u.conformation.operator;

        if (o.isIdentity) {
            boundaryHelper.radiusSphere(invariantBoundary.sphere);
        } else {
            Sphere3D.transform(tmpSphere, invariantBoundary.sphere, o.matrix);
            boundaryHelper.radiusSphere(tmpSphere);
        }
    }

    return { box: boundaryHelper.getBox(), sphere: boundaryHelper.getSphere() };
}