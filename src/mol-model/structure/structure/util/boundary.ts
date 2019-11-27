/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D, Sphere3D } from '../../../../mol-math/geometry';
import { BoundaryHelper } from '../../../../mol-math/geometry/boundary-helper';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import Structure from '../structure';
import { CentroidHelper } from '../../../../mol-math/geometry/centroid-helper';

export type Boundary = { box: Box3D, sphere: Sphere3D }

const tmpBox = Box3D.empty();
const tmpSphere = Sphere3D.zero();
const tmpVec = Vec3.zero();

const boundaryHelper = new BoundaryHelper();
const centroidHelper = new CentroidHelper();

// TODO: show this be enabled? does not solve problem with "renderable bounding spheres"
const CENTROID_COMPUTATION_THRESHOLD = -1;

export function computeStructureBoundary(s: Structure): Boundary {
    if (s.elementCount <= CENTROID_COMPUTATION_THRESHOLD && s.isAtomic) return computeStructureBoundaryFromElements(s);
    return computeStructureBoundaryFromUnits(s);
}

export function computeStructureBoundaryFromElements(s: Structure): Boundary {
    centroidHelper.reset();

    const min = Vec3.create(Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE)
    const max = Vec3.create(-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE)

    const v = tmpVec;
    for (const unit of s.units) {
        const { x, y, z } = unit.conformation;
        const elements = unit.elements;
        for (let j = 0, _j = elements.length; j < _j; j++) {
            const e = elements[j];
            Vec3.set(v, x(e), y(e), z(e));
            centroidHelper.includeStep(v);
            Vec3.min(min, min, v);
            Vec3.max(max, max, v);
        };
    }

    centroidHelper.finishedIncludeStep();

    for (const unit of s.units) {
        const { x, y, z } = unit.conformation;
        const elements = unit.elements;
        for (let j = 0, _j = elements.length; j < _j; j++) {
            const e = elements[j];
            Vec3.set(v, x(e), y(e), z(e));
            centroidHelper.radiusStep(v);
        };
    }

    return { box: { min, max }, sphere: centroidHelper.getSphere() };
}

export function computeStructureBoundaryFromUnits(s: Structure): Boundary {
    const min = Vec3.create(Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE)
    const max = Vec3.create(-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE)

    const { units } = s;

    boundaryHelper.reset(0);

    for (let i = 0, _i = units.length; i < _i; i++) {
        const u = units[i];
        const invariantBoundary = u.lookup3d.boundary;
        const o = u.conformation.operator;

        if (o.isIdentity) {
            Vec3.min(min, min, invariantBoundary.box.min);
            Vec3.max(max, max, invariantBoundary.box.max);

            boundaryHelper.boundaryStep(invariantBoundary.sphere.center, invariantBoundary.sphere.radius);
        } else {
            Box3D.transform(tmpBox, invariantBoundary.box, o.matrix);
            Vec3.min(min, min, tmpBox.min);
            Vec3.max(max, max, tmpBox.max);

            Sphere3D.transform(tmpSphere, invariantBoundary.sphere, o.matrix);
            boundaryHelper.boundaryStep(tmpSphere.center, tmpSphere.radius);
        }
    }

    boundaryHelper.finishBoundaryStep();

    for (let i = 0, _i = units.length; i < _i; i++) {
        const u = units[i];
        const invariantBoundary = u.lookup3d.boundary;
        const o = u.conformation.operator;

        if (o.isIdentity) {
            boundaryHelper.extendStep(invariantBoundary.sphere.center, invariantBoundary.sphere.radius);
        } else {
            Sphere3D.transform(tmpSphere, invariantBoundary.sphere, o.matrix);
            boundaryHelper.extendStep(tmpSphere.center, tmpSphere.radius);
        }
    }

    return { box: { min, max }, sphere: boundaryHelper.getSphere() };
}