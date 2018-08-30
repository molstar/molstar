/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Structure from '../structure'
import Unit from '../unit';
import { Box3D, Sphere3D, SymmetryOperator } from 'mol-math/geometry';
import { Vec3 } from 'mol-math/linear-algebra';
import { SortedArray } from 'mol-data/int';
import { ElementIndex } from '../../model/indexing';

export type Boundary = { box: Box3D, sphere: Sphere3D }

function computeElementsPositionBoundary(elements: SortedArray<ElementIndex>, position: SymmetryOperator.CoordinateMapper): Boundary {
    const min = Vec3.create(Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE)
    const max = Vec3.create(-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE)
    const center = Vec3.zero()

    let radiusSq = 0
    let size = 0

    const p = Vec3.zero()

    size += elements.length
    for (let j = 0, _j = elements.length; j < _j; j++) {
        position(elements[j], p)
        Vec3.min(min, min, p)
        Vec3.max(max, max, p)
        Vec3.add(center, center, p)
    }

    if (size > 0) Vec3.scale(center, center, 1/size)

    for (let j = 0, _j = elements.length; j < _j; j++) {
        position(elements[j], p)
        const d = Vec3.squaredDistance(p, center)
        if (d > radiusSq) radiusSq = d;
    }

    return {
        box: { min, max },
        sphere: { center, radius: Math.sqrt(radiusSq) }
    };
}

function computeInvariantUnitBoundary(u: Unit): Boundary {
    return computeElementsPositionBoundary(u.elements, u.conformation.invariantPosition)
}

export function computeUnitBoundary(u: Unit): Boundary {
    return computeElementsPositionBoundary(u.elements, u.conformation.position)
}

const tmpBox = Box3D.empty()
const tmpSphere = Sphere3D.zero()

export function computeStructureBoundary(s: Structure): Boundary {
    const min = Vec3.create(Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE)
    const max = Vec3.create(-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE)
    const center = Vec3.zero()

    const { units } = s;

    const boundaryMap: Map<number, Boundary> = new Map()
    function getInvariantBoundary(u: Unit) {
        let boundary: Boundary
        if (boundaryMap.has(u.invariantId)) {
            boundary = boundaryMap.get(u.invariantId)!
        } else {
            boundary = computeInvariantUnitBoundary(u)
            boundaryMap.set(u.invariantId, boundary)
        }
        return boundary
    }

    let radiusSq = 0;
    let size = 0;

    for (let i = 0, _i = units.length; i < _i; i++) {
        size += 1
        const u = units[i]
        const invariantBoundary = getInvariantBoundary(u)
        const m = u.conformation.operator.matrix
        Box3D.transform(tmpBox, invariantBoundary.box, m)
        Vec3.min(min, min, tmpBox.min)
        Vec3.max(max, max, tmpBox.max)
        Sphere3D.transform(tmpSphere, invariantBoundary.sphere, m)
        Vec3.add(center, center, tmpSphere.center)
    }

    if (size > 0) Vec3.scale(center, center, 1/size)

    for (let i = 0, _i = units.length; i < _i; i++) {
        const u = units[i]
        const invariantBoundary = getInvariantBoundary(u)
        const m = u.conformation.operator.matrix
        Sphere3D.transform(tmpSphere, invariantBoundary.sphere, m)
        const d = Vec3.squaredDistance(tmpSphere.center, center) + (tmpSphere.radius * tmpSphere.radius) * 4
        if (d > radiusSq) radiusSq = d
    }

    const b = {
        box: { min, max },
        sphere: { center, radius: Math.sqrt(radiusSq) }
    };
    console.log(b, computeStructureBoundary2(s))
    return b
}

export function computeStructureBoundary2(s: Structure): Boundary {
    const min = [Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE];
    const max = [-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE];

    const { units } = s;

    let cx = 0, cy = 0, cz = 0;
    let radiusSq = 0;
    let size = 0;

    for (let i = 0, _i = units.length; i < _i; i++) {
        const { x, y, z } = units[i].conformation;

        const elements = units[i].elements;
        size += elements.length;
        for (let j = 0, _j = elements.length; j < _j; j++) {
            const e = elements[j];
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

    for (let i = 0, _i = units.length; i < _i; i++) {
        const { x, y, z } = units[i].conformation;

        const elements = units[i].elements;
        for (let j = 0, _j = elements.length; j < _j; j++) {
            const e = elements[j];
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