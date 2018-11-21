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

function computeElementsPositionBoundary(elements: SortedArray<ElementIndex>, position: SymmetryOperator.CoordinateMapper<ElementIndex>): Boundary {
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
        if (d > radiusSq) radiusSq = d
    }

    return {
        box: { min, max },
        sphere: { center, radius: Math.sqrt(radiusSq) }
    }
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

    const { units } = s

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

    let radius = 0
    let size = 0

    for (let i = 0, _i = units.length; i < _i; i++) {
        const u = units[i]
        const invariantBoundary = getInvariantBoundary(u)
        const m = u.conformation.operator.matrix
        size += u.elements.length
        Box3D.transform(tmpBox, invariantBoundary.box, m)
        Vec3.min(min, min, tmpBox.min)
        Vec3.max(max, max, tmpBox.max)
        Sphere3D.transform(tmpSphere, invariantBoundary.sphere, m)
        Vec3.scaleAndAdd(center, center, tmpSphere.center, u.elements.length)
    }

    if (size > 0) Vec3.scale(center, center, 1/size)

    for (let i = 0, _i = units.length; i < _i; i++) {
        const u = units[i]
        const invariantBoundary = getInvariantBoundary(u)
        const m = u.conformation.operator.matrix
        Sphere3D.transform(tmpSphere, invariantBoundary.sphere, m)
        const d = Vec3.distance(tmpSphere.center, center) + tmpSphere.radius
        if (d > radius) radius = d
    }

    return { box: { min, max }, sphere: { center, radius } }
}