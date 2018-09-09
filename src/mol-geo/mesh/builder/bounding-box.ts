/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { Box3D } from 'mol-math/geometry';
import { MeshBuilder } from '../mesh-builder';
import { CylinderProps } from '../../primitive/cylinder';
import { addCylinder } from './cylinder';
import { addSphere } from './sphere';

const tmpStart = Vec3.zero()
const tmpEnd = Vec3.zero()
const cylinderProps: CylinderProps = {}

export function addBoundingBox(builder: MeshBuilder, box: Box3D, radius: number, detail: number, radialSegments: number) {
    const { min, max } = box

    cylinderProps.radiusTop = radius
    cylinderProps.radiusBottom = radius
    cylinderProps.radialSegments = radialSegments

    Vec3.set(tmpStart, max[0], max[1], max[2])
    addSphere(builder, tmpStart, radius, detail)
    Vec3.set(tmpEnd, max[0], max[1], min[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)
    Vec3.set(tmpEnd, max[0], min[1], max[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)
    Vec3.set(tmpEnd, min[0], max[1], max[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)

    Vec3.set(tmpStart, min[0], min[1], min[2])
    addSphere(builder, tmpStart, radius, detail)
    Vec3.set(tmpEnd, min[0], min[1], max[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)
    Vec3.set(tmpEnd, min[0], max[1], min[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)
    Vec3.set(tmpEnd, max[0], min[1], min[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)

    Vec3.set(tmpStart, max[0], min[1], min[2])
    addSphere(builder, tmpStart, radius, detail)
    Vec3.set(tmpEnd, max[0], min[1], max[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)
    Vec3.set(tmpEnd, max[0], max[1], min[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)

    Vec3.set(tmpStart, min[0], min[1], max[2])
    addSphere(builder, tmpStart, radius, detail)
    Vec3.set(tmpEnd, min[0], max[1], max[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)
    Vec3.set(tmpEnd, max[0], min[1], max[2])
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)

    Vec3.set(tmpStart, min[0], max[1], min[2])
    addSphere(builder, tmpStart, radius, detail)
    Vec3.set(tmpEnd, max[0], max[1], min[2])
    addSphere(builder, tmpEnd, radius, detail)
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)
    Vec3.set(tmpEnd, min[0], max[1], max[2])
    addSphere(builder, tmpEnd, radius, detail)
    addCylinder(builder, tmpStart, tmpEnd, 1, cylinderProps)
}