/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from '../../../../mol-math/linear-algebra';
import { Box3D, Axes3D } from '../../../../mol-math/geometry';
import { MeshBuilder } from '../mesh-builder';
import { CylinderProps } from '../../../primitive/cylinder';
import { addCylinder } from './cylinder';
import { addSphere } from './sphere';
import { createCage } from '../../../primitive/cage';

const tmpStart = Vec3.zero();
const tmpEnd = Vec3.zero();
const cylinderProps: CylinderProps = {};

export function addBoundingBox(state: MeshBuilder.State, box: Box3D, radius: number, detail: number, radialSegments: number) {
    const { min, max } = box;

    cylinderProps.radiusTop = radius;
    cylinderProps.radiusBottom = radius;
    cylinderProps.radialSegments = radialSegments;

    Vec3.set(tmpStart, max[0], max[1], max[2]);
    addSphere(state, tmpStart, radius, detail);
    Vec3.set(tmpEnd, max[0], max[1], min[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);
    Vec3.set(tmpEnd, max[0], min[1], max[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);
    Vec3.set(tmpEnd, min[0], max[1], max[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);

    Vec3.set(tmpStart, min[0], min[1], min[2]);
    addSphere(state, tmpStart, radius, detail);
    Vec3.set(tmpEnd, min[0], min[1], max[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);
    Vec3.set(tmpEnd, min[0], max[1], min[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);
    Vec3.set(tmpEnd, max[0], min[1], min[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);

    Vec3.set(tmpStart, max[0], min[1], min[2]);
    addSphere(state, tmpStart, radius, detail);
    Vec3.set(tmpEnd, max[0], min[1], max[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);
    Vec3.set(tmpEnd, max[0], max[1], min[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);

    Vec3.set(tmpStart, min[0], min[1], max[2]);
    addSphere(state, tmpStart, radius, detail);
    Vec3.set(tmpEnd, min[0], max[1], max[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);
    Vec3.set(tmpEnd, max[0], min[1], max[2]);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);

    Vec3.set(tmpStart, min[0], max[1], min[2]);
    addSphere(state, tmpStart, radius, detail);
    Vec3.set(tmpEnd, max[0], max[1], min[2]);
    addSphere(state, tmpEnd, radius, detail);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);
    Vec3.set(tmpEnd, min[0], max[1], max[2]);
    addSphere(state, tmpEnd, radius, detail);
    addCylinder(state, tmpStart, tmpEnd, 1, cylinderProps);
}

//

const tmpBoxVecCorner = Vec3();
const tmpBoxVecA = Vec3();
const tmpBoxVecB = Vec3();
const tmpBoxVecC = Vec3();
const tmpMatrix = Mat4.identity();

const tmpVertices = new Float32Array(8 * 3);
const tmpEdges = new Uint8Array([
    0, 1, 0, 3, 0, 6, 1, 2, 1, 7, 2, 3,
    2, 4, 3, 5, 4, 5, 4, 7, 5, 6, 6, 7
]);

export function addOrientedBox(state: MeshBuilder.State, axes: Axes3D, radiusScale: number, detail: number, radialSegments: number) {
    const { origin, dirA, dirB, dirC } = axes;
    const negDirA = Vec3.negate(tmpBoxVecA, dirA);
    const negDirB = Vec3.negate(tmpBoxVecB, dirB);
    const negDirC = Vec3.negate(tmpBoxVecC, dirC);

    let offset = 0;
    const addCornerHelper = function (v1: Vec3, v2: Vec3, v3: Vec3) {
        Vec3.copy(tmpBoxVecCorner, origin);
        Vec3.add(tmpBoxVecCorner, tmpBoxVecCorner, v1);
        Vec3.add(tmpBoxVecCorner, tmpBoxVecCorner, v2);
        Vec3.add(tmpBoxVecCorner, tmpBoxVecCorner, v3);
        Vec3.toArray(tmpBoxVecCorner, tmpVertices, offset);
        offset += 3;
    };
    addCornerHelper(dirA, dirB, dirC);
    addCornerHelper(dirA, dirB, negDirC);
    addCornerHelper(dirA, negDirB, negDirC);
    addCornerHelper(dirA, negDirB, dirC);
    addCornerHelper(negDirA, negDirB, negDirC);
    addCornerHelper(negDirA, negDirB, dirC);
    addCornerHelper(negDirA, dirB, dirC);
    addCornerHelper(negDirA, dirB, negDirC);

    const cage = createCage(tmpVertices, tmpEdges);
    const volume = Axes3D.volume(axes);
    const radius = (Math.cbrt(volume) / 300) * radiusScale;

    MeshBuilder.addCage(state, tmpMatrix, cage, radius, detail, radialSegments);
}