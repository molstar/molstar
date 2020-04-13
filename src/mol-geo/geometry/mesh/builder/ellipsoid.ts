/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from '../../../../mol-math/linear-algebra';
import { MeshBuilder } from '../mesh-builder';
import { getSphere } from './sphere';

const tmpEllipsoidMat = Mat4.identity();
const tmpVec = Vec3();

function setEllipsoidMat(m: Mat4, center: Vec3, dirMajor: Vec3, dirMinor: Vec3, radiusScale: Vec3) {
    Vec3.add(tmpVec, center, dirMajor);
    Mat4.targetTo(m, center, tmpVec, dirMinor);
    Mat4.setTranslation(m, center);
    return Mat4.scale(m, m, radiusScale);
}

export function addEllipsoid(state: MeshBuilder.State, center: Vec3, dirMajor: Vec3, dirMinor: Vec3, radiusScale: Vec3, detail: number) {
    MeshBuilder.addPrimitive(state, setEllipsoidMat(tmpEllipsoidMat, center, dirMajor, dirMinor, radiusScale), getSphere(detail));
}