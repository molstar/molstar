/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4 } from '../../../../mol-math/linear-algebra';
import { MeshBuilder } from '../mesh-builder';
import { Primitive } from '../../../primitive/primitive';
import { Sphere } from '../../../primitive/sphere';

const sphereMap = new Map<number, Primitive>();
const tmpSphereMat = Mat4.identity();

function setSphereMat(m: Mat4, center: Vec3, radius: number) {
    return Mat4.scaleUniformly(m, Mat4.fromTranslation(m, center), radius);
}

export function getSphere(detail: number) {
    let sphere = sphereMap.get(detail);
    if (sphere === undefined) {
        sphere = Sphere(detail);
        sphereMap.set(detail, sphere);
    }
    return sphere;
}

export function addSphere(state: MeshBuilder.State, center: Vec3, radius: number, detail: number) {
    MeshBuilder.addPrimitive(state, setSphereMat(tmpSphereMat, center, radius), getSphere(detail));
}