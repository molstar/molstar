/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Mat4, Vec3 } from '../../../../mol-math/linear-algebra';
import { Primitive } from '../../../primitive/primitive';
import { Sphere } from '../../../primitive/sphere';
import { MeshBuilder } from '../mesh-builder';

/** Full sphere stored on key `detail`; sphere subsets (ring and caps) stored on keys `-2*detail - 1` and `-2*detail - 2` */
const sphereCache = new Map<number, Primitive>();
const tmpSphereMat = Mat4.identity();

function setSphereMat(m: Mat4, center: Vec3, radius: number) {
    return Mat4.scaleUniformly(m, Mat4.fromTranslation(m, center), radius);
}

export function getSphere(detail: number) {
    let sphere = sphereCache.get(detail);
    if (sphere === undefined) {
        sphere = Sphere(detail);
        sphereCache.set(detail, sphere);
    }
    return sphere;
}

export function getSphereSubset(detail: number, subset: 'ring' | 'caps') {
    const cacheKey = -2 * detail + (subset === 'ring' ? -1 : -2);
    let sphere = sphereCache.get(cacheKey);
    if (sphere === undefined) {
        sphere = Sphere(detail, { subset });
        sphereCache.set(cacheKey, sphere);
    }
    return sphere;
}

export function addSphere(state: MeshBuilder.State, center: Vec3, radius: number, detail: number) {
    MeshBuilder.addPrimitive(state, setSphereMat(tmpSphereMat, center, radius), getSphere(detail));
}

export function addSphereSubset(state: MeshBuilder.State, center: Vec3, radius: number, detail: number, subset: 'ring' | 'caps') {
    MeshBuilder.addPrimitive(state, setSphereMat(tmpSphereMat, center, radius), getSphereSubset(detail, subset));
}
