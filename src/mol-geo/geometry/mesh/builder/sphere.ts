/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Vec3, Mat4 } from '../../../../mol-math/linear-algebra';
import { MeshBuilder } from '../mesh-builder';
import { Primitive } from '../../../primitive/primitive';
import { Sphere } from '../../../primitive/sphere';

const sphereCache = new Map<string, Primitive>();
const tmpSphereMat = Mat4.identity();

function setSphereMat(m: Mat4, center: Vec3, radius: number) {
    return Mat4.scaleUniformly(m, Mat4.fromTranslation(m, center), radius);
}

export function getSphere(detail: number, options?: { subset: 'ring' | 'caps' | undefined }) {
    const cacheKey = `${detail}/${options?.subset ?? ''}`;
    let sphere = sphereCache.get(cacheKey);
    if (sphere === undefined) {
        sphere = Sphere(detail, options);
        sphereCache.set(cacheKey, sphere);
    }
    return sphere;
}

export function addSphere(state: MeshBuilder.State, center: Vec3, radius: number, detail: number, options?: { subset: 'ring' | 'caps' | undefined }) {
    MeshBuilder.addPrimitive(state, setSphereMat(tmpSphereMat, center, radius), getSphere(detail, options));
}
