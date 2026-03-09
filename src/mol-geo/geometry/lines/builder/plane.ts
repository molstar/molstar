/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LinesBuilder } from '../lines-builder';
import { Vec3 } from '../../../../mol-math/linear-algebra/3d/vec3';
import { Mat4 } from '../../../../mol-math/linear-algebra/3d/mat4';

const _c0 = Vec3();
const _c1 = Vec3();
const _c2 = Vec3();
const _c3 = Vec3();

/**
 * Add wireframe edges of a quad to a LinesBuilder.
 */
export function addPlane(builder: LinesBuilder, corners: ArrayLike<number>, transform: Mat4, group: number) {
    Vec3.fromArray(_c0, corners, 0);
    Vec3.fromArray(_c1, corners, 3);
    Vec3.fromArray(_c2, corners, 6);
    Vec3.fromArray(_c3, corners, 9);

    Vec3.transformMat4(_c0, _c0, transform);
    Vec3.transformMat4(_c1, _c1, transform);
    Vec3.transformMat4(_c2, _c2, transform);
    Vec3.transformMat4(_c3, _c3, transform);

    builder.addVec(_c0, _c1, group);
    builder.addVec(_c1, _c2, group);
    builder.addVec(_c2, _c3, group);
    builder.addVec(_c3, _c0, group);
}
