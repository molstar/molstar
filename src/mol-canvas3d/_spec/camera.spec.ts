/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { Viewport, cameraProject, cameraUnproject } from '../camera/util';

describe('camera', () => {
    it('project/unproject', () => {
        const proj = Mat4.perspective(Mat4(), -1, 1, 1, -1, 1, 100);
        const invProj = Mat4.invert(Mat4(), proj);

        const c = Vec4();
        const po = Vec3();

        const vp = Viewport.create(0, 0, 100, 100);
        const pi = Vec3.create(0, 0, 1);
        cameraProject(c, pi, vp, proj);
        expect(Vec4.equals(c, Vec4.create(50, 50, 2.020202, -1))).toBe(true);
        cameraUnproject(po, c, vp, invProj);
        expect(Vec3.equals(po, pi)).toBe(true);

        Vec3.set(pi, 0.5, 0.5, 1);
        cameraProject(c, pi, vp, proj);
        cameraUnproject(po, c, vp, invProj);
        expect(Vec3.equals(po, pi)).toBe(true);

        Viewport.set(vp, 50, 50, 100, 100);
        Vec3.set(pi, 0.5, 0.5, 1);
        cameraProject(c, pi, vp, proj);
        cameraUnproject(po, c, vp, invProj);
        expect(Vec3.equals(po, pi)).toBe(true);
    });
});