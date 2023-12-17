/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '../../linear-algebra';
import { Box3D } from '../primitives/box3d';
import { Frustum3D } from '../primitives/frustum3d';
import { Sphere3D } from '../primitives/sphere3d';

const v3 = Vec3.create;
const s3 = Sphere3D.create;

describe('frustum3d', () => {
    it('intersectsSphere3D', () => {
        const f = Frustum3D();
        const m = Mat4.perspective(Mat4(), -1, 1, 1, -1, 1, 100);
        Frustum3D.fromProjectionMatrix(f, m);

        expect(Frustum3D.intersectsSphere3D(f, s3(v3(0, 0, 0), 0))).toBe(false);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(0, 0, 0), 0.9))).toBe(false);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(0, 0, 0), 1.1))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(0, 0, -50), 0))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(0, 0, -1.001), 0))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(-1, -1, -1.001), 0))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(-1.1, -1.1, -1.001), 0))).toBe(false);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(-1.1, -1.1, -1.001), 0.5))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(1, 1, -1.001), 0))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(1.1, 1.1, -1.001), 0))).toBe(false);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(1.1, 1.1, -1.001), 0.5))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(0, 0, -99.999), 0))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(-99.999, -99.999, -99.999), 0))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(-100.1, -100.1, -100.1), 0))).toBe(false);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(-100.1, -100.1, -100.1), 0.5))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(99.999, 99.999, -99.999), 0))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(100.1, 100.1, -100.1), 0))).toBe(false);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(100.1, 100.1, -100.1), 0.2))).toBe(true);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(0, 0, -101), 0))).toBe(false);
        expect(Frustum3D.intersectsSphere3D(f, s3(v3(0, 0, -101), 1.1))).toBe(true);
    });

    it('intersectsBox3D', () => {
        const f = Frustum3D();
        const m = Mat4.perspective(Mat4(), -1, 1, 1, -1, 1, 100);
        Frustum3D.fromProjectionMatrix(f, m);

        const b0 = Box3D.create(v3(0, 0, 0), v3(1, 1, 1));
        expect(Frustum3D.intersectsBox3D(f, b0)). toBe(false);

        const b1 = Box3D.create(v3(-1.1, -1.1, -1.1), v3(-0.1, -0.1, -0.1));
        expect(Frustum3D.intersectsBox3D(f, b1)). toBe(true);
    });

    it('containsPoint', () => {
        const f = Frustum3D();
        const m = Mat4.perspective(Mat4(), -1, 1, 1, -1, 1, 100);
        Frustum3D.fromProjectionMatrix(f, m);

        expect(Frustum3D.containsPoint(f, v3(0, 0, 0))).toBe(false);
        expect(Frustum3D.containsPoint(f, v3(0, 0, -50))).toBe(true);
        expect(Frustum3D.containsPoint(f, v3(0, 0, -1.001))).toBe(true);
        expect(Frustum3D.containsPoint(f, v3(-1, -1, -1.001))).toBe(true);
        expect(Frustum3D.containsPoint(f, v3(-1.1, -1.1, -1.001))).toBe(false);
        expect(Frustum3D.containsPoint(f, v3(1, 1, -1.001))).toBe(true);
        expect(Frustum3D.containsPoint(f, v3(1.1, 1.1, -1.001))).toBe(false);
        expect(Frustum3D.containsPoint(f, v3(0, 0, -99.999))).toBe(true);
        expect(Frustum3D.containsPoint(f, v3(-99.999, -99.999, -99.999))).toBe(true);
        expect(Frustum3D.containsPoint(f, v3(-100.1, -100.1, -100.1))).toBe(false);
        expect(Frustum3D.containsPoint(f, v3(99.999, 99.999, -99.999))).toBe(true);
        expect(Frustum3D.containsPoint(f, v3(100.1, 100.1, -100.1))).toBe(false);
        expect(Frustum3D.containsPoint(f, v3(0, 0, -101))).toBe(false);
    });
});
