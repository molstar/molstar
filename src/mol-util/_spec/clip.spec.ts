/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Plane3D } from '../../mol-math/geometry/primitives/plane3d';
import { Sphere3D } from '../../mol-math/geometry/primitives/sphere3d';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { Clip } from '../clip';

function getObjects(type: keyof typeof Clip.Type, invert: boolean, position: Vec3, scale: Vec3, transform: Mat4) {
    return Clip.getClip({
        variant: 'pixel',
        objects: [{ type, invert, position, rotation: { axis: Vec3.create(1, 0, 0), angle: 0 }, scale, transform }]
    } as Clip.Props).objects;
}

describe('clip', () => {
    it('getPlane', () => {
        const plane = Plane3D();
        const position = Vec3.create(0, 2, 0);

        Clip.getPlane(plane, getObjects('plane', false, position, Vec3.create(1, 1, 1), Mat4.identity()), 0);
        expect(Plane3D.distanceToPoint(plane, Vec3.create(5, 3, 5))).toBeCloseTo(1);

        Clip.getPlane(plane, getObjects('plane', true, position, Vec3.create(1, 1, 1), Mat4.identity()), 0);
        expect(Plane3D.distanceToPoint(plane, Vec3.create(5, 3, 5))).toBeCloseTo(-1);

        const transform = Mat4.fromTranslation(Mat4(), Vec3.create(0, 1, 0));
        Clip.getPlane(plane, getObjects('plane', false, position, Vec3.create(1, 1, 1), transform), 0);
        expect(Plane3D.distanceToPoint(plane, Vec3.create(0, 1, 0))).toBeCloseTo(0);
        expect(Plane3D.distanceToPoint(plane, Vec3.create(0, 4, 0))).toBeCloseTo(3);
    });

    it('canIntersectSphere', () => {
        const sphere = Sphere3D.create(Vec3.create(0, 0, 0), 1);
        const identity = Mat4.identity();
        const one = Vec3.create(1, 1, 1);

        expect(Clip.canIntersectSphere(getObjects('plane', false, Vec3.create(0, 0.5, 0), one, identity), 0, sphere)).toBe(true);
        expect(Clip.canIntersectSphere(getObjects('plane', false, Vec3.create(0, 2, 0), one, identity), 0, sphere)).toBe(false);

        expect(Clip.canIntersectSphere(getObjects('sphere', false, Vec3.create(0, 0, 0), Vec3.create(10, 10, 10), identity), 0, sphere)).toBe(false);
        expect(Clip.canIntersectSphere(getObjects('sphere', false, Vec3.create(5, 0, 0), Vec3.create(10, 10, 10), identity), 0, sphere)).toBe(true);
        expect(Clip.canIntersectSphere(getObjects('sphere', false, Vec3.create(20, 0, 0), Vec3.create(10, 10, 10), identity), 0, sphere)).toBe(false);

        expect(Clip.canIntersectSphere(getObjects('cube', false, Vec3.create(3, 0, 0), Vec3.create(4, 4, 4), identity), 0, sphere)).toBe(true);
        expect(Clip.canIntersectSphere(getObjects('cylinder', false, Vec3.create(20, 0, 0), Vec3.create(4, 4, 4), identity), 0, sphere)).toBe(false);
        expect(Clip.canIntersectSphere(getObjects('infiniteCone', false, Vec3.create(20, 0, 0), one, identity), 0, sphere)).toBe(true);
    });
});
