/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../linear-algebra';
import { Plane3D } from '../primitives/plane3d';

describe('plane3d', () => {
    it('fromNormalAndCoplanarPoint', () => {
        const normal = Vec3.create(1, 1, 1);
        Vec3.normalize(normal, normal);
        const p = Plane3D();
        Plane3D.fromNormalAndCoplanarPoint(p, normal, Vec3.zero());

        expect(p.normal).toEqual(normal);
        expect(p.constant).toBe(-0);
    });

    it('fromCoplanarPoints', () => {
        const a = Vec3.create(2.0, 0.5, 0.25);
        const b = Vec3.create(2.0, -0.5, 1.25);
        const c = Vec3.create(2.0, -3.5, 2.2);
        const p = Plane3D();
        Plane3D.fromCoplanarPoints(p, a, b, c);

        expect(p.normal).toEqual(Vec3.create(1, 0, 0));
        expect(p.constant).toBe(-2);
    });

    it('distanceToPoint', () => {
        const p = Plane3D.create(Vec3.create(2, 0, 0), -2);
        Plane3D.normalize(p, p);

        expect(Plane3D.distanceToPoint(p, Vec3.create(0, 0, 0))).toBe(-1);
        expect(Plane3D.distanceToPoint(p, Vec3.create(4, 0, 0))).toBe(3);
        expect(Plane3D.distanceToPoint(p, Plane3D.projectPoint(Vec3(), p, Vec3.zero()))).toBe(0);
    });
});
