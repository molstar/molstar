/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Vec3 } from '../../linear-algebra';
import { Box3D } from '../primitives/box3d';
import { Ray3D } from '../primitives/ray3d';

describe('ray3d', () => {
    it('intersectBox3D', () => {
        const box = Box3D.create(Vec3.create(-1, -1, -1), Vec3.create(1, 1, 1));
        const out = Vec3();

        // 1. Ray starts outside and hits the box frontally
        const ray1 = Ray3D.create(Vec3.create(-2, 0, 0), Vec3.create(1, 0, 0));
        expect(Ray3D.intersectBox3D(out, ray1, box)).toBe(true);
        expect(out).toEqual(Vec3.create(-1, 0, 0));

        // 2. Ray grazes along the top edge (tangential)
        const ray2 = Ray3D.create(Vec3.create(-2, 1, 0), Vec3.create(1, 0, 0));
        expect(Ray3D.intersectBox3D(out, ray2, box)).toBe(true);
        expect(out).toEqual(Vec3.create(-1, 1, 0));

        // 3. Ray starts exactly on the surface and goes inward
        const ray3 = Ray3D.create(Vec3.create(-1, 0, 0), Vec3.create(1, 0, 0));
        expect(Ray3D.intersectBox3D(out, ray3, box)).toBe(true);
        expect(out).toEqual(Vec3.create(-1, 0, 0));

        // 4. Ray grazes a corner exactly
        const ray4 = Ray3D.create(Vec3.create(-2, -2, -2), Vec3.create(1, 1, 1));
        expect(Ray3D.intersectBox3D(out, ray4, box)).toBe(true);
        expect(out).toEqual(Vec3.create(-1, -1, -1));

        // 5. Ray starts inside the box and exits
        const ray5 = Ray3D.create(Vec3.create(0, 0, 0), Vec3.create(1, 0, 0));
        expect(Ray3D.intersectBox3D(out, ray5, box)).toBe(false);

        // 6. Ray starts outside and points away (misses box completely)
        const ray6 = Ray3D.create(Vec3.create(-2, 2, 0), Vec3.create(1, 0, 0));
        expect(Ray3D.intersectBox3D(out, ray6, box)).toBe(false);
    });
});