/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../3d/vec3';
import { leastObstructedDirection } from '../3d/optimize-direction';

describe('OptimizeDirection', () => {
    it('works more or less as expected', () => {
        const points: Vec3[] = [
            Vec3.create(1, 0, 0),
            Vec3.create(-1, 0, 0),
            Vec3.create(0, 1, 0),
            Vec3.create(0, -1, 0),
            Vec3.create(0, 0, 1),
        ];
        const dir = leastObstructedDirection(points);

        console.log('dir', dir);
        expect(dir).toBeDefined();
        expect(dir[0]).toBeCloseTo(0, 6);
        expect(dir[1]).toBeCloseTo(0, 6);
        expect(dir[2]).toBeCloseTo(-1, 6);
    });
});