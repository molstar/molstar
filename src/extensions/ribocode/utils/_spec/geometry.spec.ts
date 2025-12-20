/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 */
import { Vec3, Mat4 } from '../../../../mol-math/linear-algebra';
import { computeCentroid } from '../geometry';

describe('computeCentroid', () => {
    test('computes centroid of given coordinates', () => {
        const coords = [Vec3.create(1, 2, 3), Vec3.create(4, 5, 6), Vec3.create(7, 8, 9)];
        const centroid = computeCentroid(coords);
        expect(centroid).toEqual(Vec3.create(4, 5, 6));
    });
});