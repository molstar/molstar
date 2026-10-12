/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Andy Turner <agdturner@gmail.com>
 */
import { Vec3, Mat4 } from '../../../../mol-math/linear-algebra';
import { alignDatasetUsingChains, computeCentroid } from '../geometry';

describe('computeCentroid', () => {
    test('computes centroid of given coordinates', () => {
        const coords = [Vec3.create(1, 2, 3), Vec3.create(4, 5, 6), Vec3.create(7, 8, 9)];
        const centroid = computeCentroid(coords);
        expect(centroid).toEqual(Vec3.create(4, 5, 6));
    });
});

describe('alignDatasetUsingChains', () => {
    test('preserves pure translation when coordinates differ by a fixed offset', () => {
        const movingTypes = ['P', 'P', 'P', 'P'];
        const movingChainIds = ['A', 'A', 'A', 'A'];
        const movingX = [10, 11, 10, 10];
        const movingY = [20, 20, 21, 20];
        const movingZ = [30, 30, 30, 31];

        const referenceTypes = ['P', 'P', 'P', 'P'];
        const referenceChainIds = ['B', 'B', 'B', 'B'];
        const referenceX = [0, 1, 0, 0];
        const referenceY = [0, 0, 1, 0];
        const referenceZ = [0, 0, 0, 1];

        const result = alignDatasetUsingChains(
            { P: true },
            'A',
            movingTypes,
            movingChainIds,
            movingX,
            movingY,
            movingZ,
            'B',
            referenceTypes,
            referenceChainIds,
            referenceX,
            referenceY,
            referenceZ
        );

        for (let i = 0; i < referenceX.length; i++) {
            expect(result.alignedX[i]).toBeCloseTo(referenceX[i], 6);
            expect(result.alignedY[i]).toBeCloseTo(referenceY[i], 6);
            expect(result.alignedZ[i]).toBeCloseTo(referenceZ[i], 6);
        }
    });
});