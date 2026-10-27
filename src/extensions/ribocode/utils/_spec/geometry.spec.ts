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

    test('uses translation-only alignment for single-atom fallback selections', () => {
        const movingTypes = ['N', 'C', 'P'];
        const movingChainIds = ['A', 'A', 'A'];
        const movingX = [10, 20, 30];
        const movingY = [5, 6, 7];
        const movingZ = [1, 2, 3];

        const referenceTypes = ['N', 'C', 'P'];
        const referenceChainIds = ['B', 'B', 'B'];
        const referenceX = [100, 200, 300];
        const referenceY = [50, 60, 70];
        const referenceZ = [11, 12, 13];

        const result = alignDatasetUsingChains(
            { O: true },
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

        for (let i = 0; i < movingX.length; i++) {
            expect(Number.isFinite(result.alignedX[i])).toBe(true);
            expect(Number.isFinite(result.alignedY[i])).toBe(true);
            expect(Number.isFinite(result.alignedZ[i])).toBe(true);
        }

        const dx = result.centroidReference[0] - result.centroid[0];
        const dy = result.centroidReference[1] - result.centroid[1];
        const dz = result.centroidReference[2] - result.centroid[2];

        for (let i = 0; i < movingX.length; i++) {
            expect(result.alignedX[i]).toBeCloseTo(movingX[i] + dx, 6);
            expect(result.alignedY[i]).toBeCloseTo(movingY[i] + dy, 6);
            expect(result.alignedZ[i]).toBeCloseTo(movingZ[i] + dz, 6);
        }

        expect(result.rotmat).toEqual([1, 0, 0, 0, 1, 0, 0, 0, 1]);
    });
});