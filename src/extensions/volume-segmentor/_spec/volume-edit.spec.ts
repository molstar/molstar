/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { BodyLabels } from '../labels';
import { recomputeStats, removeDustInPlace } from '../internal/volume-edit';
import { SwappedOrder, createTestVolume, offsetOf } from './test-volume';

describe('removeDustInPlace', () => {
    it('zeroes small components above the threshold and refreshes stats', () => {
        // 4x3x3 block (36 voxels) plus a single speck
        const volume = createTestVolume([10, 3, 3], (i, j, k) => (i < 4 || (i === 8 && j === 1 && k === 1)) ? 1 : 0, SwappedOrder);
        const data = volume.grid.cells.data as Float32Array;
        expect(volume.grid.stats.mean).toBeCloseTo(37 / 90, 6);
        const store = BodyLabels.ensure(volume);
        store.labels[offsetOf(volume, 0, 0, 0)] = 1;
        volume._propertyData['some-cache'] = { stale: true };

        expect(removeDustInPlace(volume, 5, 0.5)).toBe(1);
        expect(data[offsetOf(volume, 8, 1, 1)]).toBe(0);
        // caches are dropped, labels survive
        expect(volume._propertyData['some-cache']).toBeUndefined();
        expect(BodyLabels.get(volume)).toBe(store);
        expect(data[offsetOf(volume, 0, 0, 0)]).toBe(1);
        expect(volume.grid.stats.mean).toBeCloseTo(36 / 90, 6);
        expect(volume.grid.stats.max).toBe(1);

        // nothing left to remove
        expect(removeDustInPlace(volume, 5, 0.5)).toBe(0);
        // the block itself goes when the minimum is larger
        expect(removeDustInPlace(volume, 100, 0.5)).toBe(36);
        expect(volume.grid.stats.max).toBe(0);
    });

    it('recomputeStats matches the data', () => {
        const volume = createTestVolume([2, 2, 1], i => i === 0 ? 2 : 4);
        (volume.grid.cells.data as Float32Array)[offsetOf(volume, 1, 1, 0)] = 8;
        recomputeStats(volume);
        expect(volume.grid.stats.min).toBe(2);
        expect(volume.grid.stats.max).toBe(8);
        expect(volume.grid.stats.mean).toBeCloseTo((2 + 2 + 4 + 8) / 4, 6);
    });
});
