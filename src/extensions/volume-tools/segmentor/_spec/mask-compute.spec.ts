/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Grid } from '../../../../mol-model/volume';
import { Mat4, Vec3 } from '../../../../mol-math/linear-algebra';
import { SyncRuntimeContext } from '../../../../mol-task/execution/synchronous';
import { buildCroppedMaskVolume, computeBodyMask, scatterToFullBox } from '../internal/mask-compute';
import { BodyLabels } from '../labels';
import { CanonicalOrder, SwappedOrder, createTestVolume, offsetOf } from '../../_spec/test-volume';

describe('computeBodyMask', () => {
    for (const order of [CanonicalOrder, SwappedOrder]) {
        it(`dilates a single voxel into a sphere (axis order ${order})`, async () => {
            const volume = createTestVolume([11, 11, 11], () => 1, order);
            const store = BodyLabels.ensure(volume);
            store.labels[offsetOf(volume, 5, 5, 5)] = 3;

            const result = (await computeBodyMask(volume, store.labels, 3, { extend: 2, softEdge: 0, pruneBelowThreshold: true }, 0.5, SyncRuntimeContext))!;
            expect(result.voxelCount).toBe(1);
            expect(Array.from(result.box.min)).toEqual([2, 2, 2]);
            expect(Array.from(result.box.dims)).toEqual([7, 7, 7]);

            const full = scatterToFullBox(result, volume.grid.cells.space, new Float32Array(11 ** 3));
            const at = (i: number, j: number, k: number) => full[offsetOf(volume, i, j, k)];
            expect(at(5, 5, 5)).toBe(1);
            expect(at(5, 5, 7)).toBe(1); // d = 2
            expect(at(5, 5, 8)).toBeCloseTo(0, 6); // d = 3 = extend + softEdge + 1
            expect(at(5, 6, 7)).toBeGreaterThan(0); // d = sqrt(5)
            expect(at(5, 6, 7)).toBeLessThan(1);
            expect(at(0, 0, 0)).toBe(0);
        });
    }

    it('returns a binary mask when extend and softEdge are 0, honouring the threshold', async () => {
        const volume = createTestVolume([4, 4, 4], i => i);
        const store = BodyLabels.ensure(volume);
        for (let i = 0; i < 4; i++) store.labels[offsetOf(volume, i, 1, 1)] = 1;

        const pruned = (await computeBodyMask(volume, store.labels, 1, { extend: 0, softEdge: 0, pruneBelowThreshold: true }, 2, SyncRuntimeContext))!;
        expect(pruned.voxelCount).toBe(2);
        const kept = (await computeBodyMask(volume, store.labels, 1, { extend: 0, softEdge: 0, pruneBelowThreshold: false }, 2, SyncRuntimeContext))!;
        expect(kept.voxelCount).toBe(4);
        expect(Array.from(kept.data).every(v => v === 0 || v === 1)).toBe(true);
    });

    it('returns undefined for a body without voxels', async () => {
        const volume = createTestVolume([3, 3, 3], () => 1);
        const store = BodyLabels.ensure(volume);
        expect(await computeBodyMask(volume, store.labels, 1, { extend: 1, softEdge: 1, pruneBelowThreshold: true }, 0, SyncRuntimeContext)).toBeUndefined();
    });
});

describe('buildCroppedMaskVolume', () => {
    it('places the cropped grid at the box origin in world space', async () => {
        const matrix = Mat4.fromScaling(Mat4(), Vec3.create(2, 2, 2));
        Mat4.setTranslation(matrix, Vec3.create(10, 20, 30));
        const volume = createTestVolume([9, 9, 9], () => 1, SwappedOrder, matrix);
        const store = BodyLabels.ensure(volume);
        store.labels[offsetOf(volume, 6, 4, 5)] = 1;

        const result = (await computeBodyMask(volume, store.labels, 1, { extend: 1, softEdge: 0, pruneBelowThreshold: true }, 0, SyncRuntimeContext))!;
        const cropped = buildCroppedMaskVolume(volume, result, 'Body 1');

        const sourceG2C = Grid.getGridToCartesianTransform(volume.grid);
        const croppedG2C = Grid.getGridToCartesianTransform(cropped.grid);
        const expected = Vec3.transformMat4(Vec3(), result.box.min, sourceG2C);
        const actual = Vec3.transformMat4(Vec3(), Vec3.create(0, 0, 0), croppedG2C);
        expect(Array.from(actual)).toEqual(Array.from(expected));

        const local = Vec3.create(6 - result.box.min[0], 4 - result.box.min[1], 5 - result.box.min[2]);
        const center = Vec3.transformMat4(Vec3(), local, croppedG2C);
        expect(Array.from(center)).toEqual([10 + 12, 20 + 8, 30 + 10]);
        expect(cropped.grid.stats.max).toBe(1);
        expect(cropped.grid.stats.min).toBe(0);
    });
});
