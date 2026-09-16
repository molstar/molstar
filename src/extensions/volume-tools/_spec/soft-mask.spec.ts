/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { SyncRuntimeContext } from '../../../mol-task/execution/synchronous';
import { padGridBox, softMaskFromBinary, softValue, voxelBBox } from '../soft-mask';
import { CanonicalOrder, SwappedOrder, createTestVolume, offsetOf } from './test-volume';

describe('softValue', () => {
    it('is 1 within extend, 0 beyond extend + softEdge + 1 and monotone in between', () => {
        expect(softValue(0, 2, 3)).toBe(1);
        expect(softValue(2, 2, 3)).toBe(1);
        expect(softValue(6, 2, 3)).toBeCloseTo(0, 6);
        expect(softValue(9, 2, 3)).toBeCloseTo(0, 6);
        expect(softValue(4, 2, 3)).toBeCloseTo(0.5, 6);
        let prev = 1;
        for (let d = 0; d <= 8; d += 0.25) {
            const v = softValue(d, 2, 3);
            expect(v).toBeLessThanOrEqual(prev + 1e-9);
            prev = v;
        }
    });
});

describe('voxelBBox / padGridBox', () => {
    it('finds the tight box and clamps padding to the grid', () => {
        const volume = createTestVolume([6, 5, 4], () => 1, SwappedOrder);
        const space = volume.grid.cells.space;
        const values = new Uint8Array(volume.grid.cells.data.length);
        values[offsetOf(volume, 1, 2, 0)] = 1;
        values[offsetOf(volume, 4, 2, 3)] = 1;

        const box = voxelBBox(values, space, 1)!;
        expect(Array.from(box.min)).toEqual([1, 2, 0]);
        expect(Array.from(box.dims)).toEqual([4, 1, 4]);

        const padded = padGridBox(box, 2, space.dimensions);
        expect(Array.from(padded.min)).toEqual([0, 0, 0]);
        expect(Array.from(padded.dims)).toEqual([6, 5, 4]);

        expect(voxelBBox(values, space, 2)).toBeUndefined();
    });
});

describe('softMaskFromBinary', () => {
    for (const order of [CanonicalOrder, SwappedOrder]) {
        it(`matches an exact Euclidean falloff (axis order ${order})`, async () => {
            const dims: [number, number, number] = [13, 13, 13];
            const volume = createTestVolume(dims, () => 0, order);
            const space = volume.grid.cells.space;
            const selected = new Uint8Array(volume.grid.cells.data.length);
            selected[offsetOf(volume, 6, 6, 6)] = 1;

            const extend = 2, softEdge = 3;
            const out = await softMaskFromBinary(selected, space, extend, softEdge, SyncRuntimeContext);

            // brute force: distance from the single selected voxel
            for (let i = 0; i < dims[0]; i++) {
                for (let j = 0; j < dims[1]; j++) {
                    for (let k = 0; k < dims[2]; k++) {
                        const d = Math.sqrt((i - 6) ** 2 + (j - 6) ** 2 + (k - 6) ** 2);
                        expect(out[offsetOf(volume, i, j, k)]).toBeCloseTo(softValue(d, extend, softEdge), 5);
                    }
                }
            }
        });
    }

    it('is 1 across the core and never leaves [0, 1]', async () => {
        const volume = createTestVolume([16, 16, 16], () => 0);
        const space = volume.grid.cells.space;
        const selected = new Uint8Array(volume.grid.cells.data.length);
        for (let i = 4; i <= 11; i++) for (let j = 4; j <= 11; j++) for (let k = 4; k <= 11; k++) {
            selected[offsetOf(volume, i, j, k)] = 1;
        }

        const out = await softMaskFromBinary(selected, space, 1, 2, SyncRuntimeContext);
        expect(out[offsetOf(volume, 7, 7, 7)]).toBe(1);
        // one voxel out is still inside `extend`
        expect(out[offsetOf(volume, 3, 7, 7)]).toBe(1);
        for (let o = 0; o < out.length; o++) {
            expect(out[o]).toBeGreaterThanOrEqual(0);
            expect(out[o]).toBeLessThanOrEqual(1);
        }
    });

    it('an empty selection gives an all-zero mask', async () => {
        const volume = createTestVolume([8, 8, 8], () => 0);
        const out = await softMaskFromBinary(new Uint8Array(8 ** 3), volume.grid.cells.space, 2, 2, SyncRuntimeContext);
        expect(out.some(v => v !== 0)).toBe(false);
    });
});
