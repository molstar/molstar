/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Camera } from '../../../mol-canvas3d/camera';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { SyncRuntimeContext } from '../../../mol-task/execution/synchronous';
import { computeCandidates } from '../candidates';
import type { ViewMask } from '../types';
import { selectByViews } from '../view-selection';
import { CanonicalOrder, SwappedOrder, createTestVolume, offsetOf } from './test-volume';

/** Orthographic view down -z onto an 8^3 grid, canvas and viewport 100 x 100. */
function view(polygon: [number, number][], inverted = false): ViewMask {
    const camera = new Camera();
    camera.setState({
        mode: 'orthographic',
        target: Vec3.create(3.5, 3.5, 3.5),
        position: Vec3.create(3.5, 3.5, 40),
        up: Vec3.create(0, 1, 0),
        radius: 8,
    });
    return {
        id: 'v', label: 'View', polygon,
        canvasWidth: 100, canvasHeight: 100,
        viewportWidth: 100, viewportHeight: 100,
        cameraSnapshot: camera.state,
        inverted,
    };
}

describe('selectByViews', () => {
    for (const order of [CanonicalOrder, SwappedOrder]) {
        it(`selects only voxels inside the polygon (axis order ${order})`, async () => {
            // every voxel is above the threshold, so candidates cover the whole grid
            const volume = createTestVolume([8, 8, 8], () => 1, order);
            const candidates = computeCandidates(volume, 0.5);
            expect(candidates.length).toBe(8 * 8 * 8);

            // the camera looks down -z, so screen x grows with grid x and screen y falls with grid y:
            // the left half of the canvas is low x
            const out = new Uint8Array(volume.grid.cells.data.length);
            const selected = await selectByViews(volume, candidates, [view([[0, 0], [50, 0], [50, 100], [0, 100]])], out, SyncRuntimeContext);

            expect(selected).toBeGreaterThan(0);
            expect(selected).toBeLessThan(candidates.length);
            // a voxel on each side of the divide, at the same depth
            const inside = out[offsetOf(volume, 0, 4, 4)];
            const outside = out[offsetOf(volume, 7, 4, 4)];
            expect(inside).not.toBe(outside);
            // the selection is a slab through z: both ends agree with the middle
            expect(out[offsetOf(volume, 0, 4, 0)]).toBe(inside);
            expect(out[offsetOf(volume, 0, 4, 7)]).toBe(inside);
        });
    }

    it('inverting a view selects the complement', async () => {
        const volume = createTestVolume([8, 8, 8], () => 1);
        const candidates = computeCandidates(volume, 0.5);
        const polygon: [number, number][] = [[0, 0], [50, 0], [50, 100], [0, 100]];

        const plain = new Uint8Array(volume.grid.cells.data.length);
        const n = await selectByViews(volume, candidates, [view(polygon)], plain, SyncRuntimeContext);
        const flipped = new Uint8Array(volume.grid.cells.data.length);
        const m = await selectByViews(volume, candidates, [view(polygon, true)], flipped, SyncRuntimeContext);

        expect(n + m).toBe(candidates.length);
        for (let i = 0; i < candidates.length; i++) {
            const o = candidates[i];
            expect(plain[o] === flipped[o]).toBe(false);
        }
    });

    it('no views selects nothing', async () => {
        const volume = createTestVolume([4, 4, 4], () => 1);
        const out = new Uint8Array(volume.grid.cells.data.length);
        expect(await selectByViews(volume, computeCandidates(volume, 0.5), [], out, SyncRuntimeContext)).toBe(0);
        expect(out.some(v => v !== 0)).toBe(false);
    });
});
