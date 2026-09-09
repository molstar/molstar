/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Camera } from '../../../../mol-canvas3d/camera';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { SyncRuntimeContext } from '../../../../mol-task/execution/synchronous';
import { computeCandidates } from '../internal/candidates';
import { assignPolygons, assignRemainder, countUnassigned, countVoxels } from '../internal/label-ops';
import { BodyLabels } from '../labels';
import { ViewMask } from '../types';
import { CanonicalOrder, SwappedOrder, createTestVolume, offsetOf } from './test-volume';

/** Orthographic view down -z onto an 8^3 grid, canvas and viewport 100 x 100. */
function viewMask(polygon: [number, number][], inverted = false): ViewMask {
    const snapshot = Camera.createDefaultSnapshot();
    snapshot.mode = 'orthographic';
    snapshot.position = Vec3.create(3.5, 3.5, 100);
    snapshot.target = Vec3.create(3.5, 3.5, 0);
    snapshot.up = Vec3.create(0, 1, 0);
    snapshot.radius = 20;
    snapshot.radiusMax = 200;
    return {
        id: 'm', label: 'm', polygon, canvasWidth: 100, canvasHeight: 100,
        viewportWidth: 100, viewportHeight: 100, cameraSnapshot: snapshot, inverted,
    };
}

const LeftHalf: [number, number][] = [[0, 0], [50, 0], [50, 100], [0, 100]];
const TopHalf: [number, number][] = [[0, 0], [100, 0], [100, 50], [0, 50]];

function labelAt(labels: Uint8Array, volume: ReturnType<typeof createTestVolume>, i: number, j: number, k: number) {
    return labels[offsetOf(volume, i, j, k)];
}

describe('assignPolygons', () => {
    for (const order of [CanonicalOrder, SwappedOrder]) {
        it(`labels voxels projecting inside the polygon (axis order ${order})`, async () => {
            const volume = createTestVolume([8, 8, 8], () => 1, order);
            const store = BodyLabels.ensure(volume);
            const candidates = computeCandidates(volume, 0.5);
            expect(candidates.length).toBe(512);

            const changed = await assignPolygons(store.labels, candidates, volume, [viewMask(LeftHalf)], 1, 'replace', SyncRuntimeContext);
            expect(changed).toBe(256);
            for (let i = 0; i < 8; i++) for (let j = 0; j < 8; j++) for (let k = 0; k < 8; k++) {
                expect(labelAt(store.labels, volume, i, j, k)).toBe(i < 4 ? 1 : 0);
            }
        });
    }

    it('intersects several views and honours inverted polygons', async () => {
        const volume = createTestVolume([8, 8, 8], () => 1);
        const store = BodyLabels.ensure(volume);
        const candidates = computeCandidates(volume, 0.5);

        await assignPolygons(store.labels, candidates, volume, [viewMask(LeftHalf), viewMask(TopHalf)], 2, 'replace', SyncRuntimeContext);
        // canvas top = high world y
        for (let i = 0; i < 8; i++) for (let j = 0; j < 8; j++) {
            expect(labelAt(store.labels, volume, i, j, 3)).toBe(i < 4 && j >= 4 ? 2 : 0);
        }

        await assignPolygons(store.labels, candidates, volume, [viewMask(LeftHalf, true)], 3, 'replace', SyncRuntimeContext);
        for (let i = 0; i < 8; i++) expect(labelAt(store.labels, volume, i, 0, 0)).toBe(i < 4 ? 0 : 3);
    });

    it('supports unassigned-only and erase modes and skips voxels below threshold', async () => {
        const volume = createTestVolume([8, 8, 8], (i, j, k) => k === 0 ? 0 : 1);
        const store = BodyLabels.ensure(volume);
        const candidates = computeCandidates(volume, 0.5);
        expect(candidates.length).toBe(448);

        await assignPolygons(store.labels, candidates, volume, [viewMask(TopHalf)], 1, 'replace', SyncRuntimeContext);
        const changed = await assignPolygons(store.labels, candidates, volume, [viewMask(LeftHalf)], 2, 'unassigned-only', SyncRuntimeContext);
        expect(changed).toBe(4 * 4 * 7);
        expect(labelAt(store.labels, volume, 0, 7, 1)).toBe(1);
        expect(labelAt(store.labels, volume, 0, 0, 1)).toBe(2);
        expect(labelAt(store.labels, volume, 0, 0, 0)).toBe(0);

        await assignPolygons(store.labels, candidates, volume, [viewMask(LeftHalf)], 1, 'erase', SyncRuntimeContext);
        expect(labelAt(store.labels, volume, 0, 7, 1)).toBe(0);
        expect(labelAt(store.labels, volume, 7, 7, 1)).toBe(1);
        expect(labelAt(store.labels, volume, 0, 0, 1)).toBe(2);
    });

    it('changes nothing without polygons', async () => {
        const volume = createTestVolume([2, 2, 2], () => 1);
        const store = BodyLabels.ensure(volume);
        expect(await assignPolygons(store.labels, computeCandidates(volume, 0), volume, [], 1, 'replace', SyncRuntimeContext)).toBe(0);
    });
});

describe('remainder and counts', () => {
    it('assigns the remainder and counts voxels', () => {
        const volume = createTestVolume([4, 4, 4], i => i);
        const store = BodyLabels.ensure(volume);
        const candidates = computeCandidates(volume, 2); // i = 2, 3 -> 32 voxels
        store.labels[offsetOf(volume, 2, 0, 0)] = 1;
        expect(countUnassigned(store.labels, candidates)).toBe(31);

        expect(assignRemainder(store.labels, candidates, 5)).toBe(31);
        const counts = countVoxels(store.labels);
        expect(counts[5]).toBe(31);
        expect(counts[1]).toBe(1);
        expect(counts[0]).toBe(32);
        expect(countUnassigned(store.labels, candidates)).toBe(0);
        expect(assignRemainder(store.labels, candidates, 5)).toBe(0);
    });
});
