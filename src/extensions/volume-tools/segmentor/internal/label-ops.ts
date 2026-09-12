/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Grid, Volume } from '../../../../mol-model/volume';
import { Mat4, Vec3 } from '../../../../mol-math/linear-algebra';
import { RuntimeContext } from '../../../../mol-task';
import { passesAllViews } from '../../view-selection';
import { prepareMask } from '../../view-projection';
import { AssignMode, BodyId, MaxBodyId, ViewMask } from '../types';

const UpdateInterval = 1 << 18;

function assign(labels: Uint8Array, offset: number, bodyId: BodyId, mode: AssignMode): boolean {
    const old = labels[offset];
    switch (mode) {
        case 'replace':
            if (old === bodyId) return false;
            labels[offset] = bodyId;
            return true;
        case 'unassigned-only':
            if (old !== 0) return false;
            labels[offset] = bodyId;
            return true;
        case 'erase':
            if (old !== bodyId) return false;
            labels[offset] = 0;
            return true;
    }
}

const tmpCoords = [0, 0, 0];
const tmpPos = Vec3();
const tmpG2C = Mat4();
const tmpNorm: [number, number] = [0, 0];

/**
 * Labels every candidate voxel whose (integer) grid position projects inside all view
 * polygons (per-polygon `inverted` honoured). Integer coordinates are used because the
 * rendered isosurface places grid index `i` at coordinate `i`. Returns the number of
 * voxels changed.
 */
export async function assignPolygons(labels: Uint8Array, candidates: Int32Array, volume: Volume, masks: readonly ViewMask[], bodyId: BodyId, mode: AssignMode, ctx: RuntimeContext): Promise<number> {
    if (masks.length === 0) return 0;

    const { space } = volume.grid.cells;
    Mat4.copy(tmpG2C, Grid.getGridToCartesianTransform(volume.grid));
    const prepared = masks.map(prepareMask);
    let changed = 0;

    for (let i = 0, n = candidates.length; i < n; i++) {
        if (i % UpdateInterval === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Assigning voxels…', current: i, max: n });
        }
        const offset = candidates[i];
        space.getCoords(offset, tmpCoords);
        Vec3.set(tmpPos, tmpCoords[0], tmpCoords[1], tmpCoords[2]);
        Vec3.transformMat4(tmpPos, tmpPos, tmpG2C);

        if (passesAllViews(tmpPos, prepared, tmpNorm) && assign(labels, offset, bodyId, mode)) changed++;
    }
    return changed;
}

/** Assigns every unassigned candidate voxel to `bodyId`; returns the number of voxels changed. */
export function assignRemainder(labels: Uint8Array, candidates: Int32Array, bodyId: BodyId): number {
    let changed = 0;
    for (let i = 0, n = candidates.length; i < n; i++) {
        const offset = candidates[i];
        if (labels[offset] === 0) {
            labels[offset] = bodyId;
            changed++;
        }
    }
    return changed;
}

/** Number of voxels per label; index 0 counts unassigned voxels. */
export function countVoxels(labels: Uint8Array): Uint32Array {
    const counts = new Uint32Array(MaxBodyId + 1);
    for (let o = 0, n = labels.length; o < n; o++) counts[labels[o]]++;
    return counts;
}

/** Number of candidate voxels not assigned to any body. */
export function countUnassigned(labels: Uint8Array, candidates: Int32Array): number {
    let count = 0;
    for (let i = 0, n = candidates.length; i < n; i++) {
        if (labels[candidates[i]] === 0) count++;
    }
    return count;
}
