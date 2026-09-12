/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Selecting voxels with the polygons drawn over the viewport — shared by the volume tools.
 */

import { pointInPolygon2D } from '../../mol-math/geometry/polygon';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Grid, Volume } from '../../mol-model/volume';
import { RuntimeContext } from '../../mol-task';
import { PreparedMask, prepareMask, projectToNormInPlace } from './view-projection';
import type { ViewMask } from './types';

const UpdateInterval = 1 << 18;

const tmpCoords = [0, 0, 0];
const tmpPos = Vec3();
const tmpG2C = Mat4();
const tmpNorm: [number, number] = [0, 0];

/**
 * True when `pos` projects inside every view's polygon, honouring each view's `inverted` flag.
 * A voxel belongs to the selection only if it survives all of them.
 */
export function passesAllViews(pos: Vec3, prepared: readonly PreparedMask[], out: [number, number] = tmpNorm): boolean {
    for (let m = 0; m < prepared.length; m++) {
        projectToNormInPlace(pos, prepared[m], out);
        const inside = pointInPolygon2D(out[0], out[1], prepared[m].normPolygon);
        if (prepared[m].inverted ? inside : !inside) return false;
    }
    return true;
}

/**
 * Marks every candidate voxel selected by `masks` with `1` in `out` (untouched voxels keep
 * their current value, so clear it first when recomputing). With `invert`, the candidates the
 * views reject are marked instead. Returns the number of voxels selected.
 *
 * Integer grid coordinates are projected, not cell centres, because the rendered isosurface
 * places grid index `i` at coordinate `i` — a half-voxel offset here shows up as speckle when
 * the result is read back to color that surface.
 */
export async function selectByViews(volume: Volume, candidates: Int32Array, masks: readonly ViewMask[], out: Uint8Array, ctx: RuntimeContext, invert = false): Promise<number> {
    if (masks.length === 0) return 0;

    const { space } = volume.grid.cells;
    Mat4.copy(tmpG2C, Grid.getGridToCartesianTransform(volume.grid));
    const prepared = masks.map(prepareMask);
    let selected = 0;

    for (let i = 0, n = candidates.length; i < n; i++) {
        if (i % UpdateInterval === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Selecting voxels…', current: i, max: n });
        }
        const offset = candidates[i];
        space.getCoords(offset, tmpCoords);
        Vec3.set(tmpPos, tmpCoords[0], tmpCoords[1], tmpCoords[2]);
        Vec3.transformMat4(tmpPos, tmpPos, tmpG2C);

        if (passesAllViews(tmpPos, prepared) !== invert) {
            out[offset] = 1;
            selected++;
        }
    }
    return selected;
}
