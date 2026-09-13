/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Turning a binary voxel selection into a soft-edged mask — shared by the volume tools.
 */

import { squaredDistanceTransform3D } from '../../mol-math/geometry/distance-transform';
import { Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { RuntimeContext } from '../../mol-task';
import type { GridBox } from './types';

/**
 * Mask value at distance `d` (in voxels) from the selection: `1` up to `extend`, then a raised
 * cosine falling to `0` over `softEdge + 1` voxels.
 */
export function softValue(d: number, extend: number, softEdge: number): number {
    if (d <= extend) return 1;
    const width = softEdge + 1;
    return 0.5 + 0.5 * Math.cos(Math.PI * Math.min((d - extend) / width, 1));
}

/** Tight bounding box of the voxels whose value equals `match`, or `undefined` if there are none. */
export function voxelBBox(values: Uint8Array, space: Tensor.Space, match: number): GridBox | undefined {
    const c = [0, 0, 0];
    let found = false;
    const min = Vec3.create(Infinity, Infinity, Infinity);
    const max = Vec3.create(-Infinity, -Infinity, -Infinity);
    for (let o = 0, n = values.length; o < n; o++) {
        if (values[o] !== match) continue;
        found = true;
        space.getCoords(o, c);
        if (c[0] < min[0]) min[0] = c[0];
        if (c[1] < min[1]) min[1] = c[1];
        if (c[2] < min[2]) min[2] = c[2];
        if (c[0] > max[0]) max[0] = c[0];
        if (c[1] > max[1]) max[1] = c[1];
        if (c[2] > max[2]) max[2] = c[2];
    }
    if (!found) return undefined;
    return { min, dims: Vec3.create(max[0] - min[0] + 1, max[1] - min[1] + 1, max[2] - min[2] + 1) };
}

/** Grows `box` by `pad` voxels on every side, clamped to the grid. */
export function padGridBox(box: GridBox, pad: number, dimensions: ArrayLike<number>): GridBox {
    const min = Vec3();
    const dims = Vec3();
    for (let a = 0; a < 3; a++) {
        const lo = Math.max(0, box.min[a] - pad);
        const hi = Math.min(dimensions[a] - 1, box.min[a] + box.dims[a] - 1 + pad);
        min[a] = lo;
        dims[a] = hi - lo + 1;
    }
    return { min, dims };
}

/**
 * Extends a binary voxel selection by `extend` voxels and adds a raised-cosine falloff, writing
 * the result on the full grid.
 *
 * The distance transform only runs inside the selection's bounding box padded by the reach of the
 * edge, since every voxel outside that is `0` anyway — for a typical map that is a few percent of
 * the grid. Distances are exact Euclidean, so the falloff is smooth rather than quantised to the
 * few steps a chamfer transform can produce.
 */
export async function softMaskFromBinary(selected: Uint8Array, space: Tensor.Space, extend: number, softEdge: number, ctx: RuntimeContext): Promise<Float32Array> {
    const out = new Float32Array(selected.length);
    const tight = voxelBBox(selected, space, 1);
    if (!tight) return out;

    const dimensions = space.dimensions;
    const box = padGridBox(tight, extend + softEdge + 1, dimensions);
    const [bx, by, bz] = box.dims;
    const [x0, y0, z0] = box.min;
    const n = bx * by * bz;

    // gather the selection into the box in x-fastest order, which the transform expects
    const inside = new Uint8Array(n);
    for (let z = 0; z < bz; z++) {
        for (let y = 0; y < by; y++) {
            const row = y * bx + z * bx * by;
            for (let x = 0; x < bx; x++) {
                if (selected[space.dataOffset(x0 + x, y0 + y, z0 + z)]) inside[row + x] = 1;
            }
        }
    }

    await ctx.update({ message: 'Computing distance transform…' });
    const distSq = squaredDistanceTransform3D(inside, bx, by, bz);

    await ctx.update({ message: 'Applying soft edge…' });
    for (let z = 0; z < bz; z++) {
        for (let y = 0; y < by; y++) {
            const row = y * bx + z * bx * by;
            for (let x = 0; x < bx; x++) {
                const v = softValue(Math.sqrt(distSq[row + x]), extend, softEdge);
                if (v > 0) out[space.dataOffset(x0 + x, y0 + y, z0 + z)] = v;
            }
        }
    }
    return out;
}
