/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Grid, Volume } from '../../../../mol-model/volume';
import { Mat4, Tensor, Vec3 } from '../../../../mol-math/linear-algebra';
import { CustomProperties } from '../../../../mol-model/custom-property';
import { RuntimeContext } from '../../../../mol-task';
import { BodyId, BodyMaskParams, BodyMaskResult, GridBox } from '../types';
import { squaredDistanceTransform3D } from '../../../../mol-math/geometry/distance-transform';

/**
 * Soft mask value at Euclidean distance `d` (voxels) from the binary body:
 * 1 within `extend`, then a raised cosine falling to 0 over `softEdge + 1` voxels.
 */
export function softValue(d: number, extend: number, softEdge: number): number {
    if (d <= extend) return 1;
    const width = softEdge + 1;
    return 0.5 + 0.5 * Math.cos(Math.PI * Math.min((d - extend) / width, 1));
}

/** Grid-space bounding box of all voxels labelled `bodyId`, or undefined if there are none. */
export function labelBBox(labels: Uint8Array, space: Tensor.Space, bodyId: BodyId): GridBox | undefined {
    const c = [0, 0, 0];
    let found = false;
    const min = Vec3.create(Infinity, Infinity, Infinity);
    const max = Vec3.create(-Infinity, -Infinity, -Infinity);
    for (let o = 0, n = labels.length; o < n; o++) {
        if (labels[o] !== bodyId) continue;
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

/** Grows `box` by `pad` voxels on every side, clamped to the grid dimensions. */
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
 * Soft mask of one body: binary body voxels (optionally restricted to density >= threshold),
 * dilated by `extend` voxels and given a cosine soft edge, computed only inside the padded
 * bounding box of the body. Returns undefined when the body has no voxels.
 */
export async function computeBodyMask(volume: Volume, labels: Uint8Array, bodyId: BodyId, params: BodyMaskParams, thresholdAbs: number, ctx: RuntimeContext): Promise<BodyMaskResult | undefined> {
    const { space, data } = volume.grid.cells;
    const values = data as unknown as ArrayLike<number>;
    const tight = labelBBox(labels, space, bodyId);
    if (!tight) return undefined;

    const pad = params.extend + params.softEdge + 1;
    const box = padGridBox(tight, pad, space.dimensions);
    const [bx, by, bz] = box.dims;
    const [x0, y0, z0] = box.min;
    const n = bx * by * bz;

    await ctx.update({ message: 'Collecting body voxels…' });
    const inside = new Uint8Array(n);
    let voxelCount = 0;
    for (let z = 0; z < bz; z++) {
        for (let y = 0; y < by; y++) {
            const row = y * bx + z * bx * by;
            for (let x = 0; x < bx; x++) {
                const o = space.dataOffset(x0 + x, y0 + y, z0 + z);
                if (labels[o] !== bodyId) continue;
                if (params.pruneBelowThreshold && values[o] < thresholdAbs) continue;
                inside[row + x] = 1;
                voxelCount++;
            }
        }
    }

    const out = new Float32Array(n);
    if (voxelCount === 0) return { box, data: out, voxelCount };

    if (params.extend <= 0 && params.softEdge <= 0) {
        for (let i = 0; i < n; i++) out[i] = inside[i];
        return { box, data: out, voxelCount };
    }

    await ctx.update({ message: 'Computing distance transform…' });
    const distSq = squaredDistanceTransform3D(inside, bx, by, bz, out);

    await ctx.update({ message: 'Applying soft edge…' });
    for (let i = 0; i < n; i++) {
        out[i] = softValue(Math.sqrt(distSq[i]), params.extend, params.softEdge);
    }
    return { box, data: out, voxelCount };
}

function computeStats(data: Float32Array) {
    let min = Infinity, max = -Infinity, sum = 0;
    for (let i = 0; i < data.length; i++) {
        const v = data[i];
        if (v < min) min = v;
        if (v > max) max = v;
        sum += v;
    }
    const mean = data.length ? sum / data.length : 0;
    let sumSq = 0;
    for (let i = 0; i < data.length; i++) sumSq += (data[i] - mean) ** 2;
    return { min, max, mean, sigma: data.length ? Math.sqrt(sumSq / data.length) : 0 };
}

/**
 * Wraps a body mask result into a `Volume` covering only its bounding box. The grid transform
 * is the source transform composed with the box offset, so the cropped volume renders in place.
 */
export function buildCroppedMaskVolume(source: Volume, result: BodyMaskResult, label: string): Volume {
    const g2c = Grid.getGridToCartesianTransform(source.grid);
    const offset = Mat4.fromTranslation(Mat4(), result.box.min);
    const matrix = Mat4.mul(Mat4(), g2c, offset);
    const [bx, by, bz] = result.box.dims;
    const space = Tensor.Space([bx, by, bz], [2, 1, 0], Float32Array);

    return {
        label,
        entryId: source.entryId,
        grid: {
            transform: { kind: 'matrix', matrix },
            cells: Tensor.create(space, Tensor.Data1(result.data)),
            stats: computeStats(result.data),
        },
        instances: source.instances,
        sourceData: { kind: 'custom', name: 'Volume Body Mask', data: null } as any,
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
        _localPropertyData: Object.create(null),
    };
}

/** Writes a bounding-box mask result into a full-size array laid out like the source grid. */
export function scatterToFullBox(result: BodyMaskResult, space: Tensor.Space, out: Float32Array): Float32Array {
    const [bx, by, bz] = result.box.dims;
    const [x0, y0, z0] = result.box.min;
    for (let z = 0; z < bz; z++) {
        for (let y = 0; y < by; y++) {
            const row = y * bx + z * bx * by;
            for (let x = 0; x < bx; x++) {
                out[space.dataOffset(x0 + x, y0 + y, z0 + z)] = result.data[row + x];
            }
        }
    }
    return out;
}
