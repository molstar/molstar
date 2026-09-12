/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Grid, Volume } from '../../../../mol-model/volume';
import { Tensor, Mat4, Vec3 } from '../../../../mol-math/linear-algebra';
import { CustomProperties } from '../../../../mol-model/custom-property';
import { RuntimeContext } from '../../../../mol-task';
import { computeCandidates } from '../../candidates';
import { prepareMask } from '../../view-projection';
import { passesAllViews, selectByViews } from '../../view-selection';
import type { ViewMask } from '../types';

const tmpVec3 = Vec3();
const tmpG2C = Mat4();

const tmpNorm: [number, number] = [0, 0];

export interface MaskComputeParams {
    viewMasks: ViewMask[];
    threshold: Volume.IsoValue;
    /** When true, skip the density threshold check (use polygons as purely spatial filter). */
    skipThreshold?: boolean;
    /** When true, select voxels that are above threshold but OUTSIDE the polygon filters. */
    invertPolygons?: boolean;
}

export async function computeVolumeMask(
    volume: Volume,
    params: MaskComputeParams,
    ctx: RuntimeContext
): Promise<Uint8Array> {
    const { viewMasks, threshold, skipThreshold, invertPolygons } = params;
    const { cells: { space }, stats } = volume.grid;
    const [nx, ny, nz] = space.dimensions as [number, number, number];
    const selected = new Uint8Array(nx * ny * nz);

    if (skipThreshold) {
        // No density gate, so every voxel is a candidate and the grid has to be walked.
        Mat4.copy(tmpG2C, Grid.getGridToCartesianTransform(volume.grid));
        const prepared = viewMasks.map(prepareMask);
        const step = Math.max(1, Math.floor(nx / 20));
        for (let i = 0; i < nx; i++) {
            if (i % step === 0) await ctx.update({ message: 'Selecting voxels…', current: i, max: nx });
            for (let j = 0; j < ny; j++) {
                for (let k = 0; k < nz; k++) {
                    Vec3.set(tmpVec3, i, j, k);
                    Vec3.transformMat4(tmpVec3, tmpVec3, tmpG2C);
                    const passes = passesAllViews(tmpVec3, prepared, tmpNorm);
                    if (invertPolygons ? !passes : passes) selected[space.dataOffset(i, j, k)] = 1;
                }
            }
        }
    } else {
        const absThreshold = Volume.IsoValue.toAbsolute(threshold, stats).absoluteValue;
        await ctx.update({ message: 'Collecting voxels above the threshold…' });
        const candidates = computeCandidates(volume, absThreshold);
        await selectByViews(volume, candidates, viewMasks, selected, ctx, !!invertPolygons);
    }

    return selected;
}


/**
 * Mark all voxels within `radiusAngstrom` Å of any atom in `atomPositions`.
 * `atomPositions` is a flat Float32Array of world-space XYZ triplets.
 */
export function computeStructureMask(
    volume: Volume,
    atomPositions: Float32Array,
    radiusAngstrom: number
): Uint8Array {
    const { cells: { space } } = volume.grid;
    const [nx, ny, nz] = space.dimensions as [number, number, number];
    const selected = new Uint8Array(nx * ny * nz);

    const g2c = Grid.getGridToCartesianTransform(volume.grid);
    const c2g = Mat4.invert(Mat4(), g2c);

    // Estimate voxel size as average column magnitude of g2c
    const vx = Math.sqrt(g2c[0] ** 2 + g2c[1] ** 2 + g2c[2] ** 2);
    const vy = Math.sqrt(g2c[4] ** 2 + g2c[5] ** 2 + g2c[6] ** 2);
    const vz = Math.sqrt(g2c[8] ** 2 + g2c[9] ** 2 + g2c[10] ** 2);
    const avgVoxelSize = (vx + vy + vz) / 3;
    const r = Math.ceil(radiusAngstrom / avgVoxelSize);
    const r2 = r * r;

    const gp = Vec3();
    const atomCount = (atomPositions.length / 3) | 0;

    for (let a = 0; a < atomCount; a++) {
        Vec3.set(gp, atomPositions[a * 3], atomPositions[a * 3 + 1], atomPositions[a * 3 + 2]);
        Vec3.transformMat4(gp, gp, c2g);

        const ci = Math.round(gp[0]), cj = Math.round(gp[1]), ck = Math.round(gp[2]);

        for (let di = -r; di <= r; di++) {
            const ni = ci + di; if (ni < 0 || ni >= nx) continue;
            for (let dj = -r; dj <= r; dj++) {
                const nj = cj + dj; if (nj < 0 || nj >= ny) continue;
                for (let dk = -r; dk <= r; dk++) {
                    if (di * di + dj * dj + dk * dk > r2) continue;
                    const nk = ck + dk; if (nk < 0 || nk >= nz) continue;
                    selected[space.dataOffset(ni, nj, nk)] = 1;
                }
            }
        }
    }

    return selected;
}


function computeSigma(data: Float32Array, mean: number): number {
    let sumSq = 0;
    for (let i = 0; i < data.length; i++) sumSq += (data[i] - mean) ** 2;
    return Math.sqrt(sumSq / data.length);
}

export function buildMaskVolume(source: Volume, data: Uint8Array | Float32Array): Volume {
    const srcGrid = source.grid;
    const [nx, ny, nz] = srcGrid.cells.space.dimensions as [number, number, number];
    const axisOrder = srcGrid.cells.space.axisOrderSlowToFast.slice() as number[];
    const space = Tensor.Space([nx, ny, nz], axisOrder, Float32Array);
    const floatData = data instanceof Float32Array ? data : new Float32Array(data);
    let sum = 0;
    for (let i = 0; i < floatData.length; i++) sum += floatData[i];
    const mean = sum / floatData.length;

    return {
        label: 'Mask',
        entryId: source.entryId,
        grid: {
            transform: srcGrid.transform,
            cells: Tensor.create(space, Tensor.Data1(floatData)),
            stats: { min: 0, max: 1, mean, sigma: computeSigma(floatData, mean) },
        },
        instances: source.instances,
        sourceData: { kind: 'custom', name: 'Volume Mask', data: null } as any,
        customProperties: new CustomProperties(),
        _propertyData: {},
        _localPropertyData: {},
    };
}
