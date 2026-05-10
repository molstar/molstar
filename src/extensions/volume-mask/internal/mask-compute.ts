/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Camera } from '../../../mol-canvas3d/camera';
import { Viewport, cameraProject } from '../../../mol-canvas3d/camera/util';
import { Grid, Volume } from '../../../mol-model/volume';
import { Tensor, Mat4, Vec3, Vec4 } from '../../../mol-math/linear-algebra';
import { CustomProperties } from '../../../mol-model/custom-property';
import { RuntimeContext } from '../../../mol-task';
import { pointInPolygon2D as pointInPolygon } from '../../../mol-math/geometry/polygon';
import type { ViewMask } from '../types';

const tmpVec3 = Vec3();
const tmpVec4 = Vec4();
const tmpG2C = Mat4();

/** Pre-computed camera data for a single ViewMask — built once before the voxel loop. */
interface PreparedMask {
    projectionView: Mat4;
    viewport: Viewport;
    normPolygon: [number, number][];
    inverted: boolean;
}

function prepareMask(mask: ViewMask): PreparedMask {
    const viewport = Viewport.create(0, 0, mask.viewportWidth, mask.viewportHeight);
    const cam = new Camera(mask.cameraSnapshot, viewport);
    cam.update();
    // Deep-copy the matrix (Camera keeps it as a mutable property)
    const projectionView = Mat4.copy(Mat4(), cam.projectionView);
    const w = mask.canvasWidth, h = mask.canvasHeight;
    const normPolygon = mask.polygon.map(([px, py]) => [px / w, py / h] as [number, number]);
    return { projectionView, viewport, normPolygon, inverted: !!mask.inverted };
}

/** Projects a world-space point to normalised [0,1] canvas coords (Y flipped). */
const tmpNorm: [number, number] = [0, 0];
function projectToNormInPlace(worldPos: Vec3, prepared: PreparedMask): void {
    cameraProject(tmpVec4, worldPos, prepared.viewport, prepared.projectionView);
    tmpNorm[0] = tmpVec4[0] / prepared.viewport.width;
    tmpNorm[1] = 1 - tmpVec4[1] / prepared.viewport.height;
}

export interface MaskComputeParams {
    viewMasks: ViewMask[];
    threshold: Volume.IsoValue;
    dilation: number;
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
    const { viewMasks, threshold, dilation, skipThreshold, invertPolygons } = params;
    const { cells: { space, data }, stats } = volume.grid;
    const [nx, ny, nz] = space.dimensions as [number, number, number];
    const selected = new Uint8Array(nx * ny * nz);

    const absThreshold = skipThreshold ? -Infinity : Volume.IsoValue.toAbsolute(threshold, stats).absoluteValue;
    Mat4.copy(tmpG2C, Grid.getGridToCartesianTransform(volume.grid));

    // Pre-compute camera matrices and normalised polygons ONCE per mask
    const prepared = viewMasks.map(prepareMask);
    const step = Math.max(1, Math.floor(nx / 20));

    for (let i = 0; i < nx; i++) {
        if (i % step === 0) {
            await ctx.update({ message: 'Selecting voxels…', current: i, max: nx });
        }
        for (let j = 0; j < ny; j++) {
            for (let k = 0; k < nz; k++) {
                if (space.get(data, i, j, k) < absThreshold) continue;

                Vec3.set(tmpVec3, i + 0.5, j + 0.5, k + 0.5);
                Vec3.transformMat4(tmpVec3, tmpVec3, tmpG2C);

                let passes = true;
                for (let m = 0; m < prepared.length; m++) {
                    projectToNormInPlace(tmpVec3, prepared[m]);
                    const inside = pointInPolygon(tmpNorm[0], tmpNorm[1], prepared[m].normPolygon);
                    if (prepared[m].inverted ? inside : !inside) {
                        passes = false;
                        break;
                    }
                }

                if (invertPolygons ? !passes : passes) selected[space.dataOffset(i, j, k)] = 1;
            }
        }
    }

    if (dilation <= 0) return selected;
    await ctx.update({ message: 'Dilating mask…' });
    return dilate3D(selected, nx, ny, nz, dilation, space);
}

export function dilate3D(
    src: Uint8Array,
    nx: number, ny: number, nz: number,
    radius: number,
    space: Grid['cells']['space']
): Uint8Array {
    const dst = new Uint8Array(src.length);
    const r2 = radius * radius;
    for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) {
            for (let k = 0; k < nz; k++) {
                if (src[space.dataOffset(i, j, k)] === 0) continue;
                for (let di = -radius; di <= radius; di++) {
                    const ni = i + di; if (ni < 0 || ni >= nx) continue;
                    for (let dj = -radius; dj <= radius; dj++) {
                        const nj = j + dj; if (nj < 0 || nj >= ny) continue;
                        for (let dk = -radius; dk <= radius; dk++) {
                            if (di * di + dj * dj + dk * dk > r2) continue;
                            const nk = k + dk; if (nk < 0 || nk >= nz) continue;
                            dst[space.dataOffset(ni, nj, nk)] = 1;
                        }
                    }
                }
            }
        }
    }
    return dst;
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

/**
 * Cosine soft edge: for each outside voxel, distance to the nearest inside voxel is
 * computed via a 2-pass chamfer approximation, then mapped through
 *   value = 0.5 + 0.5 * cos(π * min(dist / (width+1), 1))
 * giving 1.0 inside, a smooth falloff in the transition zone, and 0.0 far outside.
 */
export function computeSoftEdge(
    maskData: Uint8Array,
    nx: number, ny: number, nz: number,
    softEdgeWidth: number,
    space: Grid['cells']['space']
): Float32Array {
    const n = nx * ny * nz;
    const W = softEdgeWidth + 1;
    const dist = new Float32Array(n).fill(W + 1);

    for (let i = 0; i < nx; i++)
        for (let j = 0; j < ny; j++)
            for (let k = 0; k < nz; k++)
                if (maskData[space.dataOffset(i, j, k)] > 0)
                    dist[space.dataOffset(i, j, k)] = 0;

    // Forward pass
    for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) {
            for (let k = 0; k < nz; k++) {
                const idx = space.dataOffset(i, j, k);
                if (dist[idx] === 0) continue;
                let d = dist[idx];
                if (i > 0) d = Math.min(d, dist[space.dataOffset(i - 1, j, k)] + 1);
                if (j > 0) d = Math.min(d, dist[space.dataOffset(i, j - 1, k)] + 1);
                if (k > 0) d = Math.min(d, dist[space.dataOffset(i, j, k - 1)] + 1);
                if (i > 0 && j > 0) d = Math.min(d, dist[space.dataOffset(i - 1, j - 1, k)] + 1.414);
                if (i > 0 && k > 0) d = Math.min(d, dist[space.dataOffset(i - 1, j, k - 1)] + 1.414);
                if (j > 0 && k > 0) d = Math.min(d, dist[space.dataOffset(i, j - 1, k - 1)] + 1.414);
                if (i > 0 && j > 0 && k > 0) d = Math.min(d, dist[space.dataOffset(i - 1, j - 1, k - 1)] + 1.732);
                dist[idx] = d;
            }
        }
    }

    // Backward pass
    for (let i = nx - 1; i >= 0; i--) {
        for (let j = ny - 1; j >= 0; j--) {
            for (let k = nz - 1; k >= 0; k--) {
                const idx = space.dataOffset(i, j, k);
                if (dist[idx] === 0) continue;
                let d = dist[idx];
                if (i < nx - 1) d = Math.min(d, dist[space.dataOffset(i + 1, j, k)] + 1);
                if (j < ny - 1) d = Math.min(d, dist[space.dataOffset(i, j + 1, k)] + 1);
                if (k < nz - 1) d = Math.min(d, dist[space.dataOffset(i, j, k + 1)] + 1);
                if (i < nx - 1 && j < ny - 1) d = Math.min(d, dist[space.dataOffset(i + 1, j + 1, k)] + 1.414);
                if (i < nx - 1 && k < nz - 1) d = Math.min(d, dist[space.dataOffset(i + 1, j, k + 1)] + 1.414);
                if (j < ny - 1 && k < nz - 1) d = Math.min(d, dist[space.dataOffset(i, j + 1, k + 1)] + 1.414);
                if (i < nx - 1 && j < ny - 1 && k < nz - 1) d = Math.min(d, dist[space.dataOffset(i + 1, j + 1, k + 1)] + 1.732);
                dist[idx] = d;
            }
        }
    }

    const result = new Float32Array(n);
    for (let i = 0; i < n; i++) {
        result[i] = 0.5 + 0.5 * Math.cos(Math.PI * Math.min(dist[i] / W, 1));
    }
    return result;
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
