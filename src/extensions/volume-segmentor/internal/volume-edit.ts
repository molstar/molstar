/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Volume } from '../../../mol-model/volume';
import { flipVolumeX, removeDust } from '../../volume-mask/internal/volume-ops';
import { BodyLabels } from '../labels';

type Cells = { [i: number]: number, length: number };

/** Recomputes min/max/mean/sigma of a volume whose voxel data was edited in place. */
export function recomputeStats(volume: Volume) {
    const cells = volume.grid.cells.data as unknown as Cells;
    let min = Infinity, max = -Infinity, sum = 0;
    for (let i = 0; i < cells.length; i++) {
        const v = cells[i];
        if (v < min) min = v;
        if (v > max) max = v;
        sum += v;
    }
    const mean = cells.length ? sum / cells.length : 0;
    let sumSq = 0;
    for (let i = 0; i < cells.length; i++) sumSq += (cells[i] - mean) ** 2;
    const stats = volume.grid.stats as { min: number, max: number, mean: number, sigma: number };
    stats.min = min;
    stats.max = max;
    stats.mean = mean;
    stats.sigma = cells.length ? Math.sqrt(sumSq / cells.length) : 0;
}

/**
 * Drops everything cached on the volume (custom properties, GPU textures, derived data) after
 * its voxel data changed, keeping only the body labels.
 */
export function clearVolumeCaches(volume: Volume) {
    const labels = BodyLabels.get(volume);
    volume.customProperties.dispose();
    for (const k of Object.keys(volume._propertyData)) delete volume._propertyData[k];
    for (const k of Object.keys(volume._localPropertyData)) delete volume._localPropertyData[k];
    if (labels) BodyLabels.set(volume, labels);
}

/**
 * Zeroes face-connected components of voxels >= `thresholdAbs` that are smaller than
 * `minVoxels`, editing the volume data in place; refreshes its stats and drops cached derived
 * data (labels are kept). Returns the number of voxels zeroed.
 */
export function removeDustInPlace(volume: Volume, minVoxels: number, thresholdAbs: number): number {
    const { cells } = volume.grid;
    const [nx, ny, nz] = cells.space.dimensions as [number, number, number];
    const zeroed = removeDust(cells.data as unknown as Cells, nx, ny, nz, cells.space, minVoxels, thresholdAbs);
    if (zeroed > 0) {
        recomputeStats(volume);
        clearVolumeCaches(volume);
    }
    return zeroed;
}

/**
 * Mirrors the volume along X, editing its data in place, to fix a map stored with the opposite
 * handedness; drops cached derived data (labels are kept). Mirroring only moves voxels, so the
 * grid stats stay valid.
 */
export function flipHandednessInPlace(volume: Volume) {
    const { cells } = volume.grid;
    const [nx, ny, nz] = cells.space.dimensions as [number, number, number];
    flipVolumeX(cells.data as unknown as Cells, nx, ny, nz, cells.space);
    clearVolumeCaches(volume);
}
