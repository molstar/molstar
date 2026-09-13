/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Volume } from '../../mol-model/volume';

/**
 * Data offsets (memory order) of all voxels with density >= `thresholdAbs`.
 * Labelling operations iterate these instead of the whole grid.
 */
export function computeCandidates(volume: Volume, thresholdAbs: number): Int32Array {
    const data = volume.grid.cells.data as unknown as ArrayLike<number>;
    const n = data.length;

    let count = 0;
    for (let i = 0; i < n; i++) {
        if (data[i] >= thresholdAbs) count++;
    }

    const offsets = new Int32Array(count);
    for (let i = 0, j = 0; i < n; i++) {
        if (data[i] >= thresholdAbs) offsets[j++] = i;
    }
    return offsets;
}
