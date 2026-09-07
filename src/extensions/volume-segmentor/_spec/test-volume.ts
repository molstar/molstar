/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { CustomProperties } from '../../../mol-model/custom-property';
import { Grid, Volume } from '../../../mol-model/volume';
import { Mat4, Tensor } from '../../../mol-math/linear-algebra';

export const CanonicalOrder = [2, 1, 0];
export const SwappedOrder = [0, 1, 2];

/** Synthetic volume; `fill(i, j, k)` gives the density at logical grid coords. */
export function createTestVolume(dimensions: [number, number, number], fill: (i: number, j: number, k: number) => number, axisOrderSlowToFast: number[] = CanonicalOrder, matrix: Mat4 = Mat4.identity()): Volume {
    const space = Tensor.Space(dimensions, axisOrderSlowToFast, Float32Array);
    const data = new Float32Array(dimensions[0] * dimensions[1] * dimensions[2]);
    for (let i = 0; i < dimensions[0]; i++) {
        for (let j = 0; j < dimensions[1]; j++) {
            for (let k = 0; k < dimensions[2]; k++) {
                data[space.dataOffset(i, j, k)] = fill(i, j, k);
            }
        }
    }
    let min = Infinity, max = -Infinity, sum = 0;
    for (let i = 0; i < data.length; i++) { min = Math.min(min, data[i]); max = Math.max(max, data[i]); sum += data[i]; }
    const mean = sum / data.length;
    let sq = 0;
    for (let i = 0; i < data.length; i++) sq += (data[i] - mean) ** 2;

    return {
        label: 'test.mrc',
        grid: {
            transform: { kind: 'matrix', matrix },
            cells: Tensor.create(space, Tensor.Data1(data)),
            stats: { min, max, mean, sigma: Math.sqrt(sq / data.length) },
        } satisfies Grid,
        instances: [{ transform: Mat4.identity() }],
        sourceData: { kind: 'test', name: 'test', data: {} } as any,
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
        _localPropertyData: Object.create(null),
    };
}

export function offsetOf(volume: Volume, i: number, j: number, k: number) {
    return volume.grid.cells.space.dataOffset(i, j, k);
}
