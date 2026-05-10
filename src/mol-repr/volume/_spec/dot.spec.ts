/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { CustomProperties } from '../../../mol-model/custom-property';
import { Grid, Volume } from '../../../mol-model/volume';
import { Mat4, Tensor } from '../../../mol-math/linear-algebra';
import { createVolumeSphereImpostor } from '../dot';

function createTestVolume(dimensions: [number, number, number], data: number[]): Volume {
    return {
        grid: {
            transform: { kind: 'matrix', matrix: Mat4.identity() },
            cells: Tensor.create(Tensor.Space(dimensions, [2, 1, 0]), Tensor.Data1(data)),
            stats: { min: 0, max: 1, mean: 0.5, sigma: 0.5 },
        } satisfies Grid,
        instances: [{ transform: Mat4.identity() }],
        sourceData: { kind: 'test', name: 'test', data: {} } as any,
        customProperties: new CustomProperties(),
        _propertyData: Object.create(null),
        _localPropertyData: Object.create(null),
    };
}

describe('volume dot representation', () => {
    it('adds sphere impostor dots in Morton order for LOD sampling', () => {
        const volume = createTestVolume([2, 2, 2], [
            1, 1,
            1, 1,
            1, 1,
            1, 1,
        ]);
        const spheres = createVolumeSphereImpostor(undefined as any, volume, 0, undefined as any, {
            isoValue: Volume.IsoValue.absolute(0.5),
            perturbPositions: false,
            lodLevels: [{ minDistance: 0, maxDistance: 0, overlap: 0, stride: 0, scaleBias: 3 }],
        } as any);

        expect(Array.from(spheres.groupBuffer.ref.value)).toEqual([0, 4, 2, 6, 1, 5, 3, 7]);
    });

    it('adds sphere impostor dots in row-major order when no LOD levels are configured', () => {
        const volume = createTestVolume([2, 2, 2], [
            1, 1,
            1, 1,
            1, 1,
            1, 1,
        ]);
        const spheres = createVolumeSphereImpostor(undefined as any, volume, 0, undefined as any, {
            isoValue: Volume.IsoValue.absolute(0.5),
            perturbPositions: false,
            lodLevels: [],
        } as any);

        expect(Array.from(spheres.groupBuffer.ref.value)).toEqual([0, 1, 2, 3, 4, 5, 6, 7]);
    });
});
