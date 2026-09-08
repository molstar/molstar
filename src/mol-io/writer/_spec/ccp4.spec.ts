/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Mat4, Tensor } from '../../../mol-math/linear-algebra';
import { volumeFromCcp4 } from '../../../mol-model-formats/volume/ccp4';
import { Grid } from '../../../mol-model/volume';
import { parse } from '../../reader/ccp4/parser';
import { CCP4Writer } from '../ccp4/ccp4';

/** Deliberately not a cube, so a transposed axis would not go unnoticed. */
const Dimensions: [number, number, number] = [4, 3, 2];

/** Encodes the grid coords into the value, making every voxel tell where it belongs. */
const density = (i: number, j: number, k: number) => i + 10 * j + 100 * k;

/** Every axis order `Tensor` supports for rank 3. */
const AxisOrders = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

function testGrid(axisOrderSlowToFast: number[]) {
    const [nx, ny, nz] = Dimensions;
    const space = Tensor.Space(Dimensions, axisOrderSlowToFast, Float32Array);
    const data = new Float32Array(nx * ny * nz);
    for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) {
            for (let k = 0; k < nz; k++) data[space.dataOffset(i, j, k)] = density(i, j, k);
        }
    }

    let min = Infinity, max = -Infinity, sum = 0;
    for (let o = 0; o < data.length; o++) {
        min = Math.min(min, data[o]);
        max = Math.max(max, data[o]);
        sum += data[o];
    }
    const mean = sum / data.length;
    let sq = 0;
    for (let o = 0; o < data.length; o++) sq += (data[o] - mean) ** 2;

    const grid: Grid = {
        transform: { kind: 'matrix', matrix: Mat4.identity() },
        cells: Tensor.create(space, Tensor.Data1(data)),
        stats: { min, max, mean, sigma: Math.sqrt(sq / data.length) },
    };
    return { grid, data };
}

describe('CCP4Writer', () => {
    for (const order of AxisOrders) {
        it(`round trips voxel values for axis order ${order}`, async () => {
            const { grid, data } = testGrid(order);
            const buffer = CCP4Writer.writeMrc(grid, data);
            expect(buffer.byteLength).toBe(1024 + 4 * 3 * 2 * 4);

            const parsed = await parse(new Uint8Array(buffer), 'test.mrc').run();
            if (parsed.isError) throw new Error(parsed.message);

            const { header } = parsed.result;
            expect([header.NC, header.NR, header.NS]).toEqual(Dimensions);
            expect(header.MODE).toBe(2);
            expect([header.MAPC, header.MAPR, header.MAPS]).toEqual([1, 2, 3]);
            expect(header.AMIN).toBe(0);
            expect(header.AMAX).toBe(123);

            const volume = await volumeFromCcp4(parsed.result).run();
            const { space, data: values } = volume.grid.cells;
            expect(Array.from(space.dimensions)).toEqual(Dimensions);
            for (let i = 0; i < Dimensions[0]; i++) {
                for (let j = 0; j < Dimensions[1]; j++) {
                    for (let k = 0; k < Dimensions[2]; k++) {
                        expect(space.get(values, i, j, k)).toBe(density(i, j, k));
                    }
                }
            }
        });
    }
});
