/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { VolumeData } from '../../mol-model/volume/data';
import { Task } from '../../mol-task';
import { SpacegroupCell, Box3D } from '../../mol-math/geometry';
import { Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { arrayMin, arrayMax, arrayMean, arrayRms } from '../../mol-util/array';
import { CubeFile } from '../../mol-io/reader/cube/parser';

export function volumeFromCube(source: CubeFile, params?: { dataIndex?: number }): Task<VolumeData> {
    return Task.create<VolumeData>('Create Volume Data', async () => {
        // TODO: support non-orthogonal axes

        const { header, values: sourceValues } = source;
        const angles = SpacegroupCell.Zero.anglesInRadians;
        const size = Vec3.create(header.basisX[0] * header.dim[0], header.basisY[1] * header.dim[1], header.basisZ[2] * header.dim[2]);
        const cell = SpacegroupCell.create('P 1', size, angles);

        if (header.basisX[1] !== 0 || header.basisX[2] !== 0
            || header.basisY[0] !== 0 || header.basisY[2] !== 0
            || header.basisZ[0] !== 0 || header.basisZ[1] !== 0) {
            throw new Error('Non-orthogonal bases not supported. (TODO)');
        }

        const origin_frac = Vec3.div(Vec3(), header.origin, size);
        const dimensions_frac = Vec3.add(Vec3(), origin_frac, Vec3.create(1, 1, 1));

        const space = Tensor.Space(header.dim, [0, 1, 2], Float64Array);

        let values: Float64Array;
        if (header.dataSetIds.length === 0) {
            values = sourceValues;
        } else {
            // get every nth value from the source values
            const [h, k, l] = header.dim;
            const nth = (params?.dataIndex || 0) + 1;

            let o = 0, s = 0;

            values = new Float64Array(h * k * l);
            for (let u = 0; u < h; u++) {
                for (let v = 0; v < k; v++) {
                    for (let w = 0; w < l; w++) {
                        values[o++] = sourceValues[s];
                        s += nth;
                    }
                }
            }
        }

        const data = Tensor.create(space, Tensor.Data1(values));

        return {
            cell,
            fractionalBox: Box3D.create(origin_frac, dimensions_frac),
            data,
            dataStats: {
                min: arrayMin(values),
                max: arrayMax(values),
                mean: arrayMean(values),
                sigma: arrayRms(values)
            }
        };
    });
}