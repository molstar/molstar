/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DxFile } from '../../mol-io/reader/dx/parser';
import { Mat4, Tensor } from '../../mol-math/linear-algebra';
import { VolumeData } from '../../mol-model/volume/data';
import { Task } from '../../mol-task';
import { arrayMax, arrayMean, arrayMin, arrayRms } from '../../mol-util/array';

export function volumeFromDx(source: DxFile, params?: { label?: string }): Task<VolumeData> {
    return Task.create<VolumeData>('Create Volume Data', async () => {
        const { header, values } = source;
        const space = Tensor.Space(header.dim, [0, 1, 2], Float64Array);
        const data = Tensor.create(space, Tensor.Data1(values));
        const matrix = Mat4.fromTranslation(Mat4(), header.min);
        const basis = Mat4.fromScaling(Mat4(), header.h);
        Mat4.mul(matrix, matrix, basis);

        return {
            label: params?.label,
            transform: { kind: 'matrix', matrix },
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