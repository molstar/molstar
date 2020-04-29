/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DxFile } from '../../mol-io/reader/dx/parser';
import { Mat4, Tensor } from '../../mol-math/linear-algebra';
import { Volume } from '../../mol-model/volume';
import { Task } from '../../mol-task';
import { arrayMax, arrayMean, arrayMin, arrayRms } from '../../mol-util/array';
import { ModelFormat } from '../format';
import { CustomProperties } from '../../mol-model/custom-property';

export function volumeFromDx(source: DxFile, params?: { label?: string }): Task<Volume> {
    return Task.create<Volume>('Create Volume', async () => {
        const { header, values } = source;
        const space = Tensor.Space(header.dim, [0, 1, 2], Float64Array);
        const data = Tensor.create(space, Tensor.Data1(values));
        const matrix = Mat4.fromTranslation(Mat4(), header.min);
        const basis = Mat4.fromScaling(Mat4(), header.h);
        Mat4.mul(matrix, matrix, basis);

        return {
            label: params?.label,
            grid: {
                transform: { kind: 'matrix', matrix },
                cells: data,
                stats: {
                    min: arrayMin(values),
                    max: arrayMax(values),
                    mean: arrayMean(values),
                    sigma: arrayRms(values)
                },
            },
            sourceData: DxFormat.create(source),
            customProperties: new CustomProperties(),
            _propertyData: Object.create(null),
        };
    });
}

//

export { DxFormat };

type DxFormat = ModelFormat<DxFile>

namespace DxFormat {
    export function is(x: ModelFormat): x is DxFormat {
        return x.kind === 'dx';
    }

    export function create(dx: DxFile): DxFormat {
        return { kind: 'dx', name: dx.name, data: dx };
    }
}