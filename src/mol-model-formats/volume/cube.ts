/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CubeFile } from '../../mol-io/reader/cube/parser';
import { Mat4, Tensor } from '../../mol-math/linear-algebra';
import { Volume } from '../../mol-model/volume';
import { Task } from '../../mol-task';
import { arrayMax, arrayMean, arrayMin, arrayRms } from '../../mol-util/array';
import { ModelFormat } from '../format';
import { CustomProperties } from '../../mol-model/custom-property';
import { clamp } from '../../mol-math/interpolate';

export function volumeFromCube(source: CubeFile, params?: { dataIndex?: number, label?: string, entryId?: string, clamp?: { min: number, max: number } }): Task<Volume> {
    return Task.create<Volume>('Create Volume', async () => {
        const { header, values: sourceValues } = source;
        const space = Tensor.Space(header.dim, [0, 1, 2], Float64Array);

        let values: Float64Array;
        if (header.dataSetIds.length === 0) {
            values = sourceValues;
        } else {
            // get every nth value from the source values
            const [h, k, l] = header.dim;
            const nth = (params?.dataIndex || 0) + 1;
            const { min, max } = params?.clamp || { min: -Infinity, max: Infinity };

            let o = 0, s = 0;

            values = new Float64Array(h * k * l);
            for (let u = 0; u < h; u++) {
                for (let v = 0; v < k; v++) {
                    for (let w = 0; w < l; w++) {
                        values[o++] = clamp(sourceValues[s], min, max);
                        s += nth;
                    }
                }
            }
        }

        const data = Tensor.create(space, Tensor.Data1(values));

        const matrix = Mat4.fromTranslation(Mat4(), header.origin);
        const basis = Mat4.fromBasis(Mat4(), header.basisX, header.basisY, header.basisZ);
        Mat4.mul(matrix, matrix, basis);

        return {
            label: params?.label,
            entryId: params?.entryId,
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
            sourceData: CubeFormat.create(source),
            customProperties: new CustomProperties(),
            _propertyData: Object.create(null),
        };
    });
}

//

export { CubeFormat };

type CubeFormat = ModelFormat<CubeFile>

namespace CubeFormat {
    export function is(x?: ModelFormat): x is CubeFormat {
        return x?.kind === 'cube';
    }

    export function create(cube: CubeFile): CubeFormat {
        return { kind: 'cube', name: cube.name, data: cube };
    }
}