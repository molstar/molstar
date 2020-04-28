/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Volume } from '../../mol-model/volume';
import { Task } from '../../mol-task';
import { SpacegroupCell, Box3D } from '../../mol-math/geometry';
import { Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { degToRad } from '../../mol-math/misc';
import { Dsn6File } from '../../mol-io/reader/dsn6/schema';
import { arrayMin, arrayMax, arrayMean, arrayRms } from '../../mol-util/array';
import { ModelFormat } from '../format';
import { CustomProperties } from '../../mol-model/custom-property';

export function volumeFromDsn6(source: Dsn6File, params?: { voxelSize?: Vec3, label?: string }): Task<Volume> {
    return Task.create<Volume>('Create Volume', async ctx => {
        const { header, values } = source;
        const size = Vec3.create(header.xlen, header.ylen, header.zlen);
        if (params && params.voxelSize) Vec3.mul(size, size, params.voxelSize);
        const angles = Vec3.create(degToRad(header.alpha), degToRad(header.beta), degToRad(header.gamma));
        const cell = SpacegroupCell.create('P 1', size, angles);

        const grid = [header.xRate, header.yRate, header.zRate];
        const extent = [header.xExtent, header.yExtent, header.zExtent];

        const gridOrigin = [header.xStart, header.yStart, header.zStart];

        const origin_frac = Vec3.create(gridOrigin[0] / grid[0], gridOrigin[1] / grid[1], gridOrigin[2] / grid[2]);
        const dimensions_frac = Vec3.create(extent[0] / grid[0], extent[1] / grid[1], extent[2] / grid[2]);

        const space = Tensor.Space(extent, [0, 1, 2], Float32Array);
        const data = Tensor.create(space, Tensor.Data1(values));

        return {
            label: params?.label,
            grid: {
                transform: { kind: 'spacegroup', cell, fractionalBox: Box3D.create(origin_frac, Vec3.add(Vec3.zero(), origin_frac, dimensions_frac)) },
                cells: data,
                stats: {
                    min: arrayMin(values),
                    max: arrayMax(values),
                    mean: arrayMean(values),
                    sigma: header.sigma !== undefined ? header.sigma : arrayRms(values)
                },
            },
            sourceData: Dsn6Format.create(source),
            customProperties: new CustomProperties(),
            _propertyData: Object.create(null),
        };
    });
}

//

export { Dsn6Format };

type Dsn6Format = ModelFormat<Dsn6File>

namespace Dsn6Format {
    export function is(x: ModelFormat): x is Dsn6Format {
        return x.kind === 'dsn6';
    }

    export function create(dsn6: Dsn6File): Dsn6Format {
        return { kind: 'dsn6', name: dsn6.name, data: dsn6 };
    }
}