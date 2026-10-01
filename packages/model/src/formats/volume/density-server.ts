/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { DensityServer_Data_Database } from '@molstar/io/reader/cif/schema/density-server';
import type { Volume } from '@molstar/model/model/volume';
import { Task } from '@molstar/core/task';
import { SpacegroupCell, Box3D } from '@molstar/core/math/geometry';
import { Mat4, Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import type { ModelFormat } from '../format.js';
import { CustomProperties } from '@molstar/model/model/custom-property';

export function volumeFromDensityServerData(source: DensityServer_Data_Database, params?: Partial<{ label: string, entryId: string }>): Task<Volume> {
    return Task.create<Volume>('Create Volume', async ctx => {
        const { volume_data_3d_info: info, volume_data_3d: values } = source;
        const cell = SpacegroupCell.create(
            info.spacegroup_number.value(0) || 'P 1',
            Vec3.ofArray(info.spacegroup_cell_size.value(0)),
            Vec3.scale(Vec3.zero(), Vec3.ofArray(info.spacegroup_cell_angles.value(0)), Math.PI / 180)
        );

        const axis_order_fast_to_slow = info.axis_order.value(0);

        const normalizeOrder = Tensor.convertToCanonicalAxisIndicesFastToSlow(axis_order_fast_to_slow);

        // sample count is in "axis order" and needs to be reordered
        const sample_count = normalizeOrder(info.sample_count.value(0));
        const tensorSpace = Tensor.Space(sample_count, Tensor.invertAxisOrder(axis_order_fast_to_slow), Float32Array);

        const data = Tensor.create(tensorSpace, Tensor.Data1(values.values.toArray({ array: Float32Array })));

        // origin and dimensions are in "axis order" and need to be reordered
        const origin = Vec3.ofArray(normalizeOrder(info.origin.value(0)));
        const dimensions = Vec3.ofArray(normalizeOrder(info.dimensions.value(0)));

        return {
            label: params?.label,
            entryId: params?.entryId,
            grid: {
                transform: { kind: 'spacegroup', cell, fractionalBox: Box3D.create(origin, Vec3.add(Vec3.zero(), origin, dimensions)) },
                cells: data,
                stats: {
                    min: info.min_sampled.value(0),
                    max: info.max_sampled.value(0),
                    mean: info.mean_sampled.value(0),
                    sigma: info.sigma_sampled.value(0)
                },
                periodicity: Vec3.isInteger(dimensions) ? 'xyz' : 'none',
            },
            instances: [{ transform: Mat4.identity() }],
            sourceData: DscifFormat.create(source),
            customProperties: new CustomProperties(),
            _propertyData: Object.create(null),
            _localPropertyData: Object.create(null),
        };
    });
}

//

export { DscifFormat };

type DscifFormat = ModelFormat<DensityServer_Data_Database>

namespace DscifFormat {
    export function is(x?: ModelFormat): x is DscifFormat {
        return x?.kind === 'dscif';
    }

    export function create(dscif: DensityServer_Data_Database): DscifFormat {
        return { kind: 'dscif', name: dscif._name, data: dscif };
    }
}