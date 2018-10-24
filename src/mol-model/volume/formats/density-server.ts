/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DensityServer_Data_Database } from 'mol-io/reader/cif/schema/density-server'
import { VolumeData } from '../data'
import { Task } from 'mol-task';
import { SpacegroupCell, Box3D } from 'mol-math/geometry';
import { Tensor, Vec3 } from 'mol-math/linear-algebra';

function parseDensityServerData(source: DensityServer_Data_Database): Task<VolumeData> {
    return Task.create<VolumeData>('Parse Volume Data', async ctx => {
        const { volume_data_3d_info: info, volume_data_3d: values } = source;
        const cell = SpacegroupCell.create(
            info.spacegroup_number.value(0),
            Vec3.ofArray(info.spacegroup_cell_size.value(0)),
            Vec3.scale(Vec3.zero(), Vec3.ofArray(info.spacegroup_cell_angles.value(0)), Math.PI / 180)
        );

        const tensorSpace = Tensor.Space(info.sample_count.value(0), info.axis_order.value(0), Float32Array);
        const data = Tensor.create(tensorSpace, Tensor.Data1(values.values.toArray()));

        const origin = Vec3.ofArray(info.origin.value(0))
        const dimensions = Vec3.ofArray(info.dimensions.value(0));

        return {
            cell,
            fractionalBox: Box3D.create(origin, Vec3.add(Vec3.zero(), origin, dimensions)),
            data,
            dataStats: {
                min: info.min_sampled.value(0),
                max: info.max_sampled.value(0),
                mean: info.mean_sampled.value(0),
                sigma: info.sigma_sampled.value(0)
            }
        };
    });
}

export { parseDensityServerData }