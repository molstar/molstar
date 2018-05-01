/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { VolumeData, VolumeIsoValue } from 'mol-model/volume'
import { Task } from 'mol-task'
import { Mesh } from '../../shape/mesh';
import { computeMarchingCubes } from '../../util/marching-cubes/algorithm';

function computeVolumeMesh(volume: VolumeData, isoValue: VolumeIsoValue) {
    return Task.create<Mesh>('Volume Surface', async ctx => {
        ctx.update({ message: 'Marching cubes...' });

        const mesh = await ctx.runChild(computeMarchingCubes({
            isoLevel: VolumeIsoValue.toAbsolute(isoValue).absoluteValue,
            scalarField: volume.data
        }));

        const transform = VolumeData.getGridToCartesianTransform(volume);
        ctx.update({ message: 'Transforming mesh...' });
        Mesh.transformImmediate(mesh, transform);

        return mesh;
    });
}

export { computeVolumeMesh }