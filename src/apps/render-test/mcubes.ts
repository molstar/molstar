/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Run } from 'mol-task'
import { compute } from 'mol-geo/util/marching-cubes/algorithm'
import { Surface } from 'mol-geo/shape/surface'
import { Tensor } from 'mol-math/linear-algebra'

function fillField(tensor: Tensor, f: (x: number, y: number, z: number) => number, min: number[], max: number[]): Tensor {
    const { space: { set, dimensions: [ii, jj, kk] }, data } = tensor;
 
    const dx = (max[0] - min[0]) / ii;
    const dy = (max[1] - min[1]) / jj;
    const dz = (max[2] - min[2]) / kk;

    for (let i = 0, x = min[0]; i < ii; i++, x += dx) {
        for (let j = 0, y = min[1]; j < jj; j++, y += dy) {
            for (let k = 0, z = min[2]; k < kk; k++, z += dz) {
                set(data, i, j, k, f(x, y, z));
            }
        }
    }

    return tensor
}

export default async function computeSurface(f: (x: number, y: number, z: number) => number, data?: { field: Tensor, surface: Surface }) {
    let field: Tensor;
    if (data) field = data.field;
    else {
        const space = Tensor.Space([30, 30, 30], [0, 1, 2]);
        field = Tensor.create(space, space.create(Float32Array));
    }

    fillField(field, f, [-1.1, -1.1, -1.1], [1.1, 1.1, 1.1]);
    const surface = await Run(compute({
        scalarField: field,
        isoLevel: 0,
        buffers: data ? {
            vertex: data.surface.vertexBuffer,
            index: data.surface.indexBuffer,
            normal: data.surface.normalBuffer
        } : void 0
    }));
    return { surface, field };
}