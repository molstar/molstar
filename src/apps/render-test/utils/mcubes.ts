/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { computeMarchingCubes } from 'mol-geo/util/marching-cubes/algorithm'
import { Mesh } from 'mol-geo/shape/mesh'
import { Tensor, Mat4, Vec3 } from 'mol-math/linear-algebra'

function fillField(tensor: Tensor, f: (x: number, y: number, z: number) => number, min: number[], max: number[]): Tensor {
    const { space: { set, dimensions: [ii, jj, kk] }, data } = tensor;

    const dx = (max[0] - min[0]) / (ii - 1);
    const dy = (max[1] - min[1]) / (jj - 1);
    const dz = (max[2] - min[2]) / (kk - 1);

    for (let i = 0, x = min[0]; i < ii; i++, x += dx) {
        for (let j = 0, y = min[1]; j < jj; j++, y += dy) {
            for (let k = 0, z = min[2]; k < kk; k++, z += dz) {
                set(data, i, j, k, f(x, y, z));
            }
        }
    }

    return tensor
}

export default async function computeSurface(f: (x: number, y: number, z: number) => number, data?: { field: Tensor, surface: Mesh }) {
    let field: Tensor;
    if (data) field = data.field;
    else {
        const space = Tensor.Space([30, 30, 30], [0, 1, 2]);
        field = Tensor.create(space, space.create(Float32Array));
    }

    const min = Vec3.create(-1.1, -1.1, -1.1), max = Vec3.create(1.1, 1.1, 1.1);

    fillField(field, f, min, max);
    const surface = await computeMarchingCubes({
        scalarField: field,
        isoLevel: 0,
        oldSurface: data ? data.surface : void 0
    }).run();

    const translation = Mat4.fromTranslation(Mat4.zero(), min);
    const grid = Vec3.zero();
    Vec3.fromArray(grid, field.space.dimensions as any, 0);
    const size = Vec3.sub(Vec3.zero(), max, min);
    const scale = Mat4.fromScaling(Mat4.zero(), Vec3.create(size[0] / (grid[0] - 1), size[1] / (grid[1] - 1), size[2] / (grid[2] - 1)));

    const transform = Mat4.mul(Mat4.zero(), translation, scale);
    Mesh.transformImmediate(surface, transform);
    Mesh.computeNormalsImmediate(surface);
    return { surface, field };
}