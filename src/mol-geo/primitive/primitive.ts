/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra';

export interface Primitive {
    vertices: ArrayLike<number>
    normals: ArrayLike<number>
    indices: ArrayLike<number>
}

const tri = [ Vec3.zero(), Vec3.zero(), Vec3.zero() ]
const n = Vec3.zero()

/** Create primitive with face normals from vertices and indices */
export function createPrimitive(_vertices: ArrayLike<number>, _indices: ArrayLike<number>): Primitive {
    const count = _indices.length
    const vertices = new Float32Array(count * 3)
    const normals = new Float32Array(count * 3)
    const indices = new Uint32Array(count)

    for (let i = 0; i < count; i += 3) {
        for (let j = 0; j < 3; ++j) {
            Vec3.fromArray(tri[j], _vertices, _indices[i + j] * 3)
            Vec3.toArray(tri[j], vertices, i * 3 + j * 3)
        }

        Vec3.triangleNormal(n, tri[0], tri[1], tri[2])

        for (let j = 0; j < 3; ++j) {
            Vec3.toArray(n, normals, i * 3 + j * 3)
            indices[i + j] = i + j
        }
    }

    return { vertices, normals, indices }
}