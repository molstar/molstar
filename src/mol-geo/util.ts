/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4, Mat3 } from 'mol-math/linear-algebra'

export function normalizeVec3Array<T extends Helpers.NumberArray> (a: T) {
    const n = a.length
    for (let i = 0; i < n; i += 3) {
        const x = a[ i ]
        const y = a[ i + 1 ]
        const z = a[ i + 2 ]
        const s = 1 / Math.sqrt(x * x + y * y + z * z)
        a[ i ] = x * s
        a[ i + 1 ] = y * s
        a[ i + 2 ] = z * s
    }
}

const tmpV = Vec3.zero()

export function transformPositionArray (t: Mat4, array: Helpers.NumberArray, offset: number, count: number) {
    for (let i = 0, il = count * 3; i < il; i += 3) {
        Vec3.fromArray(tmpV, array, offset + i)
        Vec3.transformMat4(tmpV, tmpV, t)
        Vec3.toArray(tmpV, array, offset + i)
    }
}

export function transformDirectionArray (t: Mat3, array: Helpers.NumberArray, offset: number, count: number) {
    // TODO
}

export function setArrayZero(array: Helpers.NumberArray) {
    const n = array.length
    for (let i = 0; i < n; ++i) array[i] = 0
}

/** iterate over the entire buffer and apply the radius to each vertex */
export function appplyRadius(vertices: Helpers.NumberArray, radius: number) {
    const v = Vec3.zero()
    const n = vertices.length
    for (let i = 0; i < n; i += 3) {
        Vec3.fromArray(v, vertices, i)
        Vec3.normalize(v, v)
        Vec3.scale(v, v, radius)
        Vec3.toArray(v, vertices, i)
    }
}

/**
 * indexed vertex normals weighted by triangle areas http://www.iquilezles.org/www/articles/normals/normals.htm
 * normal array must contain only zeros
 */
export function computeIndexedVertexNormals<T extends Helpers.NumberArray> (vertices: Helpers.NumberArray, indices: Helpers.NumberArray, normals: T) {
    const a = Vec3.zero()
    const b = Vec3.zero()
    const c = Vec3.zero()
    const cb = Vec3.zero()
    const ab = Vec3.zero()

    for (let i = 0, il = indices.length; i < il; i += 3) {
        const ai = indices[ i ] * 3
        const bi = indices[ i + 1 ] * 3
        const ci = indices[ i + 2 ] * 3

        Vec3.fromArray(a, vertices, ai)
        Vec3.fromArray(b, vertices, bi)
        Vec3.fromArray(c, vertices, ci)

        Vec3.sub(cb, c, b)
        Vec3.sub(ab, a, b)
        Vec3.cross(cb, cb, ab)

        normals[ ai ] += cb[ 0 ]
        normals[ ai + 1 ] += cb[ 1 ]
        normals[ ai + 2 ] += cb[ 2 ]

        normals[ bi ] += cb[ 0 ]
        normals[ bi + 1 ] += cb[ 1 ]
        normals[ bi + 2 ] += cb[ 2 ]

        normals[ ci ] += cb[ 0 ]
        normals[ ci + 1 ] += cb[ 1 ]
        normals[ ci + 2 ] += cb[ 2 ]
    }

    normalizeVec3Array(normals)
    return normals
}

/** vertex normals for unindexed triangle soup, normal array must contain only zeros */
export function computeVertexNormals<T extends Helpers.NumberArray> (vertices: Helpers.NumberArray, normals: T) {
    setArrayZero(normals)

    const a = Vec3.zero()
    const b = Vec3.zero()
    const c = Vec3.zero()
    const cb = Vec3.zero()
    const ab = Vec3.zero()

     for (let i = 0, il = vertices.length; i < il; i += 9) {
        Vec3.fromArray(a, vertices, i)
        Vec3.fromArray(b, vertices, i + 3)
        Vec3.fromArray(c, vertices, i + 6)

        Vec3.sub(cb, c, b)
        Vec3.sub(ab, a, b)
        Vec3.cross(cb, cb, ab)

        normals[ i ] = cb[ 0 ]
        normals[ i + 1 ] = cb[ 1 ]
        normals[ i + 2 ] = cb[ 2 ]

        normals[ i + 3 ] = cb[ 0 ]
        normals[ i + 4 ] = cb[ 1 ]
        normals[ i + 5 ] = cb[ 2 ]

        normals[ i + 6 ] = cb[ 0 ]
        normals[ i + 7 ] = cb[ 1 ]
        normals[ i + 8 ] = cb[ 2 ]
    }

    normalizeVec3Array(normals)
    return normals
}