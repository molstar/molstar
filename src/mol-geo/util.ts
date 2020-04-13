/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4, Mat3 } from '../mol-math/linear-algebra';
import { NumberArray } from '../mol-util/type-helpers';
import { arrayMax } from '../mol-util/array';

export function normalizeVec3Array<T extends NumberArray> (a: T, count: number) {
    for (let i = 0, il = count * 3; i < il; i += 3) {
        const x = a[i];
        const y = a[i + 1];
        const z = a[i + 2];
        const s = 1 / Math.sqrt(x * x + y * y + z * z);
        a[i] = x * s;
        a[i + 1] = y * s;
        a[i + 2] = z * s;
    }
    return a;
}

const tmpV3 = Vec3.zero();

export function transformPositionArray (t: Mat4, array: NumberArray, offset: number, count: number) {
    for (let i = 0, il = count * 3; i < il; i += 3) {
        Vec3.fromArray(tmpV3, array, offset + i);
        Vec3.transformMat4(tmpV3, tmpV3, t);
        Vec3.toArray(tmpV3, array, offset + i);
    }
}

export function transformDirectionArray (n: Mat3, array: NumberArray, offset: number, count: number) {
    for (let i = 0, il = count * 3; i < il; i += 3) {
        Vec3.fromArray(tmpV3, array, offset + i);
        Vec3.transformMat3(tmpV3, tmpV3, n);
        Vec3.toArray(tmpV3, array, offset + i);
    }
}

/** iterate over the entire array and apply the radius to each vertex */
export function appplyRadius(vertices: NumberArray, radius: number) {
    for (let i = 0, il = vertices.length; i < il; i += 3) {
        Vec3.fromArray(tmpV3, vertices, i);
        Vec3.normalize(tmpV3, tmpV3);
        Vec3.scale(tmpV3, tmpV3, radius);
        Vec3.toArray(tmpV3, vertices, i);
    }
}

const a = Vec3();
const b = Vec3();
const c = Vec3();
const cb = Vec3();
const ab = Vec3();

/**
 * indexed vertex normals weighted by triangle areas
 *      http://www.iquilezles.org/www/articles/normals/normals.htm
 * - normals array must contain only zeros
 */
export function computeIndexedVertexNormals<T extends NumberArray> (vertices: NumberArray, indices: NumberArray, normals: T, vertexCount: number, triangleCount: number) {

    for (let i = 0, il = triangleCount * 3; i < il; i += 3) {
        const ai = indices[i] * 3;
        const bi = indices[i + 1] * 3;
        const ci = indices[i + 2] * 3;

        Vec3.fromArray(a, vertices, ai);
        Vec3.fromArray(b, vertices, bi);
        Vec3.fromArray(c, vertices, ci);

        Vec3.sub(cb, c, b);
        Vec3.sub(ab, a, b);
        Vec3.cross(cb, cb, ab);

        normals[ai] += cb[0];
        normals[ai + 1] += cb[1];
        normals[ai + 2] += cb[2];

        normals[bi] += cb[0];
        normals[bi + 1] += cb[1];
        normals[bi + 2] += cb[2];

        normals[ci] += cb[0];
        normals[ci + 1] += cb[1];
        normals[ci + 2] += cb[2];
    }

    return normalizeVec3Array(normals, vertexCount);
}

/**
 * vertex normals for unindexed triangle soup
 * - normals array must contain only zeros
 */
export function computeVertexNormals<T extends NumberArray> (vertices: NumberArray, normals: T, vertexCount: number) {
    for (let i = 0, il = vertexCount * 3; i < il; i += 9) {
        Vec3.fromArray(a, vertices, i);
        Vec3.fromArray(b, vertices, i + 3);
        Vec3.fromArray(c, vertices, i + 6);

        Vec3.sub(cb, c, b);
        Vec3.sub(ab, a, b);
        Vec3.cross(cb, cb, ab);

        normals[i] = cb[0];
        normals[i + 1] = cb[1];
        normals[i + 2] = cb[2];

        normals[i + 3] = cb[0];
        normals[i + 4] = cb[1];
        normals[i + 5] = cb[2];

        normals[i + 6] = cb[0];
        normals[i + 7] = cb[1];
        normals[i + 8] = cb[2];
    }

    return normalizeVec3Array(normals, vertexCount);
}

/**
 * Maps groups to data, range for group i is offsets[i] to offsets[i + 1]
 */
export type GroupMapping = {
    /** data indices */
    readonly indices: ArrayLike<number>
    /** range for group i is offsets[i] to offsets[i + 1] */
    readonly offsets: ArrayLike<number>
}

/**
 * The `step` parameter allows to skip over repeated values in `groups`
 */
export function createGroupMapping(groups: ArrayLike<number>, dataCount: number, step = 1): GroupMapping {
    const maxId = arrayMax(groups);

    const offsets = new Int32Array(maxId + 2);
    const bucketFill = new Int32Array(dataCount);
    const bucketSizes = new Int32Array(dataCount);

    for (let i = 0, il = dataCount * step; i < il; i += step) ++bucketSizes[groups[i]];

    let offset = 0;
    for (let i = 0; i < dataCount; i++) {
        offsets[i] = offset;
        offset += bucketSizes[i];
    }
    offsets[dataCount] = offset;

    const indices = new Int32Array(offset);
    for (let i = 0, il = dataCount * step; i < il; i += step) {
        const g = groups[i];
        const og = offsets[g] + bucketFill[g];
        indices[og] = i;
        ++bucketFill[g];
    }

    return { indices, offsets };
}
