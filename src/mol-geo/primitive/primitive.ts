/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4, Mat3 } from '../../mol-math/linear-algebra';
import { NumberArray } from '../../mol-util/type-helpers';

export interface Primitive {
    vertices: ArrayLike<number>
    normals: ArrayLike<number>
    indices: ArrayLike<number>
}

const a = Vec3(), b = Vec3(), c = Vec3();

/** Create primitive with face normals from vertices and indices */
export function createPrimitive(vertices: ArrayLike<number>, indices: ArrayLike<number>): Primitive {
    const count = indices.length;
    const builder = PrimitiveBuilder(count / 3);

    for (let i = 0; i < count; i += 3) {
        Vec3.fromArray(a, vertices, indices[i] * 3);
        Vec3.fromArray(b, vertices, indices[i + 1] * 3);
        Vec3.fromArray(c, vertices, indices[i + 2] * 3);
        builder.add(a, b, c);
    }
    return builder.getPrimitive();
}

export function copyPrimitive(primitive: Primitive): Primitive {
    return {
        vertices: new Float32Array(primitive.vertices),
        normals: new Float32Array(primitive.normals),
        indices: new Uint32Array(primitive.indices)
    };
}

export interface PrimitiveBuilder {
    add(a: Vec3, b: Vec3, c: Vec3): void
    /** Shared vertices and normals, must be flat */
    addQuad(a: Vec3, b: Vec3, c: Vec3, d: Vec3): void
    getPrimitive(): Primitive
}

const vn = Vec3();

/** Builder to create primitive with face normals */
export function PrimitiveBuilder(triangleCount: number, vertexCount?: number): PrimitiveBuilder {
    if (vertexCount === undefined) vertexCount = triangleCount * 3;

    const vertices = new Float32Array(vertexCount * 3);
    const normals = new Float32Array(vertexCount * 3);
    const indices = new Uint32Array(triangleCount * 3);

    let vOffset = 0;
    let iOffset = 0;

    return {
        add: (a: Vec3, b: Vec3, c: Vec3) => {
            Vec3.toArray(a, vertices, vOffset);
            Vec3.toArray(b, vertices, vOffset + 3);
            Vec3.toArray(c, vertices, vOffset + 6);
            Vec3.triangleNormal(vn, a, b, c);
            for (let j = 0; j < 3; ++j) {
                Vec3.toArray(vn, normals, vOffset + 3 * j);
                indices[iOffset + j] = vOffset / 3 + j;
            }
            vOffset += 9;
            iOffset += 3;
        },
        addQuad: (a: Vec3, b: Vec3, c: Vec3, d: Vec3) => {
            Vec3.toArray(a, vertices, vOffset);
            Vec3.toArray(b, vertices, vOffset + 3);
            Vec3.toArray(c, vertices, vOffset + 6);
            Vec3.toArray(d, vertices, vOffset + 9);
            Vec3.triangleNormal(vn, a, b, c);
            for (let j = 0; j < 4; ++j) {
                Vec3.toArray(vn, normals, vOffset + 3 * j);
            }
            const vOffset3 = vOffset / 3;
            // a, b, c
            indices[iOffset] = vOffset3;
            indices[iOffset + 1] = vOffset3 + 1;
            indices[iOffset + 2] = vOffset3 + 2;
            // a, b, c
            indices[iOffset + 3] = vOffset3 + 2;
            indices[iOffset + 4] = vOffset3 + 3;
            indices[iOffset + 5] = vOffset3;
            vOffset += 12;
            iOffset += 6;
        },
        getPrimitive: () => ({ vertices, normals, indices })
    };
}

const tmpV = Vec3();
const tmpMat3 = Mat3();

/** Transform primitive in-place */
export function transformPrimitive(primitive: Primitive, t: Mat4) {
    const { vertices, normals } = primitive;
    const n = Mat3.directionTransform(tmpMat3, t);
    for (let i = 0, il = vertices.length; i < il; i += 3) {
        // position
        Vec3.transformMat4(tmpV, Vec3.fromArray(tmpV, vertices, i), t);
        Vec3.toArray(tmpV, vertices as NumberArray, i);
        // normal
        Vec3.transformMat3(tmpV, Vec3.fromArray(tmpV, normals, i), n);
        Vec3.toArray(tmpV, normals as NumberArray, i);
    }
    return primitive;
}