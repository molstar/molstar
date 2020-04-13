/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Mat4, Mat3 } from '../../../mol-math/linear-algebra';
import { ChunkedArray } from '../../../mol-data/util';
import { Mesh } from './mesh';
import { Primitive } from '../../primitive/primitive';
import { Cage } from '../../../mol-geo/primitive/cage';
import { addSphere } from './builder/sphere';
import { addCylinder } from './builder/cylinder';

const tmpV = Vec3.zero();
const tmpMat3 = Mat3.zero();
const tmpVecA = Vec3.zero();
const tmpVecB = Vec3.zero();
const tmpVecC = Vec3.zero();
const tmpVecD = Vec3.zero();

export namespace MeshBuilder {
    export interface State {
        currentGroup: number
        readonly vertices: ChunkedArray<number, 3>
        readonly normals: ChunkedArray<number, 3>
        readonly indices: ChunkedArray<number, 3>
        readonly groups: ChunkedArray<number, 1>
        readonly mesh?: Mesh
    }

    export function createState(initialCount = 2048, chunkSize = 1024, mesh?: Mesh): State {
        return {
            currentGroup: -1,
            vertices: ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.vertexBuffer.ref.value : initialCount),
            normals: ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.normalBuffer.ref.value : initialCount),
            indices: ChunkedArray.create(Uint32Array, 3, chunkSize * 3, mesh ? mesh.indexBuffer.ref.value : initialCount * 3),
            groups: ChunkedArray.create(Float32Array, 1, chunkSize, mesh ? mesh.groupBuffer.ref.value : initialCount),
            mesh
        };
    }

    export function addTriangle(state: State, a: Vec3, b: Vec3, c: Vec3) {
        const { vertices, normals, indices, groups, currentGroup } = state;
        const offset = vertices.elementCount;

        // positions
        ChunkedArray.add3(vertices, a[0], a[1], a[2]);
        ChunkedArray.add3(vertices, b[0], b[1], b[2]);
        ChunkedArray.add3(vertices, c[0], c[1], c[2]);

        Vec3.triangleNormal(tmpV, a, b, c);
        for (let i = 0; i < 3; ++i) {
            ChunkedArray.add3(normals, tmpV[0], tmpV[1], tmpV[2]);  // normal
            ChunkedArray.add(groups, currentGroup);  // group
        }
        ChunkedArray.add3(indices, offset, offset + 1, offset + 2);
    }

    export function addTriangleStrip(state: State, vertices: ArrayLike<number>, indices: ArrayLike<number>) {
        Vec3.fromArray(tmpVecC, vertices, indices[0] * 3);
        Vec3.fromArray(tmpVecD, vertices, indices[1] * 3);
        for (let i = 2, il = indices.length; i < il; i += 2) {
            Vec3.copy(tmpVecA, tmpVecC);
            Vec3.copy(tmpVecB, tmpVecD);
            Vec3.fromArray(tmpVecC, vertices, indices[i] * 3);
            Vec3.fromArray(tmpVecD, vertices, indices[i + 1] * 3);
            addTriangle(state, tmpVecA, tmpVecB, tmpVecC);
            addTriangle(state, tmpVecB, tmpVecD, tmpVecC);
        }
    }

    export function addTriangleFan(state: State, vertices: ArrayLike<number>, indices: ArrayLike<number>) {
        Vec3.fromArray(tmpVecA, vertices, indices[0] * 3);
        for (let i = 2, il = indices.length; i < il; ++i) {
            Vec3.fromArray(tmpVecB, vertices, indices[i - 1] * 3);
            Vec3.fromArray(tmpVecC, vertices, indices[i] * 3);
            addTriangle(state, tmpVecA, tmpVecC, tmpVecB);
        }
    }

    export function addPrimitive(state: State, t: Mat4, primitive: Primitive) {
        const { vertices: va, normals: na, indices: ia } = primitive;
        const { vertices, normals, indices, groups, currentGroup } = state;
        const offset = vertices.elementCount;
        const n = Mat3.directionTransform(tmpMat3, t);
        for (let i = 0, il = va.length; i < il; i += 3) {
            // position
            Vec3.transformMat4(tmpV, Vec3.fromArray(tmpV, va, i), t);
            ChunkedArray.add3(vertices, tmpV[0], tmpV[1], tmpV[2]);
            // normal
            Vec3.transformMat3(tmpV, Vec3.fromArray(tmpV, na, i), n);
            ChunkedArray.add3(normals, tmpV[0], tmpV[1], tmpV[2]);
            // group
            ChunkedArray.add(groups, currentGroup);
        }
        for (let i = 0, il = ia.length; i < il; i += 3) {
            ChunkedArray.add3(indices, ia[i] + offset, ia[i + 1] + offset, ia[i + 2] + offset);
        }
    }

    /** Flips triangle normals and winding order */
    export function addPrimitiveFlipped(state: State, t: Mat4, primitive: Primitive) {
        const { vertices: va, normals: na, indices: ia } = primitive;
        const { vertices, normals, indices, groups, currentGroup } = state;
        const offset = vertices.elementCount;
        const n = Mat3.directionTransform(tmpMat3, t);
        for (let i = 0, il = va.length; i < il; i += 3) {
            // position
            Vec3.transformMat4(tmpV, Vec3.fromArray(tmpV, va, i), t);
            ChunkedArray.add3(vertices, tmpV[0], tmpV[1], tmpV[2]);
            // normal
            Vec3.transformMat3(tmpV, Vec3.fromArray(tmpV, na, i), n);
            ChunkedArray.add3(normals, -tmpV[0], -tmpV[1], -tmpV[2]);
            // group
            ChunkedArray.add(groups, currentGroup);
        }
        for (let i = 0, il = ia.length; i < il; i += 3) {
            ChunkedArray.add3(indices, ia[i + 2] + offset, ia[i + 1] + offset, ia[i] + offset);
        }
    }

    export function addCage(state: State, t: Mat4, cage: Cage, radius: number, detail: number, radialSegments: number) {
        const { vertices: va, edges: ea } = cage;
        const cylinderProps = { radiusTop: radius, radiusBottom: radius, radialSegments };
        for (let i = 0, il = ea.length; i < il; i += 2) {
            Vec3.fromArray(tmpVecA, va, ea[i] * 3);
            Vec3.fromArray(tmpVecB, va, ea[i + 1] * 3);
            Vec3.transformMat4(tmpVecA, tmpVecA, t);
            Vec3.transformMat4(tmpVecB, tmpVecB, t);
            addSphere(state, tmpVecA, radius, detail);
            addSphere(state, tmpVecB, radius, detail);
            addCylinder(state, tmpVecA, tmpVecB, 1, cylinderProps);
        }
    }

    export function getMesh (state: State): Mesh {
        const { vertices, normals, indices, groups, mesh } = state;
        const vb = ChunkedArray.compact(vertices, true) as Float32Array;
        const ib = ChunkedArray.compact(indices, true) as Uint32Array;
        const nb = ChunkedArray.compact(normals, true) as Float32Array;
        const gb = ChunkedArray.compact(groups, true) as Float32Array;
        return Mesh.create(vb, ib, nb, gb, state.vertices.elementCount, state.indices.elementCount, mesh);
    }
}