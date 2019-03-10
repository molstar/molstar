/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Vec3, Mat4, Mat3 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';
import { Mesh } from './mesh';
import { getNormalMatrix } from '../../util';
import { Primitive } from '../../primitive/primitive';
import { Cage } from 'mol-geo/primitive/cage';
import { addSphere } from './builder/sphere';
import { addCylinder } from './builder/cylinder';

const tmpV = Vec3.zero()
const tmpMat3 = Mat3.zero()
const tmpVecA = Vec3.zero()
const tmpVecB = Vec3.zero()

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
        }
    }

    export function addPrimitive(state: State, t: Mat4, primitive: Primitive) {
        const { vertices: va, normals: na, indices: ia } = primitive
        const { vertices, normals, indices, groups, currentGroup } = state
        const offset = vertices.elementCount
        const n = getNormalMatrix(tmpMat3, t)
        for (let i = 0, il = va.length; i < il; i += 3) {
            // position
            Vec3.transformMat4(tmpV, Vec3.fromArray(tmpV, va, i), t)
            ChunkedArray.add3(vertices, tmpV[0], tmpV[1], tmpV[2]);
            // normal
            Vec3.transformMat3(tmpV, Vec3.fromArray(tmpV, na, i), n)
            ChunkedArray.add3(normals, tmpV[0], tmpV[1], tmpV[2]);
            // group
            ChunkedArray.add(groups, currentGroup);
        }
        for (let i = 0, il = ia.length; i < il; i += 3) {
            ChunkedArray.add3(indices, ia[i] + offset, ia[i + 1] + offset, ia[i + 2] + offset);
        }
    }

    export function addCage(state: State, t: Mat4, cage: Cage, radius: number, detail: number) {
        const { vertices: va, edges: ea } = cage
        const cylinderProps = { radiusTop: radius, radiusBottom: radius }
        for (let i = 0, il = ea.length; i < il; i += 2) {
            Vec3.fromArray(tmpVecA, va, ea[i] * 3)
            Vec3.fromArray(tmpVecB, va, ea[i + 1] * 3)
            Vec3.transformMat4(tmpVecA, tmpVecA, t)
            Vec3.transformMat4(tmpVecB, tmpVecB, t)
            addSphere(state, tmpVecA, radius, detail)
            addSphere(state, tmpVecB, radius, detail)
            addCylinder(state, tmpVecA, tmpVecB, 1, cylinderProps)
        }
    }

    export function getMesh (state: State): Mesh {
        const { vertices, normals, indices, groups, mesh } = state
        const vb = ChunkedArray.compact(vertices, true) as Float32Array
        const ib = ChunkedArray.compact(indices, true) as Uint32Array
        const nb = ChunkedArray.compact(normals, true) as Float32Array
        const gb = ChunkedArray.compact(groups, true) as Float32Array
        return {
            kind: 'mesh',
            vertexCount: state.vertices.elementCount,
            triangleCount: state.indices.elementCount,
            vertexBuffer: mesh ? ValueCell.update(mesh.vertexBuffer, vb) : ValueCell.create(vb),
            indexBuffer: mesh ? ValueCell.update(mesh.indexBuffer, ib) : ValueCell.create(ib),
            normalBuffer: mesh ? ValueCell.update(mesh.normalBuffer, nb) : ValueCell.create(nb),
            groupBuffer: mesh ? ValueCell.update(mesh.groupBuffer, gb) : ValueCell.create(gb),
            normalsComputed: true,
        }
    }
}