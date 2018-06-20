/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Vec3, Mat4, Mat3 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';

import Box, { BoxProps } from '../primitive/box';
import Cylinder, { CylinderProps } from '../primitive/cylinder';
import Icosahedron, { IcosahedronProps } from '../primitive/icosahedron';
import { Mesh } from './mesh';
import { getNormalMatrix } from '../util';

type Primitive = {
    vertices: Float32Array
    normals: Float32Array
    indices: Uint32Array
}

export interface MeshBuilder {
    add(t: Mat4, _vertices: Float32Array, _normals: Float32Array, _indices?: Uint32Array): number
    addBox(t: Mat4, props?: BoxProps): number
    addCylinder(t: Mat4, props?: CylinderProps): number
    addIcosahedron(t: Mat4, props?: IcosahedronProps): number
    setId(id: number): void
    getMesh(): Mesh
}

const tmpV = Vec3.zero()
const tmpMat3 = Mat3.zero()

// TODO cache primitives based on props

export namespace MeshBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, mesh?: Mesh): MeshBuilder {
        const vertices = ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.vertexBuffer.ref.value : initialCount);
        const normals = ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.normalBuffer.ref.value : initialCount);
        const indices = ChunkedArray.create(Uint32Array, 3, chunkSize * 3, mesh ? mesh.indexBuffer.ref.value : initialCount * 3);

        const ids = ChunkedArray.create(Float32Array, 1, chunkSize, mesh ? mesh.idBuffer.ref.value : initialCount);
        const offsets = ChunkedArray.create(Uint32Array, 1, chunkSize, mesh ? mesh.offsetBuffer.ref.value : initialCount);

        let currentId = -1

        const cylinderMap = new Map<string, Primitive>()
        const icosahedronMap = new Map<string, Primitive>()

        const add = (t: Mat4, _vertices: Float32Array, _normals: Float32Array, _indices: Uint32Array) => {
            const { elementCount, elementSize } = vertices
            const n = getNormalMatrix(tmpMat3, t)
            for (let i = 0, il = _vertices.length; i < il; i += 3) {
                // position
                Vec3.fromArray(tmpV, _vertices, i)
                Vec3.transformMat4(tmpV, tmpV, t)
                ChunkedArray.add3(vertices, tmpV[0], tmpV[1], tmpV[2]);
                // normal
                Vec3.fromArray(tmpV, _normals, i)
                Vec3.transformMat3(tmpV, tmpV, n)
                ChunkedArray.add3(normals, tmpV[0], tmpV[1], tmpV[2]);

                ChunkedArray.add(ids, currentId);
            }
            for (let i = 0, il = _indices.length; i < il; i += 3) {
                ChunkedArray.add3(indices, _indices[i] + elementCount, _indices[i + 1] + elementCount, _indices[i + 2] + elementCount);
            }
            return elementCount * elementSize
        }

        return {
            add,
            addBox: (t: Mat4, props?: BoxProps) => {
                const box = Box(props)
                return add(t, box.vertices, box.normals, box.indices)
            },
            addCylinder: (t: Mat4, props?: CylinderProps) => {
                const key = JSON.stringify(props)
                let cylinder = cylinderMap.get(key)
                if (cylinder === undefined) {
                    cylinder = Cylinder(props)
                    cylinderMap.set(key, cylinder)
                }
                return add(t, cylinder.vertices, cylinder.normals, cylinder.indices)
            },
            addIcosahedron: (t: Mat4, props: IcosahedronProps) => {
                const key = JSON.stringify(props)
                let icosahedron = icosahedronMap.get(key)
                if (icosahedron === undefined) {
                    icosahedron = Icosahedron(props)
                    icosahedronMap.set(key, icosahedron)
                }
                return add(t, icosahedron.vertices, icosahedron.normals, icosahedron.indices)
            },
            setId: (id: number) => {
                if (currentId !== id) {
                    currentId = id
                    ChunkedArray.add(offsets, vertices.elementCount)
                }
            },
            getMesh: () => {
                ChunkedArray.add(offsets, vertices.elementCount)
                const vb = ChunkedArray.compact(vertices, true) as Float32Array
                const ib = ChunkedArray.compact(indices, true) as Uint32Array
                const nb = ChunkedArray.compact(normals, true) as Float32Array
                const idb = ChunkedArray.compact(ids, true) as Float32Array
                const ob = ChunkedArray.compact(offsets, true) as Uint32Array
                return {
                    vertexCount: vertices.elementCount,
                    triangleCount: indices.elementCount,
                    offsetCount: offsets.elementCount,
                    vertexBuffer: mesh ? ValueCell.update(mesh.vertexBuffer, vb) : ValueCell.create(vb),
                    indexBuffer: mesh ? ValueCell.update(mesh.indexBuffer, ib) : ValueCell.create(ib),
                    normalBuffer: mesh ? ValueCell.update(mesh.normalBuffer, nb) : ValueCell.create(nb),
                    idBuffer: mesh ? ValueCell.update(mesh.idBuffer, idb) : ValueCell.create(idb),
                    offsetBuffer: mesh ? ValueCell.update(mesh.offsetBuffer, ob) : ValueCell.create(ob),
                    normalsComputed: true,
                    offsetsComputed: true,
                }
            }
        }
    }
}