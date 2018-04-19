/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';

import Box, { BoxProps } from '../primitive/box';
import Cylinder, { CylinderProps } from '../primitive/cylinder';
import Icosahedron, { IcosahedronProps } from '../primitive/icosahedron';
import { Mesh } from './mesh';

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

// TODO cache primitives based on props

export namespace MeshBuilder {
    export function create(initialCount = 2048, chunkSize = 1024): MeshBuilder {
        const vertices = ChunkedArray.create(Float32Array, 3, chunkSize, initialCount);
        const normals = ChunkedArray.create(Float32Array, 3, chunkSize, initialCount);
        const indices = ChunkedArray.create(Uint32Array, 3, chunkSize * 3, initialCount * 3);

        const ids = ChunkedArray.create(Float32Array, 1, chunkSize, initialCount);
        const offsets = ChunkedArray.create(Uint32Array, 1, chunkSize, initialCount);

        let currentId = -1

        const icosahedronMap = new Map<string, Primitive>()

        const add = (t: Mat4, _vertices: Float32Array, _normals: Float32Array, _indices: Uint32Array) => {
            const { elementCount, elementSize } = vertices
            for (let i = 0, il = _vertices.length; i < il; i += 3) {
                // position
                Vec3.fromArray(tmpV, _vertices, i)
                Vec3.transformMat4(tmpV, tmpV, t)
                ChunkedArray.add3(vertices, tmpV[0], tmpV[1], tmpV[2]);
                // normal
                Vec3.fromArray(tmpV, _normals, i)
                // Vec3.transformDirection(tmpV, tmpV, n)  // TODO
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
                const cylinder = Cylinder(props)
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
                const mesh = {
                    vertexCount: vertices.elementCount,
                    triangleCount: indices.elementCount,
                    offsetCount: offsets.elementCount,
                    vertexBuffer: ValueCell.create(ChunkedArray.compact(vertices, true) as Float32Array),
                    indexBuffer: ValueCell.create(ChunkedArray.compact(indices, true) as Uint32Array),
                    normalBuffer: ValueCell.create(ChunkedArray.compact(normals, true) as Float32Array),
                    idBuffer: ValueCell.create(ChunkedArray.compact(ids, true) as Float32Array),
                    offsetBuffer: ValueCell.create(ChunkedArray.compact(offsets, true) as Uint32Array),
                    normalsComputed: true,
                }
                return mesh
            }
        }
    }
}