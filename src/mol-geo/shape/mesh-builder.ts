/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';

import Box, { BoxProps } from '../primitive/box';
import { Mesh } from './mesh';

type ElementId = { id: number }

export interface MeshBuilder {
    add(t: Mat4, _vertices: Float32Array, _normals: Float32Array, _indices?: Uint32Array): number
    addBox(t: Mat4, props?: BoxProps & ElementId): number
    getMesh(): Mesh
}

const tmpV = Vec3.zero()

export namespace MeshBuilder {
    export function create(initialCount = 2048, chunkSize = 1024): MeshBuilder {
        const vertices = ChunkedArray.create(Float32Array, 3, chunkSize, initialCount);
        const normals = ChunkedArray.create(Float32Array, 3, chunkSize, initialCount);
        const indices = ChunkedArray.create(Uint32Array, 3, chunkSize * 3, initialCount * 3);

        // const offsets = ChunkedArray.create<number>(n => new Uint32Array(n), 1, 1000);
        // const elementIds = ChunkedArray.create(Uint32Array, 1, chunkSize, initialCount);

        ChunkedArray.compact(indices, true)

        const add = (t: Mat4, _vertices: Float32Array, _normals: Float32Array, _indices?: Uint32Array) => {
            const offset = vertices.elementCount * vertices.elementSize
            for (let i = 0, il = _vertices.length; i < il; i += 3) {
                Vec3.fromArray(tmpV, _vertices, i)
                Vec3.transformMat4(tmpV, tmpV, t)
                // Vec3.transformDirection(tmpV, tmpV, n)  // TODO
                ChunkedArray.add3(vertices, tmpV[0], tmpV[1], tmpV[2]);
            }
            // ChunkedArray.add(vertices, _vertices[i])
            return offset
        }

        return {
            add,
            addBox: (t: Mat4, props?: BoxProps & ElementId) => {
                const box = Box(props)
                return add(t, box.vertices, box.normals, box.indices)
            },
            getMesh: () => {
                return {
                    vertexCount: vertices.elementCount,
                    triangleCount: indices.elementCount,
                    vertexBuffer: ValueCell.create(ChunkedArray.compact(vertices, true) as Float32Array),
                    indexBuffer: ValueCell.create(ChunkedArray.compact(indices, true) as Uint32Array),
                    normalBuffer: ValueCell.create(ChunkedArray.compact(normals, true) as Float32Array),
                    normalsComputed: true,
                }
            }
        }
    }
}