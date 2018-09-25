/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Fred Ludlow <fred.ludlow@gmail.com>
 */

import { ChunkedArray } from '../../../mol-data/util';
import { ValueCell } from 'mol-util';
import { Mesh } from '../../geometry/mesh/mesh';
import { AllowedContours } from './tables';
import { LinesBuilder } from '../../geometry/lines/lines-builder';
import { Lines } from '../../geometry/lines/lines';

 export interface MarchinCubesBuilder<T> {
    addVertex(x: number, y: number, z: number): number
    addGroup(group: number): void
    addTriangle(vertList: number[], a: number, b: number, c: number, edgeFilter: number): void
    get(): T
}

export function MarchinCubesMeshBuilder(vertexChunkSize: number, mesh?: Mesh): MarchinCubesBuilder<Mesh> {
    const triangleChunkSize = Math.min(1 << 16, vertexChunkSize * 4)

    const vertices = ChunkedArray.create(Float32Array, 3, vertexChunkSize, mesh && mesh.vertexBuffer.ref.value);
    const groups = ChunkedArray.create(Float32Array, 1, vertexChunkSize, mesh && mesh.groupBuffer.ref.value);
    const indices = ChunkedArray.create(Uint32Array, 3, triangleChunkSize, mesh && mesh.indexBuffer.ref.value);

    let vertexCount = 0
    let triangleCount = 0

    return {
        addVertex: (x: number, y: number, z: number) => {
            ++vertexCount
            return ChunkedArray.add3(vertices, x, y, z );
        },
        addGroup: (group: number) => {
            ChunkedArray.add(groups, group);
        },
        addTriangle: (vertList: number[], a: number, b: number, c: number) => {
            ++triangleCount
            ChunkedArray.add3(indices, vertList[a], vertList[b], vertList[c]);
        },
        get: () => {
            const vb = ChunkedArray.compact(vertices, true) as Float32Array;
            const ib = ChunkedArray.compact(indices, true) as Uint32Array;
            const gb = ChunkedArray.compact(groups, true) as Float32Array;

            return {
                kind: 'mesh',
                vertexCount,
                triangleCount,
                vertexBuffer: mesh ? ValueCell.update(mesh.vertexBuffer, vb) : ValueCell.create(vb),
                groupBuffer: mesh ? ValueCell.update(mesh.groupBuffer, gb) : ValueCell.create(gb),
                indexBuffer: mesh ? ValueCell.update(mesh.indexBuffer, ib) : ValueCell.create(ib),
                normalBuffer: mesh ? mesh.normalBuffer : ValueCell.create(new Float32Array(0)),
                normalsComputed: false
            }
        }
    }
}

export function MarchinCubesLinesBuilder(vertexChunkSize: number, lines?: Lines): MarchinCubesBuilder<Lines> {
    const vertices = ChunkedArray.create(Float32Array, 3, vertexChunkSize);
    const groups = ChunkedArray.create(Float32Array, 1, vertexChunkSize);
    const indices = ChunkedArray.create(Float32Array, 2, vertexChunkSize);

    let linesCount = 0

    return {
        addVertex: (x: number, y: number, z: number) => {
            return ChunkedArray.add3(vertices, x, y, z);
        },
        addGroup: (group: number) => {
            ChunkedArray.add(groups, group);
        },
        addTriangle: (vertList: number[], a: number, b: number, c: number, edgeFilter: number) => {
            if (AllowedContours[a][b] & edgeFilter) {
                ++linesCount
                ChunkedArray.add2(indices, vertList[a], vertList[b])
            }
            if (AllowedContours[b][c] & edgeFilter) {
                ++linesCount
                ChunkedArray.add2(indices, vertList[b], vertList[c])
            }
            if (AllowedContours[a][c] & edgeFilter) {
                ++linesCount
                ChunkedArray.add2(indices, vertList[a], vertList[c])
            }
        },
        get: () => {
            const vb = ChunkedArray.compact(vertices, true) as Float32Array;
            const ib = ChunkedArray.compact(indices, true) as Uint32Array;
            const gb = ChunkedArray.compact(groups, true) as Float32Array;

            const builder = LinesBuilder.create(linesCount, linesCount / 10, lines)

            for (let i = 0; i < linesCount; ++i) {
                const la = ib[i * 2], lb = ib[i * 2 + 1]
                builder.add(
                    vb[la * 3], vb[la * 3 + 1], vb[la * 3 + 2],
                    vb[lb * 3], vb[lb * 3 + 1], vb[lb * 3 + 2],
                    gb[la]
                )
            }

            return builder.getLines()
        }
    }
}