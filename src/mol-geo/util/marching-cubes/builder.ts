/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Fred Ludlow <fred.ludlow@gmail.com>
 */

import { ChunkedArray } from '../../../mol-data/util';
import { noop } from '../../../mol-util';
import { Mesh } from '../../geometry/mesh/mesh';
import { AllowedContours } from './tables';
import { LinesBuilder } from '../../geometry/lines/lines-builder';
import { Lines } from '../../geometry/lines/lines';

export interface MarchinCubesBuilder<T> {
    addVertex(x: number, y: number, z: number): number
    addNormal(x: number, y: number, z: number): void
    addGroup(group: number): void
    addTriangle(vertList: number[], a: number, b: number, c: number, edgeFilter: number): void
    get(): T
}

export function MarchinCubesMeshBuilder(vertexChunkSize: number, mesh?: Mesh): MarchinCubesBuilder<Mesh> {
    const triangleChunkSize = Math.min(1 << 16, vertexChunkSize * 4);

    const vertices = ChunkedArray.create(Float32Array, 3, vertexChunkSize, mesh && mesh.vertexBuffer.ref.value);
    const normals = ChunkedArray.create(Float32Array, 3, vertexChunkSize, mesh && mesh.normalBuffer.ref.value);
    const groups = ChunkedArray.create(Float32Array, 1, vertexChunkSize, mesh && mesh.groupBuffer.ref.value);
    const indices = ChunkedArray.create(Uint32Array, 3, triangleChunkSize, mesh && mesh.indexBuffer.ref.value);

    let vertexCount = 0;
    let triangleCount = 0;

    return {
        addVertex: (x: number, y: number, z: number) => {
            ++vertexCount;
            return ChunkedArray.add3(vertices, x, y, z );
        },
        addNormal: (x: number, y: number, z: number) => {
            ChunkedArray.add3(normals, x, y, z );
        },
        addGroup: (group: number) => {
            ChunkedArray.add(groups, group);
        },
        addTriangle: (vertList: number[], a: number, b: number, c: number) => {
            const i = vertList[a], j = vertList[b], k = vertList[c];
            // vertex indices <0 mean that the vertex was ignored and is not available
            // and hence we don't add a triangle when this occurs
            if (i >= 0 && j >= 0 && k >= 0) {
                ++triangleCount;
                ChunkedArray.add3(indices, i, j, k);
            }
        },
        get: () => {
            const vb = ChunkedArray.compact(vertices, true) as Float32Array;
            const nb = ChunkedArray.compact(normals, true) as Float32Array;
            const ib = ChunkedArray.compact(indices, true) as Uint32Array;
            const gb = ChunkedArray.compact(groups, true) as Float32Array;
            return Mesh.create(vb, ib, nb, gb, vertexCount, triangleCount, mesh);
        }
    };
}

export function MarchinCubesLinesBuilder(vertexChunkSize: number, lines?: Lines): MarchinCubesBuilder<Lines> {
    const vertices = ChunkedArray.create(Float32Array, 3, vertexChunkSize);
    const groups = ChunkedArray.create(Float32Array, 1, vertexChunkSize);
    const indices = ChunkedArray.create(Float32Array, 2, vertexChunkSize);

    let linesCount = 0;

    return {
        addVertex: (x: number, y: number, z: number) => {
            return ChunkedArray.add3(vertices, x, y, z);
        },
        addNormal: () => noop,
        addGroup: (group: number) => {
            ChunkedArray.add(groups, group);
        },
        addTriangle: (vertList: number[], a: number, b: number, c: number, edgeFilter: number) => {
            const i = vertList[a], j = vertList[b], k = vertList[c];
            // vertex indices <0 mean that the vertex was ignored and is not available
            // and hence we don't add a triangle when this occurs
            if (i >= 0 && j >= 0 && k >= 0) {
                if (AllowedContours[a][b] & edgeFilter) {
                    ++linesCount;
                    ChunkedArray.add2(indices, vertList[a], vertList[b]);
                }
                if (AllowedContours[b][c] & edgeFilter) {
                    ++linesCount;
                    ChunkedArray.add2(indices, vertList[b], vertList[c]);
                }
                if (AllowedContours[a][c] & edgeFilter) {
                    ++linesCount;
                    ChunkedArray.add2(indices, vertList[a], vertList[c]);
                }
            }
        },
        get: () => {
            const vb = ChunkedArray.compact(vertices, true) as Float32Array;
            const ib = ChunkedArray.compact(indices, true) as Uint32Array;
            const gb = ChunkedArray.compact(groups, true) as Float32Array;

            const builder = LinesBuilder.create(linesCount, linesCount / 10, lines);

            for (let i = 0; i < linesCount; ++i) {
                const la = ib[i * 2], lb = ib[i * 2 + 1];
                builder.add(
                    vb[la * 3], vb[la * 3 + 1], vb[la * 3 + 2],
                    vb[lb * 3], vb[lb * 3 + 1], vb[lb * 3 + 2],
                    gb[la]
                );
            }

            return builder.getLines();
        }
    };
}