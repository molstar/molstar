/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { ChunkedArray } from 'mol-data/util';
import { Point } from './point';

export interface PointBuilder {
    add(x: number, y: number, z: number, group: number): void
    getPoint(): Point
}

export namespace PointBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, point?: Point): PointBuilder {
        const vertices = ChunkedArray.create(Float32Array, 3, chunkSize, point ? point.vertexBuffer.ref.value : initialCount);
        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, point ? point.groupBuffer.ref.value : initialCount);

        return {
            add: (x: number, y: number, z: number, group: number) => {
                ChunkedArray.add3(vertices, x, y, z);
                ChunkedArray.add(groups, group);
            },
            getPoint: () => {
                const vb = ChunkedArray.compact(vertices, true) as Float32Array
                const gb = ChunkedArray.compact(groups, true) as Float32Array
                return {
                    kind: 'point',
                    vertexCount: vertices.elementCount,
                    vertexBuffer: point ? ValueCell.update(point.vertexBuffer, vb) : ValueCell.create(vb),
                    groupBuffer: point ? ValueCell.update(point.groupBuffer, gb) : ValueCell.create(gb),
                }
            }
        }
    }
}