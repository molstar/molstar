/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { ChunkedArray } from 'mol-data/util';
import { Lines } from './lines';

export interface LinesBuilder {
    add(startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, group: number): void
    getLines(): Lines
}

export namespace LinesBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, lines?: Lines): LinesBuilder {
        const mappings = ChunkedArray.create(Float32Array, 2, chunkSize, lines ? lines.mappingBuffer.ref.value : initialCount);
        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, lines ? lines.groupBuffer.ref.value : initialCount);
        const indices = ChunkedArray.create(Uint32Array, 3, chunkSize * 3, lines ? lines.indexBuffer.ref.value : initialCount * 3);
        const starts = ChunkedArray.create(Float32Array, 3, chunkSize, lines ? lines.startBuffer.ref.value : initialCount);
        const ends = ChunkedArray.create(Float32Array, 3, chunkSize, lines ? lines.endBuffer.ref.value : initialCount);

        return {
            add: (startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, group: number) => {
                const offset = mappings.elementCount
                for (let i = 0; i < 4; ++i) {
                    ChunkedArray.add3(starts, startX, startY, startZ);
                    ChunkedArray.add3(ends, endX, endY, endZ);
                    ChunkedArray.add(groups, group);
                }
                ChunkedArray.add2(mappings, -1, 1);
                ChunkedArray.add2(mappings, -1, -1);
                ChunkedArray.add2(mappings, 1, 1);
                ChunkedArray.add2(mappings, 1, -1);
                // ChunkedArray.add3(indices, offset, offset + 1, offset + 2);
                // ChunkedArray.add3(indices, offset + 1, offset + 3, offset + 2);
                ChunkedArray.add3(indices, offset + 2, offset + 1, offset);
                ChunkedArray.add3(indices, offset + 2, offset + 3, offset + 1);
            },
            getLines: () => {
                const mb = ChunkedArray.compact(mappings, true) as Float32Array
                const ib = ChunkedArray.compact(indices, true) as Uint32Array
                const gb = ChunkedArray.compact(groups, true) as Float32Array
                const sb = ChunkedArray.compact(starts, true) as Float32Array
                const eb = ChunkedArray.compact(ends, true) as Float32Array
                console.log(indices.elementCount, mappings.elementCount, groups.elementCount)
                return {
                    kind: 'lines',
                    lineCount: indices.elementCount / 2,
                    mappingBuffer: lines ? ValueCell.update(lines.mappingBuffer, mb) : ValueCell.create(mb),
                    indexBuffer: lines ? ValueCell.update(lines.indexBuffer, ib) : ValueCell.create(ib),
                    groupBuffer: lines ? ValueCell.update(lines.groupBuffer, gb) : ValueCell.create(gb),
                    startBuffer: lines ? ValueCell.update(lines.startBuffer, sb) : ValueCell.create(sb),
                    endBuffer: lines ? ValueCell.update(lines.endBuffer, eb) : ValueCell.create(eb),
                }
            }
        }
    }
}