/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ChunkedArray } from '../../../mol-data/util';
import { Lines } from './lines';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Cage } from '../../../mol-geo/primitive/cage';

export interface LinesBuilder {
    add(startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, group: number): void
    addFixedCountDashes(start: Vec3, end: Vec3, segmentCount: number, group: number): void
    addFixedLengthDashes(start: Vec3, end: Vec3, segmentLength: number, group: number): void
    addCage(t: Mat4, cage: Cage, group: number): void
    getLines(): Lines
}

const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpDir = Vec3();

export namespace LinesBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, lines?: Lines): LinesBuilder {
        const mappings = ChunkedArray.create(Float32Array, 2, chunkSize, lines ? lines.mappingBuffer.ref.value : initialCount);
        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, lines ? lines.groupBuffer.ref.value : initialCount);
        const indices = ChunkedArray.create(Uint32Array, 3, chunkSize * 3, lines ? lines.indexBuffer.ref.value : initialCount * 3);
        const starts = ChunkedArray.create(Float32Array, 3, chunkSize, lines ? lines.startBuffer.ref.value : initialCount);
        const ends = ChunkedArray.create(Float32Array, 3, chunkSize, lines ? lines.endBuffer.ref.value : initialCount);

        const add = (startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, group: number) => {
            const offset = mappings.elementCount;
            for (let i = 0; i < 4; ++i) {
                ChunkedArray.add3(starts, startX, startY, startZ);
                ChunkedArray.add3(ends, endX, endY, endZ);
                ChunkedArray.add(groups, group);
            }
            ChunkedArray.add2(mappings, -1, 1);
            ChunkedArray.add2(mappings, -1, -1);
            ChunkedArray.add2(mappings, 1, 1);
            ChunkedArray.add2(mappings, 1, -1);
            ChunkedArray.add3(indices, offset, offset + 1, offset + 2);
            ChunkedArray.add3(indices, offset + 1, offset + 3, offset + 2);
        };

        const addFixedCountDashes = (start: Vec3, end: Vec3, segmentCount: number, group: number) => {
            const d = Vec3.distance(start, end);
            const s = Math.floor(segmentCount / 2);
            const step = 1 / segmentCount;

            Vec3.sub(tmpDir, end, start);
            for (let j = 0; j < s; ++j) {
                const f = step * (j * 2 + 1);
                Vec3.setMagnitude(tmpDir, tmpDir, d * f);
                Vec3.add(tmpVecA, start, tmpDir);
                Vec3.setMagnitude(tmpDir, tmpDir, d * step * ((j + 1) * 2));
                Vec3.add(tmpVecB, start, tmpDir);
                add(tmpVecA[0], tmpVecA[1], tmpVecA[2], tmpVecB[0], tmpVecB[1], tmpVecB[2], group);
            }
        };

        return {
            add,
            addFixedCountDashes,
            addFixedLengthDashes: (start: Vec3, end: Vec3, segmentLength: number, group: number) => {
                const d = Vec3.distance(start, end);
                addFixedCountDashes(start, end, d / segmentLength, group);
            },
            addCage: (t: Mat4, cage: Cage, group: number) => {
                const { vertices, edges } = cage;
                for (let i = 0, il = edges.length; i < il; i += 2) {
                    Vec3.fromArray(tmpVecA, vertices, edges[i] * 3);
                    Vec3.fromArray(tmpVecB, vertices, edges[i + 1] * 3);
                    Vec3.transformMat4(tmpVecA, tmpVecA, t);
                    Vec3.transformMat4(tmpVecB, tmpVecB, t);
                    add(tmpVecA[0], tmpVecA[1], tmpVecA[2], tmpVecB[0], tmpVecB[1], tmpVecB[2], group);
                }
            },
            getLines: () => {
                const mb = ChunkedArray.compact(mappings, true) as Float32Array;
                const ib = ChunkedArray.compact(indices, true) as Uint32Array;
                const gb = ChunkedArray.compact(groups, true) as Float32Array;
                const sb = ChunkedArray.compact(starts, true) as Float32Array;
                const eb = ChunkedArray.compact(ends, true) as Float32Array;
                return Lines.create(mb, ib, gb, sb, eb, indices.elementCount / 2);
            }
        };
    }
}