/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ChunkedArray } from '../../../mol-data/util';
import { Lines } from './lines';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Cage } from '../../../mol-geo/primitive/cage';

export interface LinesBuilder {
    add(startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, group: number): void
    addVec(start: Vec3, end: Vec3, group: number): void
    addFixedCountDashes(start: Vec3, end: Vec3, segmentCount: number, group: number): void
    addFixedLengthDashes(start: Vec3, end: Vec3, segmentLength: number, group: number): void
    addCage(t: Mat4, cage: Cage, group: number): void
    getLines(): Lines
}

const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpDir = Vec3();

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const caAdd = ChunkedArray.add;
const caAdd3 = ChunkedArray.add3;

export namespace LinesBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, lines?: Lines): LinesBuilder {
        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, lines ? lines.groupBuffer.ref.value : initialCount);
        const starts = ChunkedArray.create(Float32Array, 3, chunkSize, lines ? lines.startBuffer.ref.value : initialCount);
        const ends = ChunkedArray.create(Float32Array, 3, chunkSize, lines ? lines.endBuffer.ref.value : initialCount);

        const add = (startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, group: number) => {
            for (let i = 0; i < 4; ++i) {
                caAdd3(starts, startX, startY, startZ);
                caAdd3(ends, endX, endY, endZ);
                caAdd(groups, group);
            }
        };

        const addVec = (start: Vec3, end: Vec3, group: number) => {
            for (let i = 0; i < 4; ++i) {
                caAdd3(starts, start[0], start[1], start[2]);
                caAdd3(ends, end[0], end[1], end[2]);
                caAdd(groups, group);
            }
        };

        const addFixedCountDashes = (start: Vec3, end: Vec3, segmentCount: number, group: number) => {
            const d = Vec3.distance(start, end);
            const isOdd = segmentCount % 2 !== 0;
            const s = Math.floor((segmentCount + 1) / 2);
            const step = d / (segmentCount + 0.5);

            Vec3.setMagnitude(tmpDir, Vec3.sub(tmpDir, end, start), step);
            Vec3.copy(tmpVecA, start);
            for (let j = 0; j < s; ++j) {
                Vec3.add(tmpVecA, tmpVecA, tmpDir);
                if (isOdd && j === s - 1) {
                    Vec3.copy(tmpVecB, end);
                } else {
                    Vec3.add(tmpVecB, tmpVecA, tmpDir);
                }
                add(tmpVecA[0], tmpVecA[1], tmpVecA[2], tmpVecB[0], tmpVecB[1], tmpVecB[2], group);
                Vec3.add(tmpVecA, tmpVecA, tmpDir);
            }
        };

        return {
            add,
            addVec,
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
                const lineCount = groups.elementCount / 4;
                const gb = ChunkedArray.compact(groups, true) as Float32Array;
                const sb = ChunkedArray.compact(starts, true) as Float32Array;
                const eb = ChunkedArray.compact(ends, true) as Float32Array;
                const mb = lines && lineCount <= lines.lineCount ? lines.mappingBuffer.ref.value : new Float32Array(lineCount * 8);
                const ib = lines && lineCount <= lines.lineCount ? lines.indexBuffer.ref.value : new Uint32Array(lineCount * 6);
                if (!lines || lineCount > lines.lineCount) fillMappingAndIndices(lineCount, mb, ib);
                return Lines.create(mb, ib, gb, sb, eb, lineCount, lines);
            }
        };
    }
}

function fillMappingAndIndices(n: number, mb: Float32Array, ib: Uint32Array) {
    for (let i = 0; i < n; ++i) {
        const mo = i * 8;
        mb[mo] = -1; mb[mo + 1] = -1;
        mb[mo + 2] = 1; mb[mo + 3] = -1;
        mb[mo + 4] = -1; mb[mo + 5] = 1;
        mb[mo + 6] = 1; mb[mo + 7] = 1;
    }

    for (let i = 0; i < n; ++i) {
        const o = i * 4;
        const io = i * 6;
        ib[io] = o; ib[io + 1] = o + 1; ib[io + 2] = o + 2;
        ib[io + 3] = o + 1; ib[io + 4] = o + 3; ib[io + 5] = o + 2;
    }
}
