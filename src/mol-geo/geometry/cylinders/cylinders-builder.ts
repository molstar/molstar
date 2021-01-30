/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ChunkedArray } from '../../../mol-data/util';
import { Cylinders } from './cylinders';
import { Vec3 } from '../../../mol-math/linear-algebra';

export interface CylindersBuilder {
    add(startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, radiusScale: number, topCap: boolean, bottomCap: boolean, group: number): void
    addFixedCountDashes(start: Vec3, end: Vec3, segmentCount: number, radiusScale: number, topCap: boolean, bottomCap: boolean, group: number): void
    addFixedLengthDashes(start: Vec3, end: Vec3, segmentLength: number, radiusScale: number, topCap: boolean, bottomCap: boolean, group: number): void
    getCylinders(): Cylinders
}

const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpDir = Vec3();

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const caAdd = ChunkedArray.add;
const caAdd3 = ChunkedArray.add3;

export namespace CylindersBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, cylinders?: Cylinders): CylindersBuilder {
        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, cylinders ? cylinders.groupBuffer.ref.value : initialCount);
        const starts = ChunkedArray.create(Float32Array, 3, chunkSize, cylinders ? cylinders.startBuffer.ref.value : initialCount);
        const ends = ChunkedArray.create(Float32Array, 3, chunkSize, cylinders ? cylinders.endBuffer.ref.value : initialCount);
        const scales = ChunkedArray.create(Float32Array, 1, chunkSize, cylinders ? cylinders.scaleBuffer.ref.value : initialCount);
        const caps = ChunkedArray.create(Float32Array, 1, chunkSize, cylinders ? cylinders.capBuffer.ref.value : initialCount);

        const add = (startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, radiusScale: number, topCap: boolean, bottomCap: boolean, group: number) => {
            for (let i = 0; i < 6; ++i) {
                caAdd3(starts, startX, startY, startZ);
                caAdd3(ends, endX, endY, endZ);
                caAdd(groups, group);
                caAdd(scales, radiusScale);
                caAdd(caps, (topCap ? 1 : 0) + (bottomCap ? 2 : 0));
            }
        };

        const addFixedCountDashes = (start: Vec3, end: Vec3, segmentCount: number, radiusScale: number, topCap: boolean, bottomCap: boolean, group: number) => {
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
                add(tmpVecA[0], tmpVecA[1], tmpVecA[2], tmpVecB[0], tmpVecB[1], tmpVecB[2], radiusScale, topCap, bottomCap, group);
            }
        };

        return {
            add,
            addFixedCountDashes,
            addFixedLengthDashes: (start: Vec3, end: Vec3, segmentLength: number, radiusScale: number, topCap: boolean, bottomCap: boolean, group: number) => {
                const d = Vec3.distance(start, end);
                addFixedCountDashes(start, end, d / segmentLength, radiusScale, topCap, bottomCap, group);
            },
            getCylinders: () => {
                const cylinderCount = groups.elementCount / 6;
                const gb = ChunkedArray.compact(groups, true) as Float32Array;
                const sb = ChunkedArray.compact(starts, true) as Float32Array;
                const eb = ChunkedArray.compact(ends, true) as Float32Array;
                const ab = ChunkedArray.compact(scales, true) as Float32Array;
                const cb = ChunkedArray.compact(caps, true) as Float32Array;
                const mb = cylinders && cylinderCount <= cylinders.cylinderCount ? cylinders.mappingBuffer.ref.value : new Float32Array(cylinderCount * 18);
                const ib = cylinders && cylinderCount <= cylinders.cylinderCount ? cylinders.indexBuffer.ref.value : new Uint32Array(cylinderCount * 12);
                if (!cylinders || cylinderCount > cylinders.cylinderCount) fillMappingAndIndices(cylinderCount, mb, ib);
                return Cylinders.create(mb, ib, gb, sb, eb, ab, cb, cylinderCount, cylinders);
            }
        };
    }
}

function fillMappingAndIndices(n: number, mb: Float32Array, ib: Uint32Array) {
    for (let i = 0; i < n; ++i) {
        const mo = i * 18;
        mb[mo] = -1; mb[mo + 1] = 1; mb[mo + 2] = -1;
        mb[mo + 3] = -1; mb[mo + 4] = -1; mb[mo + 5] = -1;
        mb[mo + 6] = 1; mb[mo + 7] = 1; mb[mo + 8] = -1;
        mb[mo + 9] = 1; mb[mo + 10] = 1; mb[mo + 11] = 1;
        mb[mo + 12] = 1; mb[mo + 13] = -1; mb[mo + 14] = -1;
        mb[mo + 15] = 1; mb[mo + 16] = -1; mb[mo + 17] = 1;
    }

    for (let i = 0; i < n; ++i) {
        const o = i * 6;
        const io = i * 12;
        ib[io] = o; ib[io + 1] = o + 1; ib[io + 2] = o + 2;
        ib[io + 3] = o + 1; ib[io + 4] = o + 4; ib[io + 5] = o + 2;
        ib[io + 6] = o + 2; ib[io + 7] = o + 4; ib[io + 8] = o + 3;
        ib[io + 9] = o + 4; ib[io + 10] = o + 5; ib[io + 11] = o + 3;
    }
}
