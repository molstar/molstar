/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * A generic chunked array builder.
 *
 * When adding elements, the array growns by a specified number
 * of elements and no copying is done until ChunkedArray.compact
 * is called.
 */
export interface ChunkedArray<T> {
    creator: (size: number) => any;
    elementSize: number;
    chunkSize: number;
    current: any;
    currentIndex: number;

    parts: any[];
    elementCount: number;
}

export namespace ChunkedArray {
    export function is(x: any): x is ChunkedArray<any> {
        return x.creator && x.chunkSize;
    }

    export function add4<T>(array: ChunkedArray<T>, x: T, y: T, z: T, w: T) {
        if (array.currentIndex >= array.chunkSize) {
            array.currentIndex = 0;
            array.current = array.creator(array.chunkSize);
            array.parts[array.parts.length] = array.current;
        }

        array.current[array.currentIndex++] = x;
        array.current[array.currentIndex++] = y;
        array.current[array.currentIndex++] = z;
        array.current[array.currentIndex++] = w;
        return array.elementCount++;
    }

    export function add3<T>(array: ChunkedArray<T>, x: T, y: T, z: T) {
        if (array.currentIndex >= array.chunkSize) {
            array.currentIndex = 0;
            array.current = array.creator(array.chunkSize);
            array.parts[array.parts.length] = array.current;
        }

        array.current[array.currentIndex++] = x;
        array.current[array.currentIndex++] = y;
        array.current[array.currentIndex++] = z;
        return array.elementCount++;
    }

    export function add2<T>(array: ChunkedArray<T>, x: T, y: T) {
        if (array.currentIndex >= array.chunkSize) {
            array.currentIndex = 0;
            array.current = array.creator(array.chunkSize);
            array.parts[array.parts.length] = array.current;
        }

        array.current[array.currentIndex++] = x;
        array.current[array.currentIndex++] = y;
        return array.elementCount++;
    }

    export function add<T>(array: ChunkedArray<T>, x: T) {
        if (array.currentIndex >= array.chunkSize) {
            array.currentIndex = 0;
            array.current = array.creator(array.chunkSize);
            array.parts[array.parts.length] = array.current;
        }

        array.current[array.currentIndex++] = x;
        return array.elementCount++;
    }


    export function compact<T>(array: ChunkedArray<T>): T[] {
        const ret = array.creator(array.elementSize * array.elementCount)
        const offset = (array.parts.length - 1) * array.chunkSize
        let offsetInner = 0
        let part: any

        if (array.parts.length > 1) {
            if (array.parts[0].buffer) {
                for (let i = 0; i < array.parts.length - 1; i++) {
                    ret.set(array.parts[i], array.chunkSize * i);
                }
            } else {

                for (let i = 0; i < array.parts.length - 1; i++) {
                    offsetInner = array.chunkSize * i;
                    part = array.parts[i];

                    for (let j = 0; j < array.chunkSize; j++) {
                        ret[offsetInner + j] = part[j];
                    }
                }
            }
        }

        if (array.current.buffer && array.currentIndex >= array.chunkSize) {
            ret.set(array.current, array.chunkSize * (array.parts.length - 1));
        } else {
            for (let i = 0; i < array.currentIndex; i++) {
                ret[offset + i] = array.current[i];
            }
        }
        return ret as any;
    }

    export function create<T>(creator: (size: number) => any, chunkElementCount: number, elementSize: number): ChunkedArray<T> {
        chunkElementCount = chunkElementCount | 0;
        if (chunkElementCount <= 0) chunkElementCount = 1;

        let chunkSize = chunkElementCount * elementSize;
        let current = creator(chunkSize)

        return {
            elementSize,
            chunkSize,
            creator,
            current,
            parts: [current],
            currentIndex: 0,
            elementCount: 0
        } as ChunkedArray<T>
    }
}
