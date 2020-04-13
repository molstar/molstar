/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ChunkedArray } from '../chunked-array';

describe('Chunked Array', () => {
    it('creation', () => {
        const arr  = ChunkedArray.create<number, 2>(Array, 2, 2);
        ChunkedArray.add2(arr, 1, 2);
        ChunkedArray.add2(arr, 3, 4);
        expect(ChunkedArray.compact(arr)).toEqual([1, 2, 3, 4]);
    });

    it('initial', () => {
        const arr  = ChunkedArray.create(Int32Array, 2, 6, new Int32Array([1, 2, 3, 4]));
        ChunkedArray.add2(arr, 4, 3);
        ChunkedArray.add2(arr, 2, 1);
        ChunkedArray.add2(arr, 5, 6);
        expect(ChunkedArray.compact(arr)).toEqual(new Int32Array([4, 3, 2, 1, 5, 6]));
    });

    it('add many', () => {
        const arr  = ChunkedArray.create<number, 2>(Array, 2, 2);
        ChunkedArray.addMany(arr, [1, 2, 3, 4]);
        expect(ChunkedArray.compact(arr)).toEqual([1, 2, 3, 4]);
    });

    it('resize', () => {
        const arr  = ChunkedArray.create<number, 2>(Int32Array, 2, 2);
        ChunkedArray.add2(arr, 1, 2);
        ChunkedArray.add2(arr, 3, 4);
        ChunkedArray.add2(arr, 5, 6);
        ChunkedArray.add2(arr, 7, 8);
        ChunkedArray.add2(arr, 9, 10);
        expect(ChunkedArray.compact(arr)).toEqual(new Int32Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]));
    });

    it('resize-fraction', () => {
        const arr  = ChunkedArray.create<number, 2>(Int32Array, 2, 2.5);
        ChunkedArray.add2(arr, 1, 2);
        ChunkedArray.add2(arr, 3, 4);
        ChunkedArray.add2(arr, 5, 6);
        ChunkedArray.add2(arr, 7, 8);
        ChunkedArray.add2(arr, 9, 10);
        expect(arr.elementCount).toBe(5);
        expect(ChunkedArray.compact(arr)).toEqual(new Int32Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]));
    });
});