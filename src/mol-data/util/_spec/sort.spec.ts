/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { radixSort } from '../sort';

function randomUints(n: number, max: number) {
    const data = new Uint32Array(n);
    for (let i = 0; i < n; i++) data[i] = Math.floor(Math.random() * max);
    return data;
}

function expectSorted(data: ArrayLike<number>, reference: number[]) {
    expect(Array.from(data)).toEqual(reference.slice().sort((a, b) => a - b));
}

describe('radixSortUint', () => {
    it('sorts random values like a comparison sort', () => {
        const data = randomUints(1000, 2 ** 24);
        const reference = Array.from(data);
        expect(radixSort(data)).toBe(data);
        expectSorted(data, reference);
    });

    it('handles empty and single-element arrays', () => {
        expect(Array.from(radixSort(new Int32Array(0)))).toEqual([]);
        expect(Array.from(radixSort(Int32Array.of(42)))).toEqual([42]);
    });

    it('handles duplicates and zeros', () => {
        const data = Int32Array.of(5, 0, 3, 5, 0, 3, 3, 0);
        radixSort(data);
        expect(Array.from(data)).toEqual([0, 0, 0, 3, 3, 3, 5, 5]);
    });

    it('handles already sorted and reversed input', () => {
        const asc = Int32Array.of(1, 2, 3, 4, 5);
        expect(Array.from(radixSort(asc))).toEqual([1, 2, 3, 4, 5]);
        const desc = Int32Array.of(5, 4, 3, 2, 1);
        expect(Array.from(radixSort(desc))).toEqual([1, 2, 3, 4, 5]);
    });

    it('supports the full Uint32Array range', () => {
        const data = Uint32Array.of(4294967295, 0, 2147483648, 1, 2147483647);
        radixSort(data);
        expect(Array.from(data)).toEqual([0, 1, 2147483647, 2147483648, 4294967295]);
    });

    it('respects an explicit maxValue', () => {
        const data = randomUints(500, 255);
        const reference = Array.from(data);
        radixSort(data, 255);
        expectSorted(data, reference);
    });

    it('sorts only the given subarray view', () => {
        const data = Int32Array.of(9, 8, 3, 1, 2, 7, 6);
        radixSort(data.subarray(2, 5));
        expect(Array.from(data)).toEqual([9, 8, 1, 2, 3, 7, 6]);
    });

    it('works with an even and odd number of passes', () => {
        // 1 pass (max < 2^8), 2 passes (< 2^16), 3 passes (< 2^24), 4 passes (< 2^32)
        for (const max of [2 ** 8, 2 ** 16, 2 ** 24, 2 ** 32]) {
            const data = randomUints(300, max);
            const reference = Array.from(data);
            radixSort(data);
            expectSorted(data, reference);
        }
    });
});
