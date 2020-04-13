/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createRangeArray } from '../array';
import { makeBuckets } from '../buckets';

describe('buckets', () => {

    function reorder(order: ArrayLike<number>, data: any[]): any[] {
        const ret = [];
        for (const i of (order as number[])) ret[ret.length] = data[i];
        return ret;
    }

    it('full range', () => {
        const xs = [1, 1, 2, 2, 3, 1];
        const range = createRangeArray(0, xs.length - 1);
        const bs = makeBuckets(range, i => xs[i]);

        expect(reorder(range, xs)).toEqual([1, 1, 1, 2, 2, 3]);
        expect(Array.from(bs)).toEqual([0, 3, 5, 6]);
    });

    it('sort', () => {
        const xs = [3, 1, 2, 1, 2, 3];
        const range = createRangeArray(0, xs.length - 1);
        makeBuckets(range, i => xs[i], { sort: true });

        expect(reorder(range, xs)).toEqual([1, 1, 2, 2, 3, 3]);
    });

    it('subrange', () => {
        const xs = [2, 1, 2, 1, 2, 3, 1];
        const range = createRangeArray(0, xs.length - 1);
        const bs = makeBuckets(range, i => xs[i], { sort: false, start: 1, end: 5 });

        expect(reorder(range, xs)).toEqual([2, 1, 1, 2, 2, 3, 1]);
        expect(Array.from(bs)).toEqual([1, 3, 5]);
    });
});