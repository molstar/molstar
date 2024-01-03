/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { filterInPlace, range } from '../array';


describe('filterInPlace', () => {
    it('filterInPlace works', async () => {
        expect(filterInPlace([], () => true)).toEqual([]);
        expect(filterInPlace([], () => false)).toEqual([]);
        expect(filterInPlace([1, 2, 3, 4, 5], () => true)).toEqual([1, 2, 3, 4, 5]);
        expect(filterInPlace([1, 2, 3, 4, 5], () => false)).toEqual([]);
        expect(filterInPlace([1, 2, 3, 4, 5], x => x % 2 === 0)).toEqual([2, 4]);
        expect(filterInPlace([1, 2, 3, 4, 5], x => x % 2 === 1)).toEqual([1, 3, 5]);
        expect(filterInPlace([1, 2, 3, 4, 5], x => x <= 2)).toEqual([1, 2]);
        expect(filterInPlace([1, 2, 3, 4, 5], x => x > 2)).toEqual([3, 4, 5]);
    });
    it('filterInPlace works in place', async () => {
        const array = [1, 2, 3, 4, 5];
        filterInPlace(array, x => x % 2 === 1);
        expect(array).toEqual([1, 3, 5]);
    });
    it('filterInPlace big data', async () => {
        const array = range(10 ** 5);
        const expectedResult = array.filter(x => x % 7 === 0);
        expect(filterInPlace(array, x => x % 7 === 0)).toEqual(expectedResult);
    });
});

