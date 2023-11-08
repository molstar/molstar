/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { canonicalJsonString, filterInPlace, range, sortObjectKeys } from '../utils';


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

    it('sortObjectKeys', async () => {
        expect(sortObjectKeys({})).toEqual({});

        const obj1 = { c: 1, b: 2, d: undefined, a: { f: undefined, e: 4 } };
        expect(sortObjectKeys(obj1)).toEqual(obj1);
        expect(Object.keys(sortObjectKeys(obj1))).toEqual(['a', 'b', 'c']);

        const obj2 = { c: [], x: null, b: false, a: undefined };
        expect(sortObjectKeys(obj2)).toEqual(obj2);
        expect(Object.keys(sortObjectKeys(obj2))).toEqual(['b', 'c', 'x']);
    });

    it('canonicalJsonString', async () => {
        expect(canonicalJsonString({})).toEqual('{}');

        const obj1 = { c: 1, b: 2, d: undefined, a: { x: null, f: undefined, e: 4 } };
        expect(canonicalJsonString(obj1)).toEqual('{"a":{"e":4,"x":null},"b":2,"c":1}');

        const obj2 = { c: [1, { p: 'P', q: undefined }, 0], x: null, b: false, a: undefined };
        expect(canonicalJsonString(obj2)).toEqual('{"b":false,"c":[1,{"p":"P"},0],"x":null}');
    });
});

