/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { canonicalJsonString, sortObjectKeys } from '../object';


describe('object utils', () => {
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

