/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../iterator';

function iteratorToArray<T>(it: Iterator<T>): T[] {
    const ret = [];
    while (it.hasNext) {
        const v = it.move();
        ret[ret.length] = v;
    }
    return ret;
}

describe('basic iterators', () => {
    function check<T>(name: string, iter: Iterator<T>, expected: T[]) {
        it(name, () => {
            expect(iteratorToArray(iter)).toEqual(expected);
        });
    }

    check('empty', Iterator.Empty, []);
    check('singleton', Iterator.Value(10), [10]);
    check('array', Iterator.Array([1, 2, 3]), [1, 2, 3]);
    check('range', Iterator.Range(0, 3), [0, 1, 2, 3]);
    check('map', Iterator.map(Iterator.Range(0, 1), x => x + 1), [1, 2]);
    check('filter', Iterator.filter(Iterator.Range(0, 3), x => x >= 2), [2, 3]);
});