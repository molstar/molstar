/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../collections/iterator'
import IntPair from '../collections/int-pair'
import * as Sort from '../collections/sort'
import RangeSet from '../collections/range-set'

describe('basic iterators', () => {
    function check<T>(name: string, iter: Iterator<T>, expected: T[]) {
        it(name, () => {
            expect(Iterator.toArray(iter)).toEqual(expected);
        });
    }

    check('empty', Iterator.Empty, []);
    check('singleton', Iterator.Value(10), [10]);
    check('array', Iterator.Array([1, 2, 3]), [1, 2, 3]);
    check('range', Iterator.Range(0, 3), [0, 1, 2, 3]);
});

describe('int pair', () => {
    it('works', () => {
        const p = IntPair.zero();
        for (let i = 0; i < 10; i++) {
            for (let j = -10; j < 5; j++) {
                const t = IntPair.set(i, j);
                IntPair.get(t, p);
                expect(p.fst).toBe(i);
                expect(p.snd).toBe(j);
            }
        }
    })
})

function shuffle<T>(data: T, len: number, clone: (s: T) => T, swap: Sort.Swapper = Sort.arraySwap) {
    const a = clone(data);
    for (let i = len - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        swap(a, i, j);
    }
    return a;
}

function shuffleArray(data: any[]) {
    return shuffle(data, data.length, t => [...t]);
}

describe('qsort-array asc', () => {
    const data0 = [1, 2, 3, 4, 5, 6, 7, 8, 9];
    const data1 = [1, 1, 2, 2, 3, 3, 4, 4, 4, 6, 6, 6];

    function test(name: string, data: any[], randomize: boolean) {
        it(name, () => {
            // [ 3, 1, 6, 4, 4, 6, 4, 2, 6, 1, 2, 3 ];
            if (randomize) {
                for (let i = 0; i < 10; i++) {
                    expect(Sort.sortArray(shuffleArray(data))).toEqual(data);
                }
            } else {
                expect(Sort.sortArray([...data])).toEqual(data);
            }
        });
    }
    test('uniq', data0, false);
    test('uniq shuffle', data0, true);
    test('rep', data1, false);
    test('rep shuffle', data1, true);
})

describe('qsort-array generic', () => {
    const data0 = [1, 2, 3, 4, 5, 6, 7, 8, 9];
    const data1 = [1, 1, 2, 2, 3, 3, 4, 4, 4, 6, 6, 6];

    function test(name: string, data: any[], randomize: boolean) {
        it(name, () => {
            // [ 3, 1, 6, 4, 4, 6, 4, 2, 6, 1, 2, 3 ];
            if (randomize) {
                for (let i = 0; i < 10; i++) {
                    expect(Sort.sort(shuffleArray(data), data.length, Sort.arrayLess, Sort.arraySwap)).toEqual(data);
                }
            } else {
                expect(Sort.sort([...data], data.length, Sort.arrayLess, Sort.arraySwap)).toEqual(data);
            }
        });
    }
    test('uniq', data0, false);
    test('uniq shuffle', data0, true);
    test('rep', data1, false);
    test('rep shuffle', data1, true);
})

describe('qsort-dual array', () => {
    const len = 3;
    const data = { xs: [0, 1, 2], ys: ['x', 'y', 'z'] };

    const cmp: Sort.Comparer<typeof data> = (data, i, j) => data.xs[i] - data.xs[j];
    const swap: Sort.Swapper<typeof data> = (data, i, j) => { Sort.arraySwap(data.xs, i, j); Sort.arraySwap(data.ys, i, j); }
    const clone = (d: typeof data) => ({ xs: [...d.xs], ys: [...d.ys] })
 
    function test(name: string, src: typeof data, randomize: boolean) {
        it(name, () => {
            // [ 3, 1, 6, 4, 4, 6, 4, 2, 6, 1, 2, 3 ];
            if (randomize) {
                for (let i = 0; i < 10; i++) {
                    expect(Sort.sort(shuffle(src, len, clone, swap), len, cmp, swap)).toEqual(data);
                }
            } else {
                expect(Sort.sort(clone(src), len, cmp, swap)).toEqual(data);
            }
        });
    }
    test('sorted', data, false);
    test('shuffled', data, true);
})

describe('range set', () => {
    function testEq(name: string, set: RangeSet, expected: number[]) {
        it(name, () => {
            // copy the arrays to ensure "compatibility" between typed and native arrays
            expect(Array.prototype.slice.call(Iterator.toArray(set.elements()))).toEqual(Array.prototype.slice.call(expected));
        });
    }

    const empty = RangeSet.Empty;
    const singleton = RangeSet.ofSingleton(10);
    const range = RangeSet.ofRange(1, 4);
    const arr = RangeSet.ofSortedArray([1, 3, 6]);

    testEq('empty', empty, []);
    testEq('singleton', singleton, [10]);
    testEq('range', range, [1, 2, 3, 4]);
    testEq('sorted array', arr, [1, 3, 6]);

    expect(empty.has(10)).toBe(false);
    expect(empty.indexOf(10)).toBe(-1);

    expect(singleton.has(10)).toBe(true);
    expect(singleton.has(11)).toBe(false);
    expect(singleton.indexOf(10)).toBe(0);
    expect(singleton.indexOf(11)).toBe(-1);

    expect(range.has(4)).toBe(true);
    expect(range.has(5)).toBe(false);
    expect(range.indexOf(4)).toBe(3);
    expect(range.indexOf(11)).toBe(-1);

    expect(arr.has(3)).toBe(true);
    expect(arr.has(4)).toBe(false);
    expect(arr.indexOf(3)).toBe(1);
    expect(arr.indexOf(11)).toBe(-1);

    testEq('union ES', RangeSet.union(empty, singleton), [10]);
    testEq('union ER', RangeSet.union(empty, range), [1, 2, 3, 4]);
    testEq('union EA', RangeSet.union(empty, arr), [1, 3, 6]);
    testEq('union SS', RangeSet.union(singleton, RangeSet.ofSingleton(16)), [10, 16]);
    testEq('union SR', RangeSet.union(range, singleton), [1, 2, 3, 4, 10]);
    testEq('union SA', RangeSet.union(arr, singleton), [1, 3, 6, 10]);
    testEq('union SA1', RangeSet.union(arr, RangeSet.ofSingleton(3)), [1, 3, 6]);
    testEq('union RR', RangeSet.union(range, range), [1, 2, 3, 4]);
    testEq('union RR1', RangeSet.union(range, RangeSet.ofRange(6, 7)), [1, 2, 3, 4, 6, 7]);
    testEq('union RR2', RangeSet.union(range, RangeSet.ofRange(3, 5)), [1, 2, 3, 4, 5]);
    testEq('union RA', RangeSet.union(range, arr), [1, 2, 3, 4, 6]);
    testEq('union AA', RangeSet.union(arr, RangeSet.ofSortedArray([2, 4, 6, 7])), [1, 2, 3, 4, 6, 7]);
    testEq('union AA1', RangeSet.union(arr, RangeSet.ofSortedArray([2, 3, 4, 6, 7])), [1, 2, 3, 4, 6, 7]);

    testEq('intersect ES', RangeSet.intersect(empty, singleton), []);
    testEq('intersect ER', RangeSet.intersect(empty, range), []);
    testEq('intersect EA', RangeSet.intersect(empty, arr), []);
    testEq('intersect SS', RangeSet.intersect(singleton, RangeSet.ofSingleton(16)), []);
    testEq('intersect SS', RangeSet.intersect(singleton, singleton), [10]);
    testEq('intersect SR', RangeSet.intersect(range, singleton), []);
    testEq('intersect RR', RangeSet.intersect(range, range), [1, 2, 3, 4]);
    testEq('intersect RR2', RangeSet.intersect(range, RangeSet.ofRange(3, 5)), [3, 4]);
    testEq('intersect RA', RangeSet.intersect(range, arr), [1, 3]);
    testEq('intersect AA', RangeSet.intersect(arr, RangeSet.ofSortedArray([2, 3, 4, 6, 7])), [3, 6]);
});