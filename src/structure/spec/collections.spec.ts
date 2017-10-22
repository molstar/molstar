/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../collections/iterator'
import IntPair from '../collections/int-pair'
import * as Sort from '../collections/sort'
import OrderedSet from '../collections/ordered-set'
import LinkedIndex from '../collections/linked-index'
import MultiSet from '../collections/multi-set'


function iteratorToArray<T>(it: Iterator<T>): T[] {
    const ret = [];
    for (let v = it.move(); !it.done; v = it.move()) ret[ret.length] = v;
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
});

describe('int pair', () => {
    it('works', () => {
        const p = IntPair.zero();
        for (let i = 0; i < 10; i++) {
            for (let j = -10; j < 5; j++) {
                const t = IntPair.set1(i, j);
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
    const data0 = new Array(50);
    for (let i = 0; i < data0.length; i++) data0[i] = i;
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
    const data0 = new Array(50);
    for (let i = 0; i < data0.length; i++) data0[i] = i;
    const data1 = [1, 1, 2, 2, 3, 3, 4, 4, 4, 6, 6, 6];

    function test(name: string, data: any[], randomize: boolean) {
        it(name, () => {
            // [ 3, 1, 6, 4, 4, 6, 4, 2, 6, 1, 2, 3 ];
            if (randomize) {
                for (let i = 0; i < 10; i++) {
                    expect(Sort.sort(shuffleArray(data), 0, data.length, Sort.arrayLess, Sort.arraySwap)).toEqual(data);
                }
            } else {
                expect(Sort.sort([...data], 0, data.length, Sort.arrayLess, Sort.arraySwap)).toEqual(data);
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
                    expect(Sort.sort(shuffle(src, len, clone, swap), 0, len, cmp, swap)).toEqual(data);
                }
            } else {
                expect(Sort.sort(clone(src), 0, len, cmp, swap)).toEqual(data);
            }
        });
    }
    test('sorted', data, false);
    test('shuffled', data, true);
})

describe('range set', () => {
    function testEq(name: string, set: OrderedSet, expected: number[]) {
        it(name, () => {
            // copy the arrays to ensure "compatibility" between typed and native arrays
            expect(Array.prototype.slice.call(iteratorToArray(set.elements()))).toEqual(Array.prototype.slice.call(expected));
        });
    }

    const empty = OrderedSet.Empty;
    const singleton = OrderedSet.ofSingleton(10);
    const range = OrderedSet.ofRange(1, 4);
    const arr = OrderedSet.ofSortedArray([1, 3, 6]);

    testEq('empty', empty, []);
    testEq('singleton', singleton, [10]);
    testEq('range', range, [1, 2, 3, 4]);
    testEq('sorted array', arr, [1, 3, 6]);

    it('equality', () => {
        expect(OrderedSet.areEqual(empty, singleton)).toBe(false);
        expect(OrderedSet.areEqual(singleton, singleton)).toBe(true);
        expect(OrderedSet.areEqual(range, singleton)).toBe(false);
        expect(OrderedSet.areEqual(arr, OrderedSet.ofSortedArray([1, 3, 6]))).toBe(true);
        expect(OrderedSet.areEqual(arr, OrderedSet.ofSortedArray([1, 4, 6]))).toBe(false);
    });

    it('areIntersecting', () => {
        expect(OrderedSet.areIntersecting(range, arr)).toBe(true);
        expect(OrderedSet.areIntersecting(empty, empty)).toBe(true);
        expect(OrderedSet.areIntersecting(empty, singleton)).toBe(false);
        expect(OrderedSet.areIntersecting(empty, range)).toBe(false);
        expect(OrderedSet.areIntersecting(empty, arr)).toBe(false);
    });

    it('isSubset', () => {
        expect(OrderedSet.isSubset(singleton, empty)).toBe(true);
        expect(OrderedSet.isSubset(range, empty)).toBe(true);
        expect(OrderedSet.isSubset(arr, empty)).toBe(true);
        expect(OrderedSet.isSubset(empty, empty)).toBe(true);
        expect(OrderedSet.isSubset(empty, singleton)).toBe(false);
        expect(OrderedSet.isSubset(empty, range)).toBe(false);
        expect(OrderedSet.isSubset(empty, arr)).toBe(false);

        expect(OrderedSet.isSubset(singleton, range)).toBe(false);
        expect(OrderedSet.isSubset(range, OrderedSet.ofRange(2, 3))).toBe(true);
        expect(OrderedSet.isSubset(arr, range)).toBe(false);
        expect(OrderedSet.isSubset(arr, arr)).toBe(true);
        expect(OrderedSet.isSubset(arr, OrderedSet.ofSortedArray([1, 3]))).toBe(true);
        expect(OrderedSet.isSubset(arr, OrderedSet.ofSortedArray([1, 3, 7]))).toBe(false);
        expect(OrderedSet.isSubset(arr, OrderedSet.ofSortedArray([1, 3, 10, 45]))).toBe(false);
    });

    it('access/membership', () => {
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
    });

    testEq('union ES', OrderedSet.union(empty, singleton), [10]);
    testEq('union ER', OrderedSet.union(empty, range), [1, 2, 3, 4]);
    testEq('union EA', OrderedSet.union(empty, arr), [1, 3, 6]);
    testEq('union SS', OrderedSet.union(singleton, OrderedSet.ofSingleton(16)), [10, 16]);
    testEq('union SR', OrderedSet.union(range, singleton), [1, 2, 3, 4, 10]);
    testEq('union SA', OrderedSet.union(arr, singleton), [1, 3, 6, 10]);
    testEq('union SA1', OrderedSet.union(arr, OrderedSet.ofSingleton(3)), [1, 3, 6]);
    testEq('union RR', OrderedSet.union(range, range), [1, 2, 3, 4]);
    testEq('union RR1', OrderedSet.union(range, OrderedSet.ofRange(6, 7)), [1, 2, 3, 4, 6, 7]);
    testEq('union RR2', OrderedSet.union(range, OrderedSet.ofRange(3, 5)), [1, 2, 3, 4, 5]);
    testEq('union RA', OrderedSet.union(range, arr), [1, 2, 3, 4, 6]);
    testEq('union AA', OrderedSet.union(arr, OrderedSet.ofSortedArray([2, 4, 6, 7])), [1, 2, 3, 4, 6, 7]);
    testEq('union AA1', OrderedSet.union(arr, OrderedSet.ofSortedArray([2, 3, 4, 6, 7])), [1, 2, 3, 4, 6, 7]);

    testEq('intersect ES', OrderedSet.intersect(empty, singleton), []);
    testEq('intersect ER', OrderedSet.intersect(empty, range), []);
    testEq('intersect EA', OrderedSet.intersect(empty, arr), []);
    testEq('intersect SS', OrderedSet.intersect(singleton, OrderedSet.ofSingleton(16)), []);
    testEq('intersect SS', OrderedSet.intersect(singleton, singleton), [10]);
    testEq('intersect SR', OrderedSet.intersect(range, singleton), []);
    testEq('intersect RR', OrderedSet.intersect(range, range), [1, 2, 3, 4]);
    testEq('intersect RR2', OrderedSet.intersect(range, OrderedSet.ofRange(3, 5)), [3, 4]);
    testEq('intersect RA', OrderedSet.intersect(range, arr), [1, 3]);
    testEq('intersect AA', OrderedSet.intersect(arr, OrderedSet.ofSortedArray([2, 3, 4, 6, 7])), [3, 6]);
});

describe('linked-index', () => {
    it('initial state', () => {
        const index = LinkedIndex(2);
        expect(index.head).toBe(0);
        expect(index.has(0)).toBe(true);
        expect(index.has(1)).toBe(true);
    });

    it('singleton', () => {
        const index = LinkedIndex(1);
        expect(index.head).toBe(0);
        expect(index.has(0)).toBe(true);
        index.remove(0);
        expect(index.head).toBe(-1);
        expect(index.has(0)).toBe(false);
    });

    it('remove 0', () => {
        const index = LinkedIndex(2);
        index.remove(0);
        expect(index.head).toBe(1);
        expect(index.has(0)).toBe(false);
        expect(index.has(1)).toBe(true);
    });

    it('remove 1', () => {
        const index = LinkedIndex(2);
        index.remove(1);
        expect(index.head).toBe(0);
        expect(index.has(0)).toBe(true);
        expect(index.has(1)).toBe(false);
    });

    it('remove 01', () => {
        const index = LinkedIndex(2);
        index.remove(0);
        index.remove(1);
        expect(index.head).toBe(-1);
        expect(index.has(0)).toBe(false);
        expect(index.has(1)).toBe(false);
    });
});

describe('multiset', () => {
    const p = (i: number, j: number) => IntPair.create(i, j);
    const r = (i: number, j: number) => IntPair.set1(i, j);

    function iteratorPairsToArray(it: Iterator<IntPair>): IntPair[] {
        const ret = [];
        for (let v = it.move(); !it.done; v = it.move()) ret[ret.length] = IntPair.create(v.fst, v.snd);
        return ret;
    }

    it('singleton pair', () => {
        const set = MultiSet.create(p(10, 11));
        expect(iteratorPairsToArray(MultiSet.values(set))).toEqual([p(10, 11)]);
    });

    it('singleton number', () => {
        const set = MultiSet.create(r(10, 11));
        expect(iteratorPairsToArray(MultiSet.values(set))).toEqual([p(10, 11)]);
    });

    it('multi', () => {
        const set = MultiSet.create({
            1: OrderedSet.ofSortedArray([4, 6, 7]),
            3: OrderedSet.ofRange(0, 1),
        });
        expect(iteratorPairsToArray(MultiSet.values(set))).toEqual([p(1, 4), p(1, 6), p(1, 7), p(3, 0), p(3, 1)]);
    });

    it('packed pairs', () => {
        const set = MultiSet.create([r(1, 3), r(0, 1), r(0, 6), r(0, 2)]);
        expect(iteratorPairsToArray(MultiSet.values(set))).toEqual([p(0, 1), p(0, 2), p(0, 6), p(1, 3)]);
    });

});