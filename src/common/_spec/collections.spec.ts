/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../collections/iterator'
import IntTuple from '../collections/int-tuple'
import * as Sort from '../collections/sort'
import OrderedSet from '../collections/ordered-set'
import LinkedIndex from '../collections/linked-index'

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
        const p = IntTuple.zero();
        for (let i = 0; i < 10; i++) {
            for (let j = -10; j < 5; j++) {
                const t = IntTuple.pack(i, j);
                IntTuple.unpack(t, p);
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

describe('ordered set', () => {
    function ordSetToArray(set: OrderedSet) {
        const ret = [];
        for (let i = 0, _i = OrderedSet.size(set); i < _i; i++) ret.push(OrderedSet.getAt(set, i));
        return ret;
    }

    function testEq(name: string, set: OrderedSet, expected: number[]) {
        it(name, () => {
            // copy the arrays to ensure "compatibility" between typed and native arrays
            expect(Array.prototype.slice.call(ordSetToArray(set))).toEqual(Array.prototype.slice.call(expected));
        });
    }

    const empty = OrderedSet.Empty;
    const singleton10 = OrderedSet.ofSingleton(10);
    const range1_4 = OrderedSet.ofRange(1, 4);
    const arr136 = OrderedSet.ofSortedArray([1, 3, 6]);

    testEq('empty', empty, []);
    testEq('singleton', singleton10, [10]);
    testEq('range', range1_4, [1, 2, 3, 4]);
    testEq('sorted array', arr136, [1, 3, 6]);

    it('equality', () => {
        expect(OrderedSet.areEqual(empty, singleton10)).toBe(false);
        expect(OrderedSet.areEqual(singleton10, singleton10)).toBe(true);
        expect(OrderedSet.areEqual(range1_4, singleton10)).toBe(false);
        expect(OrderedSet.areEqual(arr136, OrderedSet.ofSortedArray([1, 3, 6]))).toBe(true);
        expect(OrderedSet.areEqual(arr136, OrderedSet.ofSortedArray([1, 4, 6]))).toBe(false);
    });

    it('areIntersecting', () => {
        expect(OrderedSet.areIntersecting(range1_4, arr136)).toBe(true);
        expect(OrderedSet.areIntersecting(empty, empty)).toBe(true);
        expect(OrderedSet.areIntersecting(empty, singleton10)).toBe(false);
        expect(OrderedSet.areIntersecting(empty, range1_4)).toBe(false);
        expect(OrderedSet.areIntersecting(empty, arr136)).toBe(false);
    });

    it('isSubset', () => {
        expect(OrderedSet.isSubset(singleton10, empty)).toBe(true);
        expect(OrderedSet.isSubset(range1_4, empty)).toBe(true);
        expect(OrderedSet.isSubset(arr136, empty)).toBe(true);
        expect(OrderedSet.isSubset(empty, empty)).toBe(true);
        expect(OrderedSet.isSubset(empty, singleton10)).toBe(false);
        expect(OrderedSet.isSubset(empty, range1_4)).toBe(false);
        expect(OrderedSet.isSubset(empty, arr136)).toBe(false);

        expect(OrderedSet.isSubset(singleton10, range1_4)).toBe(false);
        expect(OrderedSet.isSubset(range1_4, OrderedSet.ofRange(2, 3))).toBe(true);
        expect(OrderedSet.isSubset(arr136, range1_4)).toBe(false);
        expect(OrderedSet.isSubset(arr136, arr136)).toBe(true);
        expect(OrderedSet.isSubset(arr136, OrderedSet.ofSortedArray([1, 3]))).toBe(true);
        expect(OrderedSet.isSubset(arr136, OrderedSet.ofSortedArray([1, 3, 7]))).toBe(false);
        expect(OrderedSet.isSubset(OrderedSet.ofSortedArray([0, 1, 2, 3, 7, 10]), OrderedSet.ofSortedArray([1, 3, 7]))).toBe(true);
        expect(OrderedSet.isSubset(arr136, OrderedSet.ofSortedArray([1, 3, 10, 45]))).toBe(false);
    });

    it('access/membership', () => {
        expect(OrderedSet.has(empty, 10)).toBe(false);
        expect(OrderedSet.indexOf(empty, 10)).toBe(-1);

        expect(OrderedSet.has(singleton10, 10)).toBe(true);
        expect(OrderedSet.has(singleton10, 11)).toBe(false);
        expect(OrderedSet.indexOf(singleton10, 10)).toBe(0);
        expect(OrderedSet.indexOf(singleton10, 11)).toBe(-1);

        expect(OrderedSet.has(range1_4, 4)).toBe(true);
        expect(OrderedSet.has(range1_4, 5)).toBe(false);
        expect(OrderedSet.indexOf(range1_4, 4)).toBe(3);
        expect(OrderedSet.indexOf(range1_4, 11)).toBe(-1);

        expect(OrderedSet.has(arr136, 3)).toBe(true);
        expect(OrderedSet.has(arr136, 4)).toBe(false);
        expect(OrderedSet.indexOf(arr136, 3)).toBe(1);
        expect(OrderedSet.indexOf(arr136, 11)).toBe(-1);
    });

    it('interval range', () => {
        expect(OrderedSet.getIntervalRange(empty, 9, 11)).toEqual({ start: 0, end: 0 });
        expect(OrderedSet.getIntervalRange(empty, -9, -6)).toEqual({ start: 0, end: 0 });
        expect(OrderedSet.getIntervalRange(singleton10, 9, 11)).toEqual({ start: 0, end: 1 });
        expect(OrderedSet.getIntervalRange(range1_4, 2, 3)).toEqual({ start: 1, end: 3 });
        expect(OrderedSet.getIntervalRange(range1_4, -10, 2)).toEqual({ start: 0, end: 2 });
        expect(OrderedSet.getIntervalRange(range1_4, -10, 20)).toEqual({ start: 0, end: 4 });
        expect(OrderedSet.getIntervalRange(range1_4, 3, 20)).toEqual({ start: 2, end: 4 });
        expect(OrderedSet.getIntervalRange(arr136, 0, 1)).toEqual({ start: 0, end: 1 });
        expect(OrderedSet.getIntervalRange(arr136, 0, 3)).toEqual({ start: 0, end: 2 });
        expect(OrderedSet.getIntervalRange(arr136, 0, 4)).toEqual({ start: 0, end: 2 });
        expect(OrderedSet.getIntervalRange(arr136, 2, 4)).toEqual({ start: 1, end: 2 });
        expect(OrderedSet.getIntervalRange(arr136, 2, 7)).toEqual({ start: 1, end: 3 });
    })

    testEq('union ES', OrderedSet.union(empty, singleton10), [10]);
    testEq('union ER', OrderedSet.union(empty, range1_4), [1, 2, 3, 4]);
    testEq('union EA', OrderedSet.union(empty, arr136), [1, 3, 6]);
    testEq('union SS', OrderedSet.union(singleton10, OrderedSet.ofSingleton(16)), [10, 16]);
    testEq('union SR', OrderedSet.union(range1_4, singleton10), [1, 2, 3, 4, 10]);
    testEq('union SA', OrderedSet.union(arr136, singleton10), [1, 3, 6, 10]);
    testEq('union SA1', OrderedSet.union(arr136, OrderedSet.ofSingleton(3)), [1, 3, 6]);
    testEq('union RR', OrderedSet.union(range1_4, range1_4), [1, 2, 3, 4]);
    testEq('union RR1', OrderedSet.union(range1_4, OrderedSet.ofRange(6, 7)), [1, 2, 3, 4, 6, 7]);
    testEq('union RR2', OrderedSet.union(range1_4, OrderedSet.ofRange(3, 5)), [1, 2, 3, 4, 5]);
    testEq('union RA', OrderedSet.union(range1_4, arr136), [1, 2, 3, 4, 6]);
    testEq('union AA', OrderedSet.union(arr136, OrderedSet.ofSortedArray([2, 4, 6, 7])), [1, 2, 3, 4, 6, 7]);
    testEq('union AA1', OrderedSet.union(arr136, OrderedSet.ofSortedArray([2, 3, 4, 6, 7])), [1, 2, 3, 4, 6, 7]);
    testEq('union AA2', OrderedSet.union(arr136, OrderedSet.ofSortedArray([2, 4, 5, 6, 7])), [1, 2, 3, 4, 5, 6, 7]);

    testEq('intersect ES', OrderedSet.intersect(empty, singleton10), []);
    testEq('intersect ER', OrderedSet.intersect(empty, range1_4), []);
    testEq('intersect EA', OrderedSet.intersect(empty, arr136), []);
    testEq('intersect SS', OrderedSet.intersect(singleton10, OrderedSet.ofSingleton(16)), []);
    testEq('intersect SS1', OrderedSet.intersect(singleton10, singleton10), [10]);
    testEq('intersect SR', OrderedSet.intersect(range1_4, singleton10), []);
    testEq('intersect RR', OrderedSet.intersect(range1_4, range1_4), [1, 2, 3, 4]);
    testEq('intersect RR2', OrderedSet.intersect(range1_4, OrderedSet.ofRange(3, 5)), [3, 4]);
    testEq('intersect RA', OrderedSet.intersect(range1_4, arr136), [1, 3]);
    testEq('intersect AA', OrderedSet.intersect(arr136, OrderedSet.ofSortedArray([2, 3, 4, 6, 7])), [3, 6]);

    testEq('subtract ES', OrderedSet.subtract(empty, singleton10), []);
    testEq('subtract ER', OrderedSet.subtract(empty, range1_4), []);
    testEq('subtract EA', OrderedSet.subtract(empty, arr136), []);
    testEq('subtract SS', OrderedSet.subtract(singleton10, OrderedSet.ofSingleton(16)), [10]);
    testEq('subtract SS1', OrderedSet.subtract(singleton10, singleton10), []);
    testEq('subtract SR', OrderedSet.subtract(range1_4, singleton10), [1, 2, 3, 4]);
    testEq('subtract SR1', OrderedSet.subtract(range1_4, OrderedSet.ofSingleton(4)), [1, 2, 3]);
    testEq('subtract SR2', OrderedSet.subtract(range1_4, OrderedSet.ofSingleton(3)), [1, 2, 4]);
    testEq('subtract RR', OrderedSet.subtract(range1_4, range1_4), []);
    testEq('subtract RR1', OrderedSet.subtract(range1_4, OrderedSet.ofRange(3, 5)), [1, 2]);

    testEq('subtract RA', OrderedSet.subtract(range1_4, arr136), [2, 4]);
    testEq('subtract RA1', OrderedSet.subtract(range1_4, OrderedSet.ofSortedArray([0, 1, 2, 3, 4, 7])), []);
    testEq('subtract RA2', OrderedSet.subtract(range1_4, OrderedSet.ofSortedArray([0, 2, 3])), [1, 4]);

    testEq('subtract AR', OrderedSet.subtract(arr136, range1_4), [6]);
    testEq('subtract AR1', OrderedSet.subtract(arr136, OrderedSet.ofRange(0, 10)), []);
    testEq('subtract AR1', OrderedSet.subtract(arr136, OrderedSet.ofRange(2, 10)), [1]);

    testEq('subtract AA', OrderedSet.subtract(arr136, arr136), []);
    testEq('subtract AA1', OrderedSet.subtract(arr136, OrderedSet.ofSortedArray([2, 3, 4, 6, 7])), [1]);
    testEq('subtract AA2', OrderedSet.subtract(arr136, OrderedSet.ofSortedArray([0, 1, 6])), [3]);
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