/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../collections/iterator'
import IntPair from '../collections/int-pair'
import * as Sort from '../collections/sort'
import OrderedSet from '../collections/ordered-set'
import OrderedSet1 from '../collections/ordered-set.1'
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
                const t = IntPair.pack1(i, j);
                IntPair.unpack(t, p);
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

function ordSetToArray(set: OrderedSet) {
    const ret = [];
    for (let i = 0; i < set.size; i++) ret.push(set.elementAt(i));
    return ret;
}

describe('ordered set', () => {
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
        expect(empty.has(10)).toBe(false);
        expect(empty.indexOf(10)).toBe(-1);

        expect(singleton10.has(10)).toBe(true);
        expect(singleton10.has(11)).toBe(false);
        expect(singleton10.indexOf(10)).toBe(0);
        expect(singleton10.indexOf(11)).toBe(-1);

        expect(range1_4.has(4)).toBe(true);
        expect(range1_4.has(5)).toBe(false);
        expect(range1_4.indexOf(4)).toBe(3);
        expect(range1_4.indexOf(11)).toBe(-1);

        expect(arr136.has(3)).toBe(true);
        expect(arr136.has(4)).toBe(false);
        expect(arr136.indexOf(3)).toBe(1);
        expect(arr136.indexOf(11)).toBe(-1);
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

describe('ordered set 1', () => {
    function ordSet1ToArray(set: OrderedSet1) {
        const ret = [];
        for (let i = 0, _i = OrderedSet1.size(set); i < _i; i++) ret.push(OrderedSet1.elementAt(set, i));
        return ret;
    }

    function testEq(name: string, set: OrderedSet1, expected: number[]) {
        it(name, () => {
            // copy the arrays to ensure "compatibility" between typed and native arrays
            expect(Array.prototype.slice.call(ordSet1ToArray(set))).toEqual(Array.prototype.slice.call(expected));
        });
    }

    const empty = OrderedSet1.Empty;
    const singleton10 = OrderedSet1.ofSingleton(10);
    const range1_4 = OrderedSet1.ofRange(1, 4);
    const arr136 = OrderedSet1.ofSortedArray([1, 3, 6]);

    testEq('empty', empty, []);
    testEq('singleton', singleton10, [10]);
    testEq('range', range1_4, [1, 2, 3, 4]);
    testEq('sorted array', arr136, [1, 3, 6]);

    it('equality', () => {
        expect(OrderedSet1.areEqual(empty, singleton10)).toBe(false);
        expect(OrderedSet1.areEqual(singleton10, singleton10)).toBe(true);
        expect(OrderedSet1.areEqual(range1_4, singleton10)).toBe(false);
        expect(OrderedSet1.areEqual(arr136, OrderedSet1.ofSortedArray([1, 3, 6]))).toBe(true);
        expect(OrderedSet1.areEqual(arr136, OrderedSet1.ofSortedArray([1, 4, 6]))).toBe(false);
    });

    it('areIntersecting', () => {
        expect(OrderedSet1.areIntersecting(range1_4, arr136)).toBe(true);
        expect(OrderedSet1.areIntersecting(empty, empty)).toBe(true);
        expect(OrderedSet1.areIntersecting(empty, singleton10)).toBe(false);
        expect(OrderedSet1.areIntersecting(empty, range1_4)).toBe(false);
        expect(OrderedSet1.areIntersecting(empty, arr136)).toBe(false);
    });

    it('isSubset', () => {
        expect(OrderedSet1.isSubset(singleton10, empty)).toBe(true);
        expect(OrderedSet1.isSubset(range1_4, empty)).toBe(true);
        expect(OrderedSet1.isSubset(arr136, empty)).toBe(true);
        expect(OrderedSet1.isSubset(empty, empty)).toBe(true);
        expect(OrderedSet1.isSubset(empty, singleton10)).toBe(false);
        expect(OrderedSet1.isSubset(empty, range1_4)).toBe(false);
        expect(OrderedSet1.isSubset(empty, arr136)).toBe(false);

        expect(OrderedSet1.isSubset(singleton10, range1_4)).toBe(false);
        expect(OrderedSet1.isSubset(range1_4, OrderedSet1.ofRange(2, 3))).toBe(true);
        expect(OrderedSet1.isSubset(arr136, range1_4)).toBe(false);
        expect(OrderedSet1.isSubset(arr136, arr136)).toBe(true);
        expect(OrderedSet1.isSubset(arr136, OrderedSet1.ofSortedArray([1, 3]))).toBe(true);
        expect(OrderedSet1.isSubset(arr136, OrderedSet1.ofSortedArray([1, 3, 7]))).toBe(false);
        expect(OrderedSet1.isSubset(OrderedSet1.ofSortedArray([0, 1, 2, 3, 7, 10]), OrderedSet1.ofSortedArray([1, 3, 7]))).toBe(true);
        expect(OrderedSet1.isSubset(arr136, OrderedSet1.ofSortedArray([1, 3, 10, 45]))).toBe(false);
    });

    it('access/membership', () => {
        expect(OrderedSet1.has(empty, 10)).toBe(false);
        expect(OrderedSet1.indexOf(empty, 10)).toBe(-1);

        expect(OrderedSet1.has(singleton10, 10)).toBe(true);
        expect(OrderedSet1.has(singleton10, 11)).toBe(false);
        expect(OrderedSet1.indexOf(singleton10, 10)).toBe(0);
        expect(OrderedSet1.indexOf(singleton10, 11)).toBe(-1);

        expect(OrderedSet1.has(range1_4, 4)).toBe(true);
        expect(OrderedSet1.has(range1_4, 5)).toBe(false);
        expect(OrderedSet1.indexOf(range1_4, 4)).toBe(3);
        expect(OrderedSet1.indexOf(range1_4, 11)).toBe(-1);

        expect(OrderedSet1.has(arr136, 3)).toBe(true);
        expect(OrderedSet1.has(arr136, 4)).toBe(false);
        expect(OrderedSet1.indexOf(arr136, 3)).toBe(1);
        expect(OrderedSet1.indexOf(arr136, 11)).toBe(-1);
    });

    it('interval range', () => {
        expect(OrderedSet1.getIntervalRange(empty, 9, 11)).toEqual({ start: 0, end: 0 });
        expect(OrderedSet1.getIntervalRange(empty, -9, -6)).toEqual({ start: 0, end: 0 });
        expect(OrderedSet1.getIntervalRange(singleton10, 9, 11)).toEqual({ start: 0, end: 1 });
        expect(OrderedSet1.getIntervalRange(range1_4, 2, 3)).toEqual({ start: 1, end: 3 });
        expect(OrderedSet1.getIntervalRange(range1_4, -10, 2)).toEqual({ start: 0, end: 2 });
        expect(OrderedSet1.getIntervalRange(range1_4, -10, 20)).toEqual({ start: 0, end: 4 });
        expect(OrderedSet1.getIntervalRange(range1_4, 3, 20)).toEqual({ start: 2, end: 4 });
        expect(OrderedSet1.getIntervalRange(arr136, 0, 1)).toEqual({ start: 0, end: 1 });
        expect(OrderedSet1.getIntervalRange(arr136, 0, 3)).toEqual({ start: 0, end: 2 });
        expect(OrderedSet1.getIntervalRange(arr136, 0, 4)).toEqual({ start: 0, end: 2 });
        expect(OrderedSet1.getIntervalRange(arr136, 2, 4)).toEqual({ start: 1, end: 2 });
        expect(OrderedSet1.getIntervalRange(arr136, 2, 7)).toEqual({ start: 1, end: 3 });
    })

    testEq('union ES', OrderedSet1.union(empty, singleton10), [10]);
    testEq('union ER', OrderedSet1.union(empty, range1_4), [1, 2, 3, 4]);
    testEq('union EA', OrderedSet1.union(empty, arr136), [1, 3, 6]);
    testEq('union SS', OrderedSet1.union(singleton10, OrderedSet1.ofSingleton(16)), [10, 16]);
    testEq('union SR', OrderedSet1.union(range1_4, singleton10), [1, 2, 3, 4, 10]);
    testEq('union SA', OrderedSet1.union(arr136, singleton10), [1, 3, 6, 10]);
    testEq('union SA1', OrderedSet1.union(arr136, OrderedSet1.ofSingleton(3)), [1, 3, 6]);
    testEq('union RR', OrderedSet1.union(range1_4, range1_4), [1, 2, 3, 4]);
    testEq('union RR1', OrderedSet1.union(range1_4, OrderedSet1.ofRange(6, 7)), [1, 2, 3, 4, 6, 7]);
    testEq('union RR2', OrderedSet1.union(range1_4, OrderedSet1.ofRange(3, 5)), [1, 2, 3, 4, 5]);
    testEq('union RA', OrderedSet1.union(range1_4, arr136), [1, 2, 3, 4, 6]);
    testEq('union AA', OrderedSet1.union(arr136, OrderedSet1.ofSortedArray([2, 4, 6, 7])), [1, 2, 3, 4, 6, 7]);
    testEq('union AA1', OrderedSet1.union(arr136, OrderedSet1.ofSortedArray([2, 3, 4, 6, 7])), [1, 2, 3, 4, 6, 7]);
    testEq('union AA2', OrderedSet1.union(arr136, OrderedSet1.ofSortedArray([2, 4, 5, 6, 7])), [1, 2, 3, 4, 5, 6, 7]);

    testEq('intersect ES', OrderedSet1.intersect(empty, singleton10), []);
    testEq('intersect ER', OrderedSet1.intersect(empty, range1_4), []);
    testEq('intersect EA', OrderedSet1.intersect(empty, arr136), []);
    testEq('intersect SS', OrderedSet1.intersect(singleton10, OrderedSet1.ofSingleton(16)), []);
    testEq('intersect SS1', OrderedSet1.intersect(singleton10, singleton10), [10]);
    testEq('intersect SR', OrderedSet1.intersect(range1_4, singleton10), []);
    testEq('intersect RR', OrderedSet1.intersect(range1_4, range1_4), [1, 2, 3, 4]);
    testEq('intersect RR2', OrderedSet1.intersect(range1_4, OrderedSet1.ofRange(3, 5)), [3, 4]);
    testEq('intersect RA', OrderedSet1.intersect(range1_4, arr136), [1, 3]);
    testEq('intersect AA', OrderedSet1.intersect(arr136, OrderedSet1.ofSortedArray([2, 3, 4, 6, 7])), [3, 6]);

    testEq('subtract ES', OrderedSet1.subtract(empty, singleton10), []);
    testEq('subtract ER', OrderedSet1.subtract(empty, range1_4), []);
    testEq('subtract EA', OrderedSet1.subtract(empty, arr136), []);
    testEq('subtract SS', OrderedSet1.subtract(singleton10, OrderedSet1.ofSingleton(16)), [10]);
    testEq('subtract SS1', OrderedSet1.subtract(singleton10, singleton10), []);
    testEq('subtract SR', OrderedSet1.subtract(range1_4, singleton10), [1, 2, 3, 4]);
    testEq('subtract SR1', OrderedSet1.subtract(range1_4, OrderedSet1.ofSingleton(4)), [1, 2, 3]);
    testEq('subtract SR2', OrderedSet1.subtract(range1_4, OrderedSet1.ofSingleton(3)), [1, 2, 4]);
    testEq('subtract RR', OrderedSet1.subtract(range1_4, range1_4), []);
    testEq('subtract RR1', OrderedSet1.subtract(range1_4, OrderedSet1.ofRange(3, 5)), [1, 2]);

    testEq('subtract RA', OrderedSet1.subtract(range1_4, arr136), [2, 4]);
    testEq('subtract RA1', OrderedSet1.subtract(range1_4, OrderedSet1.ofSortedArray([0, 1, 2, 3, 4, 7])), []);
    testEq('subtract RA2', OrderedSet1.subtract(range1_4, OrderedSet1.ofSortedArray([0, 2, 3])), [1, 4]);

    testEq('subtract AR', OrderedSet1.subtract(arr136, range1_4), [6]);
    testEq('subtract AR1', OrderedSet1.subtract(arr136, OrderedSet1.ofRange(0, 10)), []);
    testEq('subtract AR1', OrderedSet1.subtract(arr136, OrderedSet1.ofRange(2, 10)), [1]);

    testEq('subtract AA', OrderedSet1.subtract(arr136, arr136), []);
    testEq('subtract AA1', OrderedSet1.subtract(arr136, OrderedSet1.ofSortedArray([2, 3, 4, 6, 7])), [1]);
    testEq('subtract AA2', OrderedSet1.subtract(arr136, OrderedSet1.ofSortedArray([0, 1, 6])), [3]);
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
    const r = (i: number, j: number) => IntPair.pack1(i, j);

    function setToPairs(set: MultiSet): IntPair[] {
        const ret = [];
        const it = MultiSet.values(set);
        for (let v = it.move(); !it.done; v = it.move()) ret[ret.length] = IntPair.create(v.fst, v.snd);
        return ret;
    }

    it('singleton pair', () => {
        const set = MultiSet.create(p(10, 11));
        expect(setToPairs(set)).toEqual([p(10, 11)]);
    });

    it('singleton number', () => {
        const set = MultiSet.create(r(10, 11));
        expect(setToPairs(set)).toEqual([p(10, 11)]);
    });

    it('multi', () => {
        const set = MultiSet.create({
            1: OrderedSet.ofSortedArray([4, 6, 7]),
            3: OrderedSet.ofRange(0, 1),
        });
        expect(setToPairs(set)).toEqual([p(1, 4), p(1, 6), p(1, 7), p(3, 0), p(3, 1)]);
    });

    it('packed pairs', () => {
        const set = MultiSet.create([r(1, 3), r(0, 1), r(0, 6), r(0, 2)]);
        expect(setToPairs(set)).toEqual([p(0, 1), p(0, 2), p(0, 6), p(1, 3)]);
    });

    it('equality', () => {
        const a = MultiSet.create([r(1, 3), r(0, 1), r(0, 6), r(0, 2)]);
        const b = MultiSet.create([r(1, 3), r(0, 1), r(0, 6), r(0, 2)]);
        const c = MultiSet.create([r(1, 3), r(0, 4), r(0, 6), r(0, 2)]);
        const d = MultiSet.create([r(1, 3)]);
        const e = MultiSet.create([r(1, 3)]);
        const f = MultiSet.create([r(3, 3)]);

        expect(MultiSet.areEqual(a, a)).toBe(true);
        expect(MultiSet.areEqual(a, b)).toBe(true);
        expect(MultiSet.areEqual(a, c)).toBe(false);
        expect(MultiSet.areEqual(a, d)).toBe(false);
        expect(MultiSet.areEqual(d, d)).toBe(true);
        expect(MultiSet.areEqual(d, e)).toBe(true);
        expect(MultiSet.areEqual(d, f)).toBe(false);
    });

    it('are intersecting', () => {
        const a = MultiSet.create([r(1, 3), r(0, 1), r(0, 6), r(0, 2)]);
        const b = MultiSet.create([r(1, 3), r(0, 1), r(0, 6), r(0, 2)]);
        const c = MultiSet.create([r(1, 3), r(0, 4), r(0, 6), r(0, 2)]);
        const d = MultiSet.create([r(1, 3)]);
        const e = MultiSet.create([r(1, 3)]);
        const f = MultiSet.create([r(3, 3)]);
        const g = MultiSet.create([r(10, 3), r(8, 1), r(7, 6), r(3, 2)]);

        expect(MultiSet.areIntersecting(a, a)).toBe(true);
        expect(MultiSet.areIntersecting(a, b)).toBe(true);
        expect(MultiSet.areIntersecting(a, c)).toBe(true);
        expect(MultiSet.areIntersecting(a, d)).toBe(true);
        expect(MultiSet.areIntersecting(a, g)).toBe(false);
        expect(MultiSet.areIntersecting(d, d)).toBe(true);
        expect(MultiSet.areIntersecting(d, e)).toBe(true);
        expect(MultiSet.areIntersecting(d, f)).toBe(false);
    });

    it('intersection', () => {
        const a = MultiSet.create([r(1, 3), r(0, 1), r(0, 6), r(0, 2)]);
        const b = MultiSet.create([r(10, 3), r(0, 1), r(0, 6), r(4, 2)]);
        const c = MultiSet.create([r(1, 3)]);
        const d = MultiSet.create([r(2, 3)]);
        expect(MultiSet.intersect(a, a)).toBe(a);
        expect(setToPairs(MultiSet.intersect(a, b))).toEqual([p(0, 1), p(0, 6)]);
        expect(setToPairs(MultiSet.intersect(a, c))).toEqual([p(1, 3)]);
        expect(setToPairs(MultiSet.intersect(c, d))).toEqual([]);
    });

    it('subtract', () => {
        const a = MultiSet.create([r(1, 3), r(0, 1), r(0, 6), r(0, 2)]);
        const a1 = MultiSet.create([r(1, 3), r(0, 1), r(0, 6), r(0, 2)]);
        const b = MultiSet.create([r(10, 3), r(0, 1), r(0, 6), r(4, 2)]);
        const c = MultiSet.create([r(1, 3)]);
        const d = MultiSet.create([r(2, 3)]);
        expect(setToPairs(MultiSet.subtract(a, a))).toEqual([]);
        expect(setToPairs(MultiSet.subtract(a, a1))).toEqual([]);
        expect(setToPairs(MultiSet.subtract(a, b))).toEqual([p(0, 2), p(1, 3)]);
        expect(setToPairs(MultiSet.subtract(c, d))).toEqual([p(1, 3)]);
        expect(setToPairs(MultiSet.subtract(a, c))).toEqual([p(0, 1), p(0, 2), p(0, 6)]);
        expect(setToPairs(MultiSet.subtract(c, a))).toEqual([]);
        expect(setToPairs(MultiSet.subtract(d, a))).toEqual([p(2, 3)]);
    });

    it('union', () => {
        const a = MultiSet.create([r(1, 3), r(0, 1)]);
        const a1 = MultiSet.create([r(1, 3), r(0, 1)]);
        const b = MultiSet.create([r(10, 3), r(0, 1)]);
        const c = MultiSet.create([r(1, 3)]);
        const d = MultiSet.create([r(2, 3)]);
        expect(MultiSet.unionMany([a])).toBe(a);
        expect(MultiSet.union(a, a)).toBe(a);
        expect(setToPairs(MultiSet.union(a, a))).toEqual([p(0, 1), p(1, 3)]);
        expect(setToPairs(MultiSet.union(a, a1))).toEqual([p(0, 1), p(1, 3)]);
        expect(setToPairs(MultiSet.union(a, b))).toEqual([p(0, 1), p(1, 3), p(10, 3)]);
        expect(setToPairs(MultiSet.union(c, d))).toEqual([p(1, 3), p(2, 3)]);
        expect(setToPairs(MultiSet.unionMany([a, b, c, d]))).toEqual([p(0, 1), p(1, 3), p(2, 3), p(10, 3)]);
    });
});