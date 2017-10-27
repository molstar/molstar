/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../collections/iterator'
import IntTuple from '../collections/int-tuple'
import * as Sort from '../collections/sort'
import LinkedIndex from '../collections/linked-index'
import EquivalenceClasses from '../collections/equivalence-classes'
import Interval from '../collections/interval'
import SortedArray from '../collections/sorted-array'

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
        for (let i = 0; i < 10; i++) {
            for (let j = -10; j < 5; j++) {
                const t = IntTuple.create(i, j);
                expect(IntTuple.fst(t)).toBe(i);
                expect(IntTuple.snd(t)).toBe(j);
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

describe('interval', () => {
    function testI(name: string, a: Interval, b: Interval) {
        it(name, () => expect(Interval.areEqual(a, b)).toBe(true));
    }

    function test(name: string, a: any, b: any) {
        it(name, () => expect(a).toEqual(b));
    }

    const e = Interval.Empty;
    const r05 = Interval.ofRange(0, 5);
    const se05 = Interval.ofBounds(0, 5);

    test('size', Interval.size(e), 0);
    test('size', Interval.size(r05), 6);
    test('size', Interval.size(se05), 5);

    test('min/max', [Interval.min(e), Interval.max(e)], [0, -1]);
    test('min/max', [Interval.min(r05), Interval.max(r05)], [0, 5]);
    test('min/max', [Interval.min(se05), Interval.max(se05)], [0, 4]);

    test('start/end', [Interval.start(e), Interval.end(e)], [0, 0]);
    test('start/end', [Interval.start(r05), Interval.end(r05)], [0, 6]);
    test('start/end', [Interval.start(se05), Interval.end(se05)], [0, 5]);

    test('has', Interval.has(e, 5), false);
    test('has', Interval.has(r05, 5), true);
    test('has', Interval.has(r05, 6), false);
    test('has', Interval.has(r05, -1), false);
    test('has', Interval.has(se05, 5), false);
    test('has', Interval.has(se05, 4), true);

    test('indexOf', Interval.indexOf(e, 5), -1);
    test('indexOf', Interval.indexOf(r05, 5), 5);
    test('indexOf', Interval.indexOf(r05, 6), -1);

    test('getAt', Interval.getAt(r05, 5), 5);

    test('areEqual', Interval.areEqual(r05, se05), false);
    test('areIntersecting1', Interval.areIntersecting(r05, se05), true);
    test('areIntersecting2', Interval.areIntersecting(r05, e), false);
    test('areIntersecting3', Interval.areIntersecting(e, r05), false);
    test('areIntersecting4', Interval.areIntersecting(e, e), true);

    test('areIntersecting5', Interval.areIntersecting(Interval.ofRange(0, 5), Interval.ofRange(-4, 3)), true);
    test('areIntersecting6', Interval.areIntersecting(Interval.ofRange(0, 5), Interval.ofRange(-4, -3)), false);
    test('areIntersecting7', Interval.areIntersecting(Interval.ofRange(0, 5), Interval.ofRange(1, 2)), true);
    test('areIntersecting8', Interval.areIntersecting(Interval.ofRange(0, 5), Interval.ofRange(3, 6)), true);

    test('isSubInterval', Interval.isSubInterval(Interval.ofRange(0, 5), Interval.ofRange(3, 6)), false);
    test('isSubInterval', Interval.isSubInterval(Interval.ofRange(0, 5), Interval.ofRange(3, 5)), true);

    testI('intersect', Interval.intersect(Interval.ofRange(0, 5), Interval.ofRange(-4, 3)), Interval.ofRange(0, 3));
    testI('intersect1', Interval.intersect(Interval.ofRange(0, 5), Interval.ofRange(1, 3)), Interval.ofRange(1, 3));
    testI('intersect2', Interval.intersect(Interval.ofRange(0, 5), Interval.ofRange(3, 5)), Interval.ofRange(3, 5));
    testI('intersect3', Interval.intersect(Interval.ofRange(0, 5), Interval.ofRange(-4, -3)), Interval.Empty);

    test('predIndex1', Interval.findPredecessorIndex(r05, 5), 5);
    test('predIndex2', Interval.findPredecessorIndex(r05, -1), 0);
    test('predIndex3', Interval.findPredecessorIndex(r05, 6), 6);
    test('predIndexInt', Interval.findPredecessorIndexInInterval(r05, 0, Interval.ofRange(2, 3)), 2);
    test('predIndexInt1', Interval.findPredecessorIndexInInterval(r05, 4, Interval.ofRange(2, 3)), 4);

    testI('findRange', Interval.findRange(r05, 2, 3), Interval.ofRange(2, 3));
});

describe('sortedArray', () => {
    function testI(name: string, a: Interval, b: Interval) {
        it(name, () => expect(Interval.areEqual(a, b)).toBe(true));
    }

    function test(name: string, a: any, b: any) {
        it(name, () => expect(a).toEqual(b));
    }

    const a1234 = SortedArray.ofSortedArray([1, 2, 3, 4]);
    const a2468 = SortedArray.ofSortedArray([2, 4, 6, 8]);

    test('size', SortedArray.size(a1234), 4);

    test('min/max', [SortedArray.min(a1234), SortedArray.max(a1234)], [1, 4]);
    test('start/end', [SortedArray.start(a1234), SortedArray.end(a1234)], [1, 5]);

    test('has', SortedArray.has(a1234, 5), false);
    test('has', SortedArray.has(a1234, 4), true);

    it('has-all', () => {
        for (let i = 1; i <= 4; i++) expect(SortedArray.has(a1234, i)).toBe(true);
    });

    test('indexOf', SortedArray.indexOf(a2468, 5), -1);
    test('indexOf', SortedArray.indexOf(a2468, 2), 0);

    test('getAt', SortedArray.getAt(a2468, 1), 4);

    test('areEqual', SortedArray.areEqual(a2468, a2468), true);
    test('areEqual1', SortedArray.areEqual(a2468, SortedArray.create([4, 2, 8, 6])), true);
    test('areEqual2', SortedArray.areEqual(a1234, a2468), false);

    test('predIndex1', SortedArray.findPredecessorIndex(a1234, 5), 4);
    test('predIndex2', SortedArray.findPredecessorIndex(a1234, 2), 1);
    test('predIndex3', SortedArray.findPredecessorIndex(a2468, 4), 1);
    test('predIndex4', SortedArray.findPredecessorIndex(a2468, 3), 1);
    test('predIndexInt', SortedArray.findPredecessorIndexInInterval(a1234, 0, Interval.ofRange(2, 3)), 2);

    testI('findRange', SortedArray.findRange(a2468, 2, 4), Interval.ofRange(0, 1));
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

describe('equiv-classes', () => {
    it('integer mod classes', () => {
        const cls = EquivalenceClasses<number, number>(x => x % 2, (a, b) => (a - b) % 2 === 0);
        for (let i = 0; i < 6; i++) cls.add(i, i);

        expect(cls.groups.length).toBe(2);
        expect(cls.groups[0]).toEqual([0, 2, 4]);
        expect(cls.groups[1]).toEqual([1, 3, 5]);
    });
});