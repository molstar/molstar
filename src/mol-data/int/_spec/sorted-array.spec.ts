/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Interval from '../interval'
import SortedArray from '../sorted-array'

describe('sortedArray', () => {
    function testI(name: string, a: Interval, b: Interval) {
        it(name, () => expect(Interval.areEqual(a, b)).toBe(true));
    }

    function test(name: string, a: any, b: any) {
        it(name, () => expect(a).toEqual(b));
    }

    function compareArrays(a: ArrayLike<number>, b: ArrayLike<number>) {
        expect(a.length).toBe(b.length);
        for (let i = 0; i < a.length; i++) expect(a[i]).toBe(b[i]);
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

    test('getAt', a2468[1], 4);

    test('areEqual', SortedArray.areEqual(a2468, a2468), true);
    test('areEqual1', SortedArray.areEqual(a2468, SortedArray.ofUnsortedArray([4, 2, 8, 6])), true);
    test('areEqual2', SortedArray.areEqual(a1234, a2468), false);

    test('predIndex1', SortedArray.findPredecessorIndex(a1234, 5), 4);
    test('predIndex2', SortedArray.findPredecessorIndex(a1234, 2), 1);
    test('predIndex3', SortedArray.findPredecessorIndex(a2468, 4), 1);
    test('predIndex4', SortedArray.findPredecessorIndex(a2468, 3), 1);
    test('predIndexInt', SortedArray.findPredecessorIndexInInterval(a1234, 0, Interval.ofRange(2, 3)), 2);

    testI('findRange', SortedArray.findRange(a2468, 2, 4), Interval.ofRange(0, 1));

    it('deduplicate', () => {
        compareArrays(SortedArray.deduplicate(SortedArray.ofSortedArray([1, 1, 1, 1])), [1]);
        compareArrays(SortedArray.deduplicate(SortedArray.ofSortedArray([1, 1, 2, 2, 3, 4])), [1, 2, 3, 4]);
        compareArrays(SortedArray.deduplicate(SortedArray.ofSortedArray([1, 2, 3])), [1, 2, 3]);
    });

    it('indicesOf', () => {
        compareArrays(SortedArray.indicesOf(SortedArray.ofSortedArray([10, 11, 12]), SortedArray.ofSortedArray([10, 12, 14])), [0, 2]);
    })

    test('intersectionSize', SortedArray.intersectionSize(a1234, a2468), 2);

    it('union1', () => {
        compareArrays(
            SortedArray.union(SortedArray.ofSortedArray([830, 831, 832, 833, 834, 836, 837, 838, 839, 840, 841, 842, 843]), SortedArray.ofSortedArray([835])),
            SortedArray.ofSortedArray([830, 831, 832, 833, 834, 835, 836, 837, 838, 839, 840, 841, 842, 843])
        )
    })

    it('union2', () => {
        compareArrays(
            SortedArray.union(SortedArray.ofSortedArray([830, 832, 833]), SortedArray.ofSortedArray([831])),
            SortedArray.ofSortedArray([830, 831, 832, 833])
        )
    })

    // console.log(Interval.findPredecessorIndexInInterval(Interval.ofBounds(0, 3), 2, Interval.ofBounds(0, 3)))
    // console.log(SortedArray.findPredecessorIndexInInterval(SortedArray.ofSortedArray([0, 1, 2]), 2, Interval.ofBounds(0, 3)))
});