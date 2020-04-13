/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Interval from '../interval';

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

    test('predIndexInt2', Interval.findPredecessorIndex(Interval.ofRange(3, 10), 5), 2);
    test('predIndexInt3', Interval.findPredecessorIndexInInterval(Interval.ofRange(3, 10), 5, Interval.ofRange(2, 6)), 2);

    testI('findRange', Interval.findRange(r05, 2, 3), Interval.ofRange(2, 3));

    test('intersectionSize1', Interval.intersectionSize(Interval.ofRange(0, 5), Interval.ofRange(0, 5)), 6);
    test('intersectionSize2', Interval.intersectionSize(Interval.ofRange(0, 5), Interval.ofRange(1, 2)), 2);
    test('intersectionSize3', Interval.intersectionSize(Interval.ofRange(1, 2), Interval.ofRange(0, 5)), 2);
    test('intersectionSize4', Interval.intersectionSize(Interval.ofRange(0, 5), Interval.ofRange(3, 8)), 3);
    test('intersectionSize5', Interval.intersectionSize(Interval.ofRange(0, 5), Interval.ofRange(6, 8)), 0);
});