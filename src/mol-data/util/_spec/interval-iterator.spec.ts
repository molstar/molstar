/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Interval, OrderedSet, SortedArray } from '../../int';
import { IntervalIterator } from '../interval-iterator';

describe('interval', () => {
    function testIterator(name: string, interval: Interval, set: OrderedSet, expectedValues: { index: number[], start: number[], end: number[]}) {
        it(`iterator, ${name}`, () => {
            const intervalIt = new IntervalIterator(interval, set);
            const { index, start, end } = expectedValues;

            let i = 0;
            while (intervalIt.hasNext) {
                const segment = intervalIt.move();
                expect(segment.index).toBe(index[i]);
                expect(segment.start).toBe(start[i]);
                expect(segment.end).toBe(end[i]);
                ++i;
            }
            expect(i).toBe(index.length);
        });
    }

    testIterator('basic',
        Interval.ofRange(0, 5),
        SortedArray.ofSortedArray([1, 3, 7, 8]),
        { index: [1, 3], start: [0, 1], end: [1, 2] }
    );
});