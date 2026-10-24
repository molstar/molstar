/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Interval } from '../../int/interval.js';
import { OrderedSet } from '../../int/ordered-set.js';
import { SortedArray } from '../../int/sorted-array.js';
import { IntervalIterator } from '../interval-iterator.js';

describe('interval', () => {
  function testIterator(
    name: string,
    interval: Interval,
    set: OrderedSet,
    expectedValues: { index: number[]; start: number[]; end: number[] },
  ) {
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

  testIterator('basic', Interval.ofRange(0, 5), SortedArray.ofSortedArray([1, 3, 7, 8]), {
    index: [1, 3],
    start: [0, 1],
    end: [1, 2],
  });

  // perf-audit-2 01/B7: an interval not starting at 0 tested its indices for membership, not its values
  testIterator('offset', Interval.ofRange(10, 15), SortedArray.ofSortedArray([10, 12, 13, 20]), {
    index: [0, 2, 3],
    start: [0, 1, 2],
    end: [1, 2, 3],
  });
  testIterator('first value not in the set', Interval.ofRange(5, 9), SortedArray.ofSortedArray([2, 6, 9]), {
    index: [1, 4],
    start: [1, 2],
    end: [2, 3],
  });
  testIterator('disjoint', Interval.ofRange(5, 9), SortedArray.ofSortedArray([1, 2, 12]), {
    index: [],
    start: [],
    end: [],
  });
});
