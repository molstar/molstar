/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, test, vi } from 'vitest';
import { SortedRanges } from '../sorted-ranges.js';
import { OrderedSet } from '../ordered-set.js';
import { SortedArray } from '../sorted-array.js';

describe('rangesArray', () => {
  function test(name: string, a: any, b: any) {
    it(name, () => expect(a).toEqual(b));
  }

  function testIterator(
    name: string,
    ranges: SortedRanges,
    set: OrderedSet,
    expectedValues: { index: number[]; start: number[]; end: number[] },
  ) {
    it(`iterator, ${name}`, () => {
      const rangesIt = SortedRanges.transientSegments(ranges, set);
      const { index, start, end } = expectedValues;

      let i = 0;
      while (rangesIt.hasNext) {
        const segment = rangesIt.move();
        expect(segment.index).toBe(index[i]);
        expect(segment.start).toBe(start[i]);
        expect(segment.end).toBe(end[i]);
        ++i;
      }
      expect(i).toBe(index.length);
    });
  }

  const a1234 = SortedRanges.ofSortedRanges([1, 2, 3, 4]);
  const a1134 = SortedRanges.ofSortedRanges([1, 1, 3, 4]);

  test('size', SortedRanges.size(a1234), 4);
  test('size', SortedRanges.size(a1134), 3);

  test('min/max', [SortedRanges.min(a1234), SortedRanges.max(a1234)], [1, 4]);
  test('start/end', [SortedRanges.start(a1234), SortedRanges.end(a1234)], [1, 5]);

  testIterator('two ranges', SortedRanges.ofSortedRanges([1, 2, 3, 4]), OrderedSet.ofBounds(1, 5), {
    index: [0, 1],
    start: [0, 2],
    end: [2, 4],
  });
  testIterator('first range', SortedRanges.ofSortedRanges([1, 2, 3, 4]), SortedArray.ofSortedArray([2]), {
    index: [0],
    start: [0],
    end: [1],
  });
  testIterator('second range', SortedRanges.ofSortedRanges([1, 2, 3, 4]), SortedArray.ofSortedArray([4]), {
    index: [1],
    start: [0],
    end: [1],
  });
  testIterator('set not in ranges', SortedRanges.ofSortedRanges([1, 2, 3, 4]), SortedArray.ofSortedArray([10]), {
    index: [],
    start: [],
    end: [],
  });
  testIterator(
    'set in second range and beyond',
    SortedRanges.ofSortedRanges([1, 2, 3, 4]),
    SortedArray.ofSortedArray([3, 10]),
    { index: [1], start: [0], end: [1] },
  );
  testIterator('length 1 range', SortedRanges.ofSortedRanges([1, 1, 3, 4]), SortedArray.ofSortedArray([0, 1, 10]), {
    index: [0],
    start: [1],
    end: [2],
  });
});

describe('firstIntersectionIndexFrom (perf-audit-2 01/P7)', () => {
  it('finds the first intersecting range by search, as a scan does', () => {
    let seed = 5;
    const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
    for (let n = 0; n < 100; n++) {
      const bounds: number[] = [];
      for (let v = Math.floor(rnd() * 5), k = 0; k < 12; k++) {
        const len = Math.floor(rnd() * 4);
        bounds.push(v, v + len);
        v += len + 1 + Math.floor(rnd() * 6);
      }
      const ranges = SortedRanges.ofSortedRanges(bounds);
      const xs = Array.from(new Set(Array.from({ length: 6 }, () => Math.floor(rnd() * 80)))).sort((a, b) => a - b);
      const set = OrderedSet.ofSortedArray(xs);
      for (let from = 0; from < SortedRanges.count(ranges); from++) {
        let expected = -1;
        for (let i = from; i < SortedRanges.count(ranges) && expected < 0; i++) {
          if (xs.some((x) => x >= bounds[2 * i] && x <= bounds[2 * i + 1])) expected = i;
        }
        expect(SortedRanges.firstIntersectionIndexFrom(ranges, set, from)).toBe(expected);
      }
    }
  });

  it('does not test the ranges before the set one by one', () => {
    const bounds: number[] = [];
    for (let i = 0; i < 1000; i++) bounds.push(10 * i, 10 * i + 4);
    const ranges = SortedRanges.ofSortedRanges(bounds);
    const intersect = vi.spyOn(OrderedSet, 'areIntersecting');
    const find = vi.spyOn(OrderedSet, 'findPredecessorIndex');
    expect(SortedRanges.firstIntersectionIndex(ranges, OrderedSet.ofSingleton(9002))).toBe(900);
    expect(intersect.mock.calls.length + find.mock.calls.length).toBeLessThan(5);
    intersect.mockRestore();
    find.mockRestore();
  });
});
