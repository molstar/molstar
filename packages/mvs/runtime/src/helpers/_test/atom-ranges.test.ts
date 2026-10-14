/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ElementRanges, ElementRangesBuffer } from '../element-ranges.js';
import { IntBuffer } from '../int-buffer.js';
import { Sorting } from '../indexing.js';
import { SortedArray } from '@molstar/core/data/int';

describe('union', () => {
  it('union non-overlapping', async () => {
    const a = {
      from: [0, 20, 40, 60, 80],
      to: [10, 30, 50, 70, 90],
    } as ElementRanges;
    const b = {
      from: [11, 37, 51, 205],
      to: [15, 39, 55, 210],
    } as ElementRanges;
    const c = {
      from: [-10, 200, 300],
      to: [-5, 202, 305],
    } as ElementRanges;
    const result = {
      from: [-10, 0, 11, 20, 37, 40, 51, 60, 80, 200, 205, 300],
      to: [-5, 10, 15, 30, 39, 50, 55, 70, 90, 202, 210, 305],
    } as ElementRanges;
    expect(ElementRanges.union([a, b, c])).toEqual(result);
  });
  it('union overlapping', async () => {
    const a = {
      from: [0, 20, 40, 60, 80],
      to: [10, 30, 50, 70, 90],
    } as ElementRanges;
    const b = {
      from: [10, 37, 51, 84, 205],
      to: [15, 40, 55, 88, 220],
    } as ElementRanges;
    const c = {
      from: [-10, 67, 200, 300],
      to: [5, 80, 210, 305],
    } as ElementRanges;
    const result = {
      from: [-10, 20, 37, 51, 60, 200, 300],
      to: [15, 30, 50, 55, 90, 220, 305],
    } as ElementRanges;
    expect(ElementRanges.union([a, b, c])).toEqual(result);
  });
});

describe('union (sweep)', () => {
  const empty = { from: [], to: [] } as ElementRanges;
  it('handles empty and undefined inputs', () => {
    expect(ElementRanges.union([])).toEqual(empty);
    expect(ElementRanges.union([undefined, empty])).toEqual(empty);
  });
  it('merges touching and nested ranges', () => {
    const a = { from: [0, 10], to: [10, 20] } as ElementRanges; // touching
    const b = { from: [5, 12], to: [8, 14] } as ElementRanges; // nested in a
    const c = { from: [30], to: [31] } as ElementRanges;
    expect(ElementRanges.union([c, b, a])).toEqual({ from: [0, 30], to: [20, 31] });
  });
  it('matches brute force on random input', () => {
    let seed = 1;
    const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    for (let iter = 0; iter < 200; iter++) {
      const parts: ElementRanges[] = [];
      const covered = new Set<number>();
      for (let p = 0; p < 1 + Math.floor(rnd() * 4); p++) {
        const part = { from: [], to: [] } as unknown as ElementRanges;
        let pos = Math.floor(rnd() * 20);
        for (let r = 0; r < Math.floor(rnd() * 8); r++) {
          const from = pos;
          const to = from + 1 + Math.floor(rnd() * 10);
          part.from.push(from as any);
          part.to.push(to as any);
          for (let x = from; x < to; x++) covered.add(x);
          pos = to + Math.floor(rnd() * 10);
        }
        parts.push(part);
      }
      const result = ElementRanges.union(parts);
      const got = new Set<number>();
      for (let i = 0; i < result.from.length; i++) {
        if (i > 0) expect(result.from[i]).toBeGreaterThan(result.to[i - 1]); // disjoint, sorted, not touching
        for (let x = result.from[i]; x < result.to[i]; x++) got.add(x);
      }
      expect([...got].sort((a, b) => a - b)).toEqual([...covered].sort((a, b) => a - b));
    }
  });
});

describe('ElementRangesBuffer', () => {
  it('merges touching ranges, rejects overlaps, grows, and clears without losing capacity', () => {
    const b = new ElementRangesBuffer();
    b.add(0, 5).add(5, 8).add(10, 12);
    expect(b.toRanges()).toEqual({ from: [0, 10], to: [8, 12] });
    expect(() => b.add(11, 20)).toThrow();
    for (let i = 0; i < 1000; i++) b.add(100 + 2 * i, 101 + 2 * i);
    expect(b.count).toBe(1002);
    const capacity = b.from.length;
    b.clear();
    expect(b.count).toBe(0);
    expect(b.from.length).toBe(capacity);
    expect(b.toRanges()).toEqual({ from: [], to: [] });
  });
  it('append keeps order of insertion, normalize sorts and merges it', () => {
    const b = new ElementRangesBuffer();
    b.append({ from: [20, 0], to: [30, 10] } as ElementRanges).append({ from: [10], to: [15] } as ElementRanges);
    expect(ElementRanges.normalize(b)).toEqual({ from: [0, 20], to: [15, 30] });
  });
});

describe('IntBuffer', () => {
  it('grows, keeps contents, and views only the filled part', () => {
    const b = new IntBuffer(2);
    for (let i = 0; i < 10; i++) b.push(i);
    expect(Array.from(b.view())).toEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    b.clear();
    b.push(42);
    expect(Array.from(b.view())).toEqual([42]);
  });
});

describe('Sorting.findIndexRange', () => {
  const keys = [10, 11, 12, 13, 14, 15];
  const valueList = [1, 3, 3, 5, 7, 7];
  const values = SortedArray.ofSortedArray(valueList);
  const sorting = { keys, values } as Sorting<number, number>;
  const find = (min?: any, max?: any) => {
    const r = Sorting.findIndexRange(sorting, min, max, { from: 0, to: 0 });
    return keys.slice(r.from, r.to);
  };
  /** The definition: keys whose value is in `[min, max]` (undefined = unbounded) */
  const bruteForce = (min?: number, max?: number) =>
    keys.filter((_, i) => (min === undefined || valueList[i] >= min) && (max === undefined || valueList[i] <= max));

  it('finds closed intervals, open ends, empty and inverted intervals', () => {
    expect(find(3, 3)).toEqual([11, 12]);
    expect(find(2, 6)).toEqual([11, 12, 13]);
    expect(find(undefined, 3)).toEqual([10, 11, 12]);
    expect(find(5, undefined)).toEqual([13, 14, 15]);
    expect(find(undefined, undefined)).toEqual(keys);
    expect(find(4, 4)).toEqual([]);
    expect(find(7, 3)).toEqual([]);
    expect(find(100, 200)).toEqual([]);
  });
  it('agrees with the definition for every pair of bounds, integer or not', () => {
    const bounds = [undefined, -1, 0, 0.5, 1, 2, 2.5, 3, 4, 5, 6.99, 7, 7.01, 8, 100];
    for (const min of bounds) for (const max of bounds) expect(find(min, max)).toEqual(bruteForce(min, max));
  });
  it('compares bounds as numbers: numeric strings work as their numbers, other strings select nothing', () => {
    expect(find('3', '3')).toEqual([11, 12]);
    expect(find('5', undefined)).toEqual([13, 14, 15]);
    expect(find(undefined, '5')).toEqual([10, 11, 12, 13]);
    expect(find('2', '6')).toEqual([11, 12, 13]); // (`'6' + 1` would be '61')
    expect(find('x', 5)).toEqual([]);
    expect(find(1, NaN)).toEqual([]);
  });
  it('getKeysWithValueInRange gives the keys of the same range as a new array', () => {
    for (const [min, max] of [
      [3, 5],
      [undefined, 3],
      [5, undefined],
      [4, 4],
      [7, 3],
    ] as const) {
      const result = Sorting.getKeysWithValueInRange(sorting, min, max);
      expect(result).toEqual(bruteForce(min, max));
      expect(result).not.toBe(keys); // (newly allocated)
    }
  });
});
