/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import {
  arrayMax,
  arrayMean,
  arrayMin,
  arrayRms,
  arrayStats,
  filterInPlace,
  getSharedSerial,
  isSharedSerial,
  range,
} from '../array.js';

describe('filterInPlace', () => {
  it('filterInPlace works', async () => {
    expect(filterInPlace([], () => true)).toEqual([]);
    expect(filterInPlace([], () => false)).toEqual([]);
    expect(filterInPlace([1, 2, 3, 4, 5], () => true)).toEqual([1, 2, 3, 4, 5]);
    expect(filterInPlace([1, 2, 3, 4, 5], () => false)).toEqual([]);
    expect(filterInPlace([1, 2, 3, 4, 5], (x) => x % 2 === 0)).toEqual([2, 4]);
    expect(filterInPlace([1, 2, 3, 4, 5], (x) => x % 2 === 1)).toEqual([1, 3, 5]);
    expect(filterInPlace([1, 2, 3, 4, 5], (x) => x <= 2)).toEqual([1, 2]);
    expect(filterInPlace([1, 2, 3, 4, 5], (x) => x > 2)).toEqual([3, 4, 5]);
  });
  it('filterInPlace works in place', async () => {
    const array = [1, 2, 3, 4, 5];
    filterInPlace(array, (x) => x % 2 === 1);
    expect(array).toEqual([1, 3, 5]);
  });
  it('filterInPlace big data', async () => {
    const array = range(10 ** 5);
    const expectedResult = array.filter((x) => x % 7 === 0);
    expect(filterInPlace(array, (x) => x % 7 === 0)).toEqual(expectedResult);
  });
});

describe('arrayStats', () => {
  it('equals arrayMin, arrayMax, arrayMean and arrayRms bit for bit (perf-audit-2 04/P18)', () => {
    let seed = 7;
    const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647) * 200 - 37;
    for (const a of [
      Float32Array.from({ length: 1001 }, rnd),
      Array.from({ length: 77 }, rnd),
      [3],
      [-1, Number.NaN, 2],
    ]) {
      expect(arrayStats(a)).toEqual({ min: arrayMin(a), max: arrayMax(a), mean: arrayMean(a), rms: arrayRms(a) });
    }
  });
});

describe('getSharedSerial', () => {
  it('is 0 to n - 1, views of one array grown as asked, the views before kept', () => {
    const a = getSharedSerial(5);
    expect(Array.from(a)).toEqual([0, 1, 2, 3, 4]);
    const b = getSharedSerial(1000);
    expect(b.length).toBe(1000);
    for (let i = 0; i < 1000; ++i) expect(b[i]).toBe(i);
    expect(Array.from(a)).toEqual([0, 1, 2, 3, 4]);
    // a view of the same array where it is long enough
    expect(getSharedSerial(10).buffer).toBe(b.buffer);
  });

  it('tells its views from other arrays', () => {
    const a = getSharedSerial(8);
    expect(isSharedSerial(a)).toBe(true);
    expect(isSharedSerial(getSharedSerial(0))).toBe(true);
    // not from the start, a copy, another type
    expect(isSharedSerial(a.subarray(2))).toBe(false);
    expect(isSharedSerial(a.slice())).toBe(false);
    expect(isSharedSerial(Uint32Array.from([0, 1, 2]))).toBe(false);
    expect(isSharedSerial([0, 1, 2])).toBe(false);
  });
});
