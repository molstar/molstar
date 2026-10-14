/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { NumberMap } from '../number-map.js';

describe('NumberMap with an offset and strict keys', () => {
  it('keeps the keys from the offset in the array, others in the map', () => {
    const map = new NumberMap<number, number>(5, 11); // 11..15
    for (let k = 11; k <= 15; k++) map.set(k, k - 11);
    expect([11, 13, 15].map((k) => map.get(k))).toEqual([0, 2, 4]);
    expect(map.overflowSize).toBe(0);
    map.set(10, 7);
    map.set(16, 8);
    map.set(1000, 9);
    expect([10, 16, 1000, 17].map((k) => map.get(k))).toEqual([7, 8, 9, undefined]);
    expect(map.overflowSize).toBe(3);
  });
  it('a key which is not a number has no value (an array would convert it)', () => {
    const map = new NumberMap<number, number>(5, 0);
    map.set(3, 1);
    expect(map.get('3' as any)).toBeUndefined();
    expect(map.get(3)).toBe(1);
  });
});

describe('NumberMap', () => {
  it('keeps non-negative integer values for keys in and out of the limit, and has no value for the others', () => {
    const map = new NumberMap<number, number>(5);
    expect(map.get(0)).toBeUndefined();
    map.set(0, 0); // (0 is a value, not a missing one)
    map.set(3, 7);
    map.set(10, 2); // (out of the limit)
    map.set(-1, 4);
    expect([map.get(0), map.get(3), map.get(4), map.get(10), map.get(-1), map.get(11)]).toEqual([
      0,
      7,
      undefined,
      2,
      4,
      undefined,
    ]);
    map.set(3, 9);
    expect(map.get(3)).toBe(9);
  });
});
