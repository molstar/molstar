/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { countingSortIndices, radixSort, radixSortIndices } from '../sort.js';

function randomUints(n: number, max: number) {
  const data = new Uint32Array(n);
  for (let i = 0; i < n; i++) data[i] = Math.floor(Math.random() * max);
  return data;
}

function expectSorted(data: ArrayLike<number>, reference: number[]) {
  expect(Array.from(data)).toEqual(reference.slice().sort((a, b) => a - b));
}

describe('radixSortUint', () => {
  it('sorts random values like a comparison sort', () => {
    const data = randomUints(1000, 2 ** 24);
    const reference = Array.from(data);
    expect(radixSort(data)).toBe(data);
    expectSorted(data, reference);
  });

  it('handles empty and single-element arrays', () => {
    expect(Array.from(radixSort(new Int32Array(0)))).toEqual([]);
    expect(Array.from(radixSort(Int32Array.of(42)))).toEqual([42]);
  });

  it('handles duplicates and zeros', () => {
    const data = Int32Array.of(5, 0, 3, 5, 0, 3, 3, 0);
    radixSort(data);
    expect(Array.from(data)).toEqual([0, 0, 0, 3, 3, 3, 5, 5]);
  });

  it('handles already sorted and reversed input', () => {
    const asc = Int32Array.of(1, 2, 3, 4, 5);
    expect(Array.from(radixSort(asc))).toEqual([1, 2, 3, 4, 5]);
    const desc = Int32Array.of(5, 4, 3, 2, 1);
    expect(Array.from(radixSort(desc))).toEqual([1, 2, 3, 4, 5]);
  });

  it('supports the full Uint32Array range', () => {
    const data = Uint32Array.of(4294967295, 0, 2147483648, 1, 2147483647);
    radixSort(data);
    expect(Array.from(data)).toEqual([0, 1, 2147483647, 2147483648, 4294967295]);
  });

  it('respects an explicit maxValue', () => {
    const data = randomUints(500, 255);
    const reference = Array.from(data);
    radixSort(data, 255);
    expectSorted(data, reference);
  });

  it('sorts only the given subarray view', () => {
    const data = Int32Array.of(9, 8, 3, 1, 2, 7, 6);
    radixSort(data.subarray(2, 5));
    expect(Array.from(data)).toEqual([9, 8, 1, 2, 3, 7, 6]);
  });

  it('works with an even and odd number of passes', () => {
    // 1 pass (max < 2^8), 2 passes (< 2^16), 3 passes (< 2^24), 4 passes (< 2^32)
    for (const max of [2 ** 8, 2 ** 16, 2 ** 24, 2 ** 32]) {
      const data = randomUints(300, max);
      const reference = Array.from(data);
      radixSort(data);
      expectSorted(data, reference);
    }
  });
});

describe('radixSortIndices', () => {
  /** keys of a fixed sequence (a linear congruential generator): deterministic */
  function keysOf(n: number, max: number) {
    const keys = new Uint32Array(n);
    let x = 12345;
    for (let i = 0; i < n; i++) {
      x = (Math.imul(x, 1103515245) + 12345) >>> 0;
      keys[i] = x % max;
    }
    return keys;
  }

  it('orders the indices by key, equal keys by index, for 1 to 4 passes', () => {
    for (const max of [2 ** 8, 2 ** 16, 2 ** 24, 2 ** 31]) {
      const keys = keysOf(500, max);
      const reference = Array.from(keys.keys()).sort((a, b) => keys[a] - keys[b] || a - b);
      expect(Array.from(radixSortIndices(keys))).toEqual(reference);
    }
  });

  it('orders many keys by bucket as by index, skewed ones too', () => {
    const n = 70_000;
    const check = (keys: Uint32Array, max?: number) => {
      const reference = Array.from(keys.keys()).sort((a, b) => keys[a] - keys[b] || a - b);
      expect(Array.from(radixSortIndices(keys, max))).toEqual(reference);
    };
    for (const max of [2 ** 17, 2 ** 24, 2 ** 31, 2 ** 32 - 1]) check(keysOf(n, max));
    // all in one top bucket but one, and few distinct keys (many equal)
    const skewed = keysOf(n, 2 ** 10);
    skewed[17] = 2 ** 27 - 1;
    check(skewed, 2 ** 27 - 1);
  });

  it('keeps equal keys in index order, the keys unchanged', () => {
    const keys = Int32Array.of(3, 1, 3, 0, 1, 3);
    expect(Array.from(radixSortIndices(keys, 3))).toEqual([3, 1, 4, 0, 2, 5]);
    expect(Array.from(keys)).toEqual([3, 1, 3, 0, 1, 3]);
  });

  it('handles empty, single and all-zero keys', () => {
    expect(Array.from(radixSortIndices(new Uint32Array(0)))).toEqual([]);
    expect(Array.from(radixSortIndices(Uint32Array.of(7)))).toEqual([0]);
    expect(Array.from(radixSortIndices(new Uint32Array(4)))).toEqual([0, 1, 2, 3]);
  });
});

describe('countingSortIndices', () => {
  it('orders the indices as radixSortIndices, equal keys by index', () => {
    const keys = new Uint32Array(5000);
    let x = 4321;
    for (let i = 0; i < keys.length; i++) {
      x = (Math.imul(x, 1103515245) + 12345) >>> 0;
      keys[i] = x % 65536;
    }
    const out = new Uint32Array(keys.length);
    countingSortIndices(keys, keys.length, 65535, out);
    expect(Array.from(out)).toEqual(Array.from(radixSortIndices(keys, 65535)));
  });

  it('sorts the first n keys into out from an offset, keys of any array', () => {
    const out = new Int32Array(7).fill(-1);
    countingSortIndices([3, 1, 3, 0, 1, 3, 9], 6, 3, out, 1);
    expect(Array.from(out)).toEqual([-1, 3, 1, 4, 0, 2, 5]);
    // a smaller range after a larger one: its counts reset
    countingSortIndices([1, 0], 2, 1, out);
    expect(Array.from(out.subarray(0, 2))).toEqual([1, 0]);
  });

  it('handles no keys and one', () => {
    const out = new Uint32Array(1);
    countingSortIndices([], 0, 10, out);
    countingSortIndices([5], 1, 10, out);
    expect(Array.from(out)).toEqual([0]);
  });
});
