/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ChunkedArray } from '../chunked-array.js';

describe('Chunked Array', () => {
  it('creation', () => {
    const arr = ChunkedArray.create<number, 2>(Array, 2, 2);
    ChunkedArray.add2(arr, 1, 2);
    ChunkedArray.add2(arr, 3, 4);
    expect(ChunkedArray.compact(arr)).toEqual([1, 2, 3, 4]);
  });

  it('initial', () => {
    const arr = ChunkedArray.create(Int32Array, 2, 6, new Int32Array([1, 2, 3, 4]));
    ChunkedArray.add2(arr, 4, 3);
    ChunkedArray.add2(arr, 2, 1);
    ChunkedArray.add2(arr, 5, 6);
    expect(ChunkedArray.compact(arr)).toEqual(new Int32Array([4, 3, 2, 1, 5, 6]));
  });

  it('add many', () => {
    const arr = ChunkedArray.create<number, 2>(Array, 2, 2);
    ChunkedArray.addMany(arr, [1, 2, 3, 4]);
    expect(ChunkedArray.compact(arr)).toEqual([1, 2, 3, 4]);
  });

  it('resize', () => {
    const arr = ChunkedArray.create<number, 2>(Int32Array, 2, 2);
    ChunkedArray.add2(arr, 1, 2);
    ChunkedArray.add2(arr, 3, 4);
    ChunkedArray.add2(arr, 5, 6);
    ChunkedArray.add2(arr, 7, 8);
    ChunkedArray.add2(arr, 9, 10);
    expect(ChunkedArray.compact(arr)).toEqual(new Int32Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]));
  });

  it('resize-fraction', () => {
    const arr = ChunkedArray.create<number, 2>(Int32Array, 2, 2.5);
    ChunkedArray.add2(arr, 1, 2);
    ChunkedArray.add2(arr, 3, 4);
    ChunkedArray.add2(arr, 5, 6);
    ChunkedArray.add2(arr, 7, 8);
    ChunkedArray.add2(arr, 9, 10);
    expect(arr.elementCount).toBe(5);
    expect(ChunkedArray.compact(arr)).toEqual(new Int32Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]));
  });

  it('copies the used part of an underfilled last chunk', () => {
    const typed = ChunkedArray.create<number, 1>(Int32Array, 1, 4);
    for (let i = 1; i <= 5; i++) ChunkedArray.add(typed, i);
    expect(ChunkedArray.compact(typed)).toEqual(new Int32Array([1, 2, 3, 4, 5]));
    const plain = ChunkedArray.create<number, 1>(Array, 1, 4);
    for (let i = 1; i <= 5; i++) ChunkedArray.add(plain, i);
    expect(ChunkedArray.compact(plain)).toEqual([1, 2, 3, 4, 5]);
  });

  it('trims an underfilled single chunk', () => {
    const arr = ChunkedArray.create<number, 3>(Float32Array, 3, 4, 10);
    for (let i = 0; i < 7; i++) ChunkedArray.add3(arr, i, i, i);
    // 7 of 10 used: below 0.8, a copy of the used part
    const trimmed = ChunkedArray.compactTrimmed(arr);
    expect(trimmed.length).toBe(21);
    expect(trimmed).toEqual(ChunkedArray.compact(arr, true).slice(0, 21));
    // at least `minFill` used: the chunk as it is
    expect(ChunkedArray.compactTrimmed(arr, 0.7)).toBe(arr.chunks[0]);
    ChunkedArray.add3(arr, 7, 7, 7);
    expect(ChunkedArray.compactTrimmed(arr)).toBe(arr.chunks[0]);
    // several chunks: compacted to the used size as before
    for (let i = 8; i < 12; i++) ChunkedArray.add3(arr, i, i, i);
    expect(ChunkedArray.compactTrimmed(arr).length).toBe(36);
  });

  it('is (perf-audit-2 01/B5)', () => {
    expect(ChunkedArray.is(ChunkedArray.create(Float32Array, 3, 4))).toBe(true);
    expect(ChunkedArray.is(ChunkedArray.create(Array, 1, 4))).toBe(true);
    expect(ChunkedArray.is([1, 2, 3])).toBe(false);
    expect(ChunkedArray.is({})).toBe(false);
    expect(ChunkedArray.is(undefined)).toBe(false);
  });
});
