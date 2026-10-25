/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { SegmentedArray } from '../segmented-array.js';

function create() {
  // elements 0-2 of a, a gap of one, 4-5 of b biased by 10, zeros up to 8
  const a = new Uint32Array([1, 2, 3, 99]);
  const b = new Uint32Array([4, 5]);
  return SegmentedArray.create(
    Uint32Array,
    [
      { array: a, offset: 0, count: 3, bias: 0, version: 0 },
      { array: b, offset: 4, count: 2, bias: 10, version: 0 },
    ],
    8,
  );
}

describe('segmented array', () => {
  it('composes ranges from its parts, gaps zero, biased where a part has a bias', () => {
    const s = create();
    expect(s.length).toBe(8);
    expect(s.byteLength).toBe(32);
    expect(Array.from(SegmentedArray.flatten(s))).toEqual([1, 2, 3, 0, 14, 15, 0, 0]);
    expect(Array.from(SegmentedArray.read(s, 2, 5))).toEqual([3, 0, 14]);
  });

  it('writes all parts first, then the parts whose version changed', () => {
    const s = create();
    const writes: [number[], number][] = [];
    const write = (data: Uint32Array, offset: number) => writes.push([Array.from(data), offset]);
    const uploaded = SegmentedArray.writeChanged(s, undefined, write);
    expect(writes).toEqual([
      [[1, 2, 3], 0],
      [[14, 15], 4],
    ]);

    writes.length = 0;
    expect(SegmentedArray.writeChanged(s, uploaded, write)).toBe(uploaded);
    expect(writes).toEqual([]);

    (s.parts[1] as SegmentedArray.Part).array = new Uint32Array([6, 7]);
    (s.parts[1] as SegmentedArray.Part).version += 1;
    SegmentedArray.writeChanged(s, uploaded, write);
    expect(writes).toEqual([[[16, 17], 4]]);
  });

  it('hands out the ranges of changed parts, all for another array', () => {
    const s = create();
    const uploaded = SegmentedArray.uploaded(s);
    expect(SegmentedArray.changedRanges(s, uploaded)).toEqual([]);
    (s.parts[0] as SegmentedArray.Part).version += 1;
    expect(SegmentedArray.changedRanges(s, uploaded)).toEqual([0, 3]);
    expect(SegmentedArray.changedRanges(s, uploaded)).toEqual([]);
    expect(SegmentedArray.changedRanges(create(), uploaded)).toBeUndefined();
  });
});
