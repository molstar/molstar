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

  it('reads any range of many parts, empty parts and gaps among them', () => {
    // part i holds i at offset 3i, one element (none where i % 5 === 0), gaps of two
    const parts: SegmentedArray.Part<Uint32Array>[] = [];
    for (let i = 0; i < 1000; ++i) {
      parts.push({ array: new Uint32Array([i]), offset: 3 * i, count: i % 5 === 0 ? 0 : 1, bias: 0, version: 0 });
    }
    const s = SegmentedArray.create(Uint32Array, parts);
    const all = new Uint32Array(s.length);
    for (let i = 0; i < 1000; ++i) if (i % 5 !== 0) all[3 * i] = i;
    for (const [start, end] of [
      [0, 1],
      [2, 4],
      [3, 4],
      [1500, 1503],
      [2995, s.length],
      [0, s.length],
    ]) {
      expect(Array.from(SegmentedArray.read(s, start, end))).toEqual(Array.from(all.subarray(start, end)));
    }
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

  it('writes the parts an extension added alone, also of a line of them; another array all', () => {
    const s = create();
    const writes: [number[], number][] = [];
    const write = (data: Uint32Array, offset: number) => writes.push([Array.from(data), offset]);
    let uploaded = SegmentedArray.writeChanged(s, undefined, write);
    expect(writes).toHaveLength(2);
    writes.length = 0;
    // a part after the others, within the length: the same parts first, the same length
    const c = { array: new Uint32Array([7]), offset: 6, count: 1, bias: 0, version: 0 };
    const e1 = SegmentedArray.extend(s, [c]);
    expect(e1.length).toBe(s.length);
    expect(e1.parts.slice(0, 2)).toEqual(s.parts);
    expect(e1.from).toBe(s);
    expect(Array.from(SegmentedArray.flatten(e1))).toEqual([1, 2, 3, 0, 14, 15, 7, 0]);
    uploaded = SegmentedArray.writeChanged(e1, uploaded, write);
    expect(writes).toEqual([[[7], 6]]);
    writes.length = 0;
    // two extensions since the last upload: both parts
    const e3 = SegmentedArray.extend(SegmentedArray.extend(e1, [{ ...c, offset: 7 }]), []);
    expect(SegmentedArray.changedRanges(e3, { array: e1, versions: [0, 0, 0] })).toEqual([7, 8]);
    // not after the others, or past the length: refused
    expect(() => SegmentedArray.extend(s, [{ ...c, offset: 5 }])).toThrow();
    expect(() => SegmentedArray.extend(s, [{ ...c, offset: 8 }])).toThrow();
    // longer: room after the parts, which stay where they were
    const longer = SegmentedArray.extend(e1, [{ ...c, offset: 9 }], 12);
    expect(longer.length).toBe(12);
    expect(longer.byteLength).toBe(48);
    expect(longer.from).toBe(e1);
    expect(Array.from(SegmentedArray.flatten(longer))).toEqual([1, 2, 3, 0, 14, 15, 7, 0, 0, 7, 0, 0]);
    expect(SegmentedArray.changedRanges(longer, { array: e1, versions: [0, 0, 0] })).toEqual([9, 10]);
    expect(() => SegmentedArray.extend(e1, [], 4)).toThrow();
    // a long line: cut behind its last extensions, an uploader one behind still writes the added part alone
    let line: SegmentedArray = s;
    for (let k = 0; k < 40; ++k) line = SegmentedArray.extend(line, []);
    let depth = 0;
    for (let x = line.from; x; x = x.from) depth += 1;
    expect(depth).toBeLessThanOrEqual(16);
    const next = SegmentedArray.extend(line, [c]);
    expect(SegmentedArray.changedRanges(next, { array: line, versions: [0, 0] })).toEqual([6, 7]);
    // another array: all
    uploaded = SegmentedArray.writeChanged(create(), uploaded, write);
    expect(writes).toHaveLength(2);
  });
});
