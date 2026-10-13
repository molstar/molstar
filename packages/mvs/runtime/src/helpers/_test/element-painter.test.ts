/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ElementPainter } from '../element-painter.js';
import { GlobalResidueIndices } from '../indexing.js';
import { Column } from '@molstar/core/data/db';

function paintAll(n: number, ranges: [number, number, number][]) {
  const painter = new ElementPainter();
  painter.begin(n);
  for (const [from, to, value] of ranges) painter.paint([from], [to], 1, value);
  return painter;
}

describe('ElementPainter', () => {
  it('returns null when nothing was painted and when there are no elements', () => {
    expect(paintAll(10, []).result()).toBeNull();
    const p = new ElementPainter();
    p.begin(0);
    expect(p.isComplete).toBe(true);
    p.paint([0], [0], 1, 3);
    expect(p.result()).toBeNull();
  });
  it('the first painted range wins on overlaps, unpainted elements are -1', () => {
    const p = paintAll(10, [
      [2, 6, 7],
      [4, 9, 8],
      [0, 3, 9],
    ]);
    expect(Array.from(p.result()!)).toEqual([9, 9, 7, 7, 7, 7, 8, 8, 8, -1]);
    expect(p.isComplete).toBe(false);
  });
  it('knows when every element is painted, and then ignores further ranges', () => {
    const p = paintAll(6, [
      [0, 4, 1],
      [4, 6, 2],
    ]);
    expect(p.isComplete).toBe(true);
    p.paint([0], [6], 1, 99);
    expect(Array.from(p.result()!)).toEqual([1, 1, 1, 1, 2, 2]);
  });
  it('begin() forgets the previous result and the shared buffer is reset', () => {
    const p = new ElementPainter();
    p.begin(5);
    p.paint([0], [5], 1, 1);
    expect(p.isComplete).toBe(true);
    p.begin(8);
    expect(p.isComplete).toBe(false);
    expect(p.result()).toBeNull();
    p.paint([2], [4], 1, 5);
    expect(Array.from(p.result()!)).toEqual([-1, -1, 5, 5, -1, -1, -1, -1]);
  });
  it('matches brute force (first range wins) on random ranges', () => {
    let seed = 7;
    const rnd = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    for (let iter = 0; iter < 200; iter++) {
      const n = 1 + Math.floor(rnd() * 60);
      const ranges: [number, number, number][] = Array.from({ length: Math.floor(rnd() * 12) }, (_, i) => {
        const from = Math.floor(rnd() * n);
        return [from, Math.min(n, from + Math.floor(rnd() * 20)), i];
      });
      const expected = new Array(n).fill(-1);
      for (const [from, to, v] of ranges) for (let i = from; i < to; i++) if (expected[i] < 0) expected[i] = v;
      const p = paintAll(n, ranges);
      expect(Array.from(p.result() ?? new Array(n).fill(-1))).toEqual(expected);
      expect(p.isComplete).toBe(expected.every((v) => v >= 0));
    }
  });
});

describe('GlobalResidueIndices', () => {
  const hierarchy = (sourceIndex: number[]) =>
    ({
      residues: { _rowCount: sourceIndex.length },
      residueSourceIndex: Column.ofIntArray(sourceIndex),
    }) as any;

  it('maps residue_index to the residue when it is unique', () => {
    const map = new GlobalResidueIndices(hierarchy([10, 11, 12, 40])).residueBySourceIndex!;
    expect(map.get(12)).toBe(2);
    expect(map.get(40)).toBe(3);
    expect(map.get(13)).toBeUndefined();
  });
  it('is null when residue_index repeats (the per-chain indices are used then)', () => {
    expect(new GlobalResidueIndices(hierarchy([0, 1, 1, 2])).residueBySourceIndex).toBeNull();
  });
});
