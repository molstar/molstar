/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createQuadIndexArray, createQuadIndices, QuadIndexPattern } from '../quad-indices.js';

describe('quad indices', () => {
  it('index the six-vertex quads by four of their ids, as their two triangles', () => {
    const indices = createQuadIndexArray(3);
    expect(Array.from(indices)).toEqual([0, 1, 2, 1, 4, 2, 6, 7, 8, 7, 10, 8, 12, 13, 14, 13, 16, 14]);
    // the corners the shaders map per id (`m = id % 6`), the triangles of the unindexed quad
    const corner = (id: number) =>
      [
        [-1, 1],
        [-1, -1],
        [1, 1],
        [-1, -1],
        [1, -1],
        [1, 1],
      ][id % 6];
    const triangles = (ids: readonly number[]) => [ids.slice(0, 3).map(corner), ids.slice(3, 6).map(corner)];
    expect(triangles(QuadIndexPattern)).toEqual(triangles([0, 1, 2, 3, 4, 5]));
  });

  it('are the same for any item, in a new array per call (none kept between uploads)', () => {
    const small = createQuadIndexArray(2);
    const large = createQuadIndexArray(100_000);
    expect(large.buffer).not.toBe(small.buffer);
    expect(large.length).toBe(600_000);
    expect(Array.from(large.subarray(0, small.length))).toEqual(Array.from(small));
    const i = 99_999;
    expect(Array.from(large.subarray(i * 6, i * 6 + 6))).toEqual(QuadIndexPattern.map((k) => i * 6 + k));
  });

  it('share one buffer, grown by doubling, released when replaced, none beyond the limit', () => {
    const created: number[] = [];
    const released: number[] = [];
    const quads = createQuadIndices(
      (indices) => {
        created.push(indices.length / 6);
        return created.length;
      },
      (buffer) => released.push(buffer),
      24 * 50_000,
    );
    const first = quads.get(10)!;
    expect(first.itemCount).toBe(4096);
    expect(quads.get(4096)).toBe(first);
    const second = quads.get(5000)!;
    expect(second.itemCount).toBe(8192);
    expect(second.id).not.toBe(first.id);
    expect(released).toEqual([first.buffer]);
    expect(quads.byteCount).toBe(8192 * 24);
    // a request past twice the buffer as asked, doubling capped by the limit, past it none
    expect(quads.get(30_000)!.itemCount).toBe(30_000);
    expect(quads.get(40_000)!.itemCount).toBe(50_000);
    expect(quads.get(50_001)).toBeUndefined();
    quads.dispose();
    expect(released.length).toBe(created.length);
    expect(quads.byteCount).toBe(0);
  });
});
