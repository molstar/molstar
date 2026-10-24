/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { ChunkedArray } from '@molstar/core/data/util';
import { MarchingCubesLinesBuilder } from '../builder.js';

describe('MarchingCubesLinesBuilder', () => {
  // a Float32Array rounds ids above 2^24: a wireframe of more than 16.7M vertices would join the wrong ones (perf-audit-2
  // 08/B4); exact ids need a vertex array that large, so the test checks the kind the ids are stored in
  it('stores the line vertex ids as integers', () => {
    const create = vi.spyOn(ChunkedArray, 'create');
    MarchingCubesLinesBuilder(1024);
    const pairs = create.mock.calls.filter((c) => c[1] === 2);
    expect(pairs.length).toBe(1);
    expect(pairs[0][0]).toBe(Uint32Array);
    expect(Math.fround(2 ** 24 + 1)).not.toBe(2 ** 24 + 1);
    expect(new Uint32Array([2 ** 24 + 1])[0]).toBe(2 ** 24 + 1);
  });
});
