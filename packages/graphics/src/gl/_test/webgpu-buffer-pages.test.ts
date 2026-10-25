/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createRangeAllocator } from '../range-allocator.js';

describe('range allocator', () => {
  it('allocates aligned ranges first fit, until full', () => {
    const a = createRangeAllocator(64, 16);
    expect(a.allocate(1)).toBe(0);
    expect(a.allocate(17)).toBe(16);
    expect(a.allocate(16)).toBe(48);
    expect(a.used).toBe(64);
    expect(a.allocate(1)).toBe(-1);
  });

  it('merges freed ranges with their free neighbors', () => {
    const a = createRangeAllocator(64, 16);
    const offsets = [a.allocate(16), a.allocate(16), a.allocate(16), a.allocate(16)];
    expect(offsets).toEqual([0, 16, 32, 48]);
    a.free(16, 16);
    a.free(48, 16);
    // not adjacent: no 32 bytes in one piece
    expect(a.allocate(32)).toBe(-1);
    a.free(32, 16);
    // 16..64 merged
    expect(a.allocate(48)).toBe(16);
    a.free(16, 48);
    a.free(0, 16);
    expect(a.used).toBe(0);
    expect(a.allocate(64)).toBe(0);
  });

  it('reuses the first fitting gap', () => {
    const a = createRangeAllocator(128, 16);
    for (let i = 0; i < 8; ++i) a.allocate(16);
    a.free(32, 32);
    a.free(96, 16);
    expect(a.allocate(16)).toBe(32);
    expect(a.allocate(16)).toBe(48);
    expect(a.allocate(16)).toBe(96);
    expect(a.allocate(16)).toBe(-1);
  });
});
