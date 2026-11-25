/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { grown } from '../draw-lists.js';

describe('grown', () => {
  it('keeps what was written while a buffer grows as it fills', () => {
    let a: Float32Array = new Float32Array(4);
    for (let i = 0; i < 100; i++) {
      a = grown(a, i + 1);
      a[i] = i;
    }
    expect(Array.from(a.subarray(0, 100))).toEqual(Array.from({ length: 100 }, (_, i) => i));
  });

  it('returns the same array when it is large enough, and keeps its type', () => {
    const a = new Int32Array(8);
    expect(grown(a, 8)).toBe(a);
    expect(grown(a, 9)).toBeInstanceOf(Int32Array);
  });
});
