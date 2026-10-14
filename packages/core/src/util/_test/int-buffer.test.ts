/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { IntBuffer } from '../int-buffer.js';

describe('IntBuffer', () => {
  it('grows, keeps contents, and views only the filled part', () => {
    const b = new IntBuffer(2);
    for (let i = 0; i < 10; i++) b.push(i);
    expect(Array.from(b.view())).toEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    b.clear();
    b.push(42);
    expect(Array.from(b.view())).toEqual([42]);
  });
  it('clear keeps the capacity, trim frees it only above the limit', () => {
    const b = new IntBuffer(4);
    for (let i = 0; i < 1000; i++) b.push(i);
    const capacity = b.data.length;
    expect(capacity).toBeGreaterThanOrEqual(1000);
    b.clear();
    expect(b.n).toBe(0);
    expect(b.data.length).toBe(capacity);
    b.trim(capacity);
    expect(b.data.length).toBe(capacity);
    b.trim(capacity - 1);
    expect(b.data.length).toBeLessThan(capacity);
    expect(b.n).toBe(0);
  });
  it('ensureCapacity keeps the filled values', () => {
    const b = new IntBuffer(2);
    b.push(7);
    b.push(8);
    b.ensureCapacity(1000);
    expect(b.data.length).toBeGreaterThanOrEqual(1000);
    expect(Array.from(b.view())).toEqual([7, 8]);
  });
});
