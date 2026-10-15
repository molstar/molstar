/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { IntBuffer, MaxRetainedBufferCapacity, ScratchMemory } from '../int-buffer.js';

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

describe('ScratchMemory', () => {
  it('releases the holders of a context, or all of them, until they are unregistered', () => {
    const a = vi.fn();
    const b = vi.fn();
    const c = vi.fn();
    const unregisterA = ScratchMemory.register('test-a', a);
    ScratchMemory.register('test-a', b);
    const unregisterC = ScratchMemory.register('test-c', c);
    ScratchMemory.release('test-a');
    expect([a.mock.calls.length, b.mock.calls.length, c.mock.calls.length]).toEqual([1, 1, 0]);
    ScratchMemory.release('test-c');
    expect([a.mock.calls.length, b.mock.calls.length, c.mock.calls.length]).toEqual([1, 1, 1]);
    ScratchMemory.release('nothing'); // (no such context: nothing to do)
    ScratchMemory.release('test-a', 'test-c', 'nothing'); // (several contexts)
    expect([a.mock.calls.length, b.mock.calls.length, c.mock.calls.length]).toEqual([2, 2, 2]);
    ScratchMemory.release();
    expect([a.mock.calls.length, b.mock.calls.length, c.mock.calls.length]).toEqual([3, 3, 3]);
    unregisterA();
    ScratchMemory.release('test-a');
    expect([a.mock.calls.length, b.mock.calls.length]).toEqual([3, 4]);
    unregisterC();
    ScratchMemory.release();
    expect(c.mock.calls.length).toBe(3);
  });
  it('keeps what is small and frees what a large call grew', () => {
    const b = new IntBuffer(4);
    for (let i = 0; i < MaxRetainedBufferCapacity + 1; i++) b.push(i);
    expect(b.data.length).toBeGreaterThan(MaxRetainedBufferCapacity);
    b.trim(MaxRetainedBufferCapacity);
    expect(b.data.length).toBeLessThanOrEqual(MaxRetainedBufferCapacity);
    for (let i = 0; i < 1000; i++) b.push(i);
    const kept = b.data;
    b.trim(MaxRetainedBufferCapacity);
    expect(b.data).toBe(kept);
  });
});
