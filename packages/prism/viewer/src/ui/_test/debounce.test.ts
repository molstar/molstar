/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { debounce, Debouncer } from '../debounce.js';

// The form-edit debounce: the last call wins, flush runs it now, cancel and dispose drop it.

describe('debounce', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => vi.useRealTimers());

  it('runs once, with the last arguments, after the calls stop', () => {
    const fn = vi.fn();
    const d = debounce(fn, 100);
    d(1);
    vi.advanceTimersByTime(60);
    d(2);
    vi.advanceTimersByTime(60);
    expect(fn).not.toHaveBeenCalled();
    expect(d.pending).toBe(true);
    vi.advanceTimersByTime(60);
    expect(fn).toHaveBeenCalledTimes(1);
    expect(fn).toHaveBeenCalledWith(2);
    expect(d.pending).toBe(false);
  });

  it('flushes now, cancels, and is cleared on dispose', () => {
    const fn = vi.fn();
    const d = debounce(fn, 100);
    d('a');
    d.flush();
    expect(fn).toHaveBeenCalledWith('a');
    vi.advanceTimersByTime(500);
    expect(fn).toHaveBeenCalledTimes(1);

    d('b');
    d.cancel();
    d('c');
    d.dispose();
    vi.advanceTimersByTime(500);
    expect(fn).toHaveBeenCalledTimes(1);
    expect(d.pending).toBe(false);
  });
});

describe('Debouncer', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => vi.useRealTimers());

  it('waits per key', () => {
    const run = vi.fn();
    const d = new Debouncer<string>(run, 100);
    d.schedule('a');
    vi.advanceTimersByTime(50);
    d.schedule('b');
    expect(d.has('a') && d.has('b')).toBe(true);
    vi.advanceTimersByTime(60);
    expect(run.mock.calls).toEqual([['a']]);
    d.cancel('b');
    vi.advanceTimersByTime(200);
    expect(run).toHaveBeenCalledTimes(1);
    d.schedule('c');
    d.flushAll();
    expect(run.mock.calls).toEqual([['a'], ['c']]);
    d.schedule('d');
    d.dispose();
    vi.advanceTimersByTime(200);
    expect(run).toHaveBeenCalledTimes(2);
    expect(d.waiting).toBe(false);
  });
});
