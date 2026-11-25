/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { requestFrameCallback } from '../frame.js';

describe('dom: the frame lane', () => {
  afterEach(() => {
    vi.useRealTimers();
  });

  it('asks again for a tick asked on a clock that is gone', () => {
    // a tick asked on fake timers that are uninstalled before it runs never runs; the lane must not wait for it
    vi.useFakeTimers();
    requestFrameCallback(() => {});
    vi.useRealTimers();

    vi.useFakeTimers();
    const ran = vi.fn();
    requestFrameCallback(ran);
    vi.advanceTimersByTime(50);
    expect(ran).toHaveBeenCalledTimes(1);
  });
});
