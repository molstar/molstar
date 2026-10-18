/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterEach, beforeEach, describe, expect, it, vi, type MockInstance } from 'vitest';
import { isDebugMode, setDebugMode } from '@molstar/core/util/debug';
import { logGpuFallback } from '../util.js';

describe('gpu fallback log', () => {
  const debug = isDebugMode;
  let log: MockInstance;
  beforeEach(() => {
    log = vi.spyOn(console, 'log').mockImplementation(() => {});
  });
  afterEach(() => {
    log.mockRestore();
    setDebugMode(debug);
  });

  it('logs nothing outside of debug mode', () => {
    setDebugMode(false);
    logGpuFallback('Test surface', 'quiet reason');
    expect(log).not.toHaveBeenCalled();
  });

  it('logs each reason once in debug mode', () => {
    setDebugMode(true);
    logGpuFallback('Test surface', 'first reason');
    logGpuFallback('Test surface', 'first reason');
    logGpuFallback('Test surface', 'second reason');
    expect(log).toHaveBeenCalledTimes(2);
    expect(log).toHaveBeenNthCalledWith(1, 'Test surface: computed on the CPU, not the GPU: first reason');
  });
});
