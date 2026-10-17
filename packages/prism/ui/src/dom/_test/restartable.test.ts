/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { runInScope, Scope } from '../scope.js';
import { restartable } from '../scoped.js';

beforeEach(() => vi.useFakeTimers());
afterEach(() => vi.useRealTimers());

describe('restartable', () => {
  it('runs after the wait, and a start replaces the one in progress', () => {
    const fn = vi.fn();
    const t = runInScope(new Scope(), () => restartable(fn));
    t.start(100);
    vi.advanceTimersByTime(60);
    t.start(100);
    vi.advanceTimersByTime(60);
    expect(fn).not.toHaveBeenCalled();
    vi.advanceTimersByTime(40);
    expect(fn).toHaveBeenCalledTimes(1);
  });

  it('does not run while the scope is paused, runs in full after it resumes, and never after it ends', () => {
    const fn = vi.fn();
    const scope = new Scope();
    const t = runInScope(scope, () => restartable(fn));
    t.start(100);
    scope.setActive(false);
    vi.advanceTimersByTime(500);
    expect(fn).not.toHaveBeenCalled();
    t.start(100);
    vi.advanceTimersByTime(500);
    expect(fn).not.toHaveBeenCalled();
    scope.setActive(true);
    vi.advanceTimersByTime(99);
    expect(fn).not.toHaveBeenCalled();
    vi.advanceTimersByTime(1);
    expect(fn).toHaveBeenCalledTimes(1);
    t.start(100);
    scope.dispose();
    vi.advanceTimersByTime(500);
    expect(fn).toHaveBeenCalledTimes(1);
  });

  it('cancels', () => {
    const fn = vi.fn();
    const scope = new Scope();
    const t = runInScope(scope, () => restartable(fn));
    t.start(10);
    t.cancel();
    scope.setActive(false);
    scope.setActive(true);
    vi.advanceTimersByTime(100);
    expect(fn).not.toHaveBeenCalled();
  });
});
