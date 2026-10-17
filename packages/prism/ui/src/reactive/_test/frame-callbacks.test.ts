/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterAll, beforeAll, describe, expect, it, vi } from 'vitest';
import { createRoot, requireScope, type Scope } from '../../dom/scope.js';
import { frame, frames } from '../../dom/scoped.js';
import { cancelFrameCallback, pendingFrameCallbacks, requestFrameCallback } from '../frame-state.js';

// Frame callbacks on the page's one tick (rendering.md §2.3). Outside a browser the tick is a 16 ms timer: faked here.

beforeAll(() => {
  vi.useFakeTimers();
});
afterAll(() => {
  vi.useRealTimers();
});
const tick = () => vi.advanceTimersByTime(16);

describe('frame callbacks', () => {
  it('one cancelled by another in the same frame does not run', () => {
    const ran: string[] = [];
    let b = 0;
    requestFrameCallback(() => {
      ran.push('a');
      cancelFrameCallback(b);
    });
    b = requestFrameCallback(() => ran.push('b'));
    tick();
    expect(ran).toEqual(['a']);
    tick();
    expect(ran).toEqual(['a']);
  });

  it('frames().cancel() from another callback in the same frame stops it', () => {
    const ran: string[] = [];
    const dispose = createRoot((d) => {
      const later = frames(() => ran.push('later'));
      const first = frames(() => {
        ran.push('first');
        later.cancel();
      });
      first.request();
      later.request();
      return d;
    });
    tick();
    expect(ran).toEqual(['first']);
    dispose();
  });

  it('frame() paused again and again waits with one request, not one a pause', () => {
    let scope!: Scope;
    let runs = 0;
    const dispose = createRoot((d) => {
      scope = requireScope('test');
      frame(() => runs++);
      return d;
    });
    for (let i = 0; i < 3; i++) {
      scope.setActive(false);
      tick();
      scope.setActive(true);
    }
    expect(pendingFrameCallbacks()).toBe(1);
    tick();
    expect(runs).toBe(1);
    expect(pendingFrameCallbacks()).toBe(0);
    dispose();
  });
});
