/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subject } from 'rxjs';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { watchedState } from '../watched.js';

// A source that follows its events only while it is watched, and reads a burst of them once per frame.

describe('watchedState', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => vi.useRealTimers());

  const make = (coalesce: 'frame' | 'none' = 'frame') => {
    const events = new Subject<void>();
    let n = 0;
    let reads = 0;
    const state = watchedState(
      () => {
        reads++;
        return n;
      },
      () => events,
      { coalesce },
    );
    return {
      events,
      state,
      reads: () => reads,
      bump: () => {
        n++;
        events.next();
      },
    };
  };

  it('does not listen until it is subscribed, and reads `value` fresh', () => {
    const s = make();
    expect(s.events.observed).toBe(false);
    s.bump();
    expect(s.events.observed).toBe(false);
    expect(s.state.value).toBe(1);
    const seen: number[] = [];
    const sub = s.state.subscribe((x) => seen.push(x));
    expect(s.events.observed).toBe(true);
    expect(seen).toEqual([1]);
    sub.unsubscribe();
    expect(s.events.observed).toBe(false);
  });

  it('reads a burst of events once per frame, and `value` never waits for it', () => {
    const s = make();
    const seen: number[] = [];
    const sub = s.state.subscribe((x) => seen.push(x));
    const before = s.reads();
    s.bump();
    s.bump();
    s.bump();
    expect(s.reads()).toBe(before);
    // Asked for in the meantime, it is current (one read), and the frame has nothing left to do.
    expect(s.state.value).toBe(3);
    expect(s.reads()).toBe(before + 1);
    vi.advanceTimersByTime(40);
    expect(s.reads()).toBe(before + 1);
    expect(seen).toEqual([0, 3]);

    s.bump();
    s.bump();
    vi.advanceTimersByTime(40);
    expect(s.reads()).toBe(before + 2);
    expect(seen).toEqual([0, 3, 5]);
    sub.unsubscribe();
  });

  it('reads at once where it is told not to wait, and stops with the last subscriber', () => {
    const s = make('none');
    const a = s.state.subscribe();
    const b = s.state.subscribe();
    s.bump();
    expect(s.state.value).toBe(1);
    a.unsubscribe();
    expect(s.events.observed).toBe(true);
    b.unsubscribe();
    expect(s.events.observed).toBe(false);
    // Whatever happened while nobody looked is read when someone does.
    s.bump();
    const seen: number[] = [];
    s.state.subscribe((x) => seen.push(x)).unsubscribe();
    expect(seen).toEqual([2]);
  });
});
