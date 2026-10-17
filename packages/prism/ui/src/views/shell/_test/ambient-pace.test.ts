/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { ActiveMs, AmbientPace, IdleStepMs } from '../ambient-pace.js';

// The splash's pace: full rate for a while after input, the idle pace otherwise, and full rate while the glass animates
// (a shatter, a reassembly that outlasts the window).
describe('AmbientPace', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => vi.useRealTimers());

  it('steps at the idle pace after the active window, unless something animates', () => {
    vi.setSystemTime(10_000);
    const pace = new AmbientPace();
    pace.input();
    pace.stepped();
    vi.advanceTimersByTime(ActiveMs + 100);
    pace.stepped();
    vi.advanceTimersByTime(IdleStepMs / 2);
    expect(pace.due()).toBe(false);
    expect(pace.due(true)).toBe(true);
    vi.advanceTimersByTime(IdleStepMs);
    expect(pace.due()).toBe(true);
  });
});
