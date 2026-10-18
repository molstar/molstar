/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { ActiveMs, AmbientPace, FadeMs, RestAfterMs } from '../ambient-pace.js';

// The splash's rest and wake (prism AGENTS.md, idle is idle): awake while input comes, the glass animates or the intro
// plays; then the particles and the wave fade out and the pace sleeps (the field asks for no more frames); input wakes
// it and they fade back in.
describe('AmbientPace', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => vi.useRealTimers());

  /** Frames of 16 ms for `ms`, as the field's tick advances the pace. */
  const run = (pace: AmbientPace, ms: number, animating = false): number => {
    let level = pace.level;
    for (let t = 0; t < ms; t += 16) {
      vi.advanceTimersByTime(16);
      level = pace.advance(animating);
    }
    return level;
  };

  /** A pace as the field starts it: input at the start. */
  const started = (options?: ConstructorParameters<typeof AmbientPace>[0]): AmbientPace => {
    const pace = new AmbientPace(options);
    pace.input();
    pace.advance(false);
    return pace;
  };

  it('lets the intro play: at least ten seconds before it may rest', () => {
    expect(RestAfterMs).toBeGreaterThanOrEqual(10_000);
    vi.setSystemTime(10_000);
    const pace = started();
    expect(run(pace, RestAfterMs - 500)).toBe(1);
    expect(pace.asleep).toBe(false);
    run(pace, 500 + FadeMs / 2);
    expect(pace.level).toBeLessThan(1);
    expect(pace.level).toBeGreaterThan(0);
  });

  it('rests after the intro and the active window: fades out, then sleeps', () => {
    vi.setSystemTime(10_000);
    const pace = started();
    run(pace, RestAfterMs + FadeMs / 2);
    expect(pace.level).toBeLessThan(1);
    expect(run(pace, FadeMs)).toBe(0);
    expect(pace.asleep).toBe(true);
  });

  it('stays awake under input, and while the glass animates', () => {
    vi.setSystemTime(10_000);
    const pace = started();
    // Input every second, past the intro: never rests.
    for (let i = 0; i < RestAfterMs / 1000 + 3; i++) {
      pace.input();
      expect(run(pace, 1000)).toBe(1);
    }
    // The glass animating (a shatter, a gathering): never rests, however long.
    expect(run(pace, ActiveMs + 5000, true)).toBe(1);
    expect(pace.asleep).toBe(false);
    // Done animating: rests (the active window long passed).
    run(pace, FadeMs + 100);
    expect(pace.asleep).toBe(true);
  });

  it('wakes on input: fades back in from where it was', () => {
    vi.setSystemTime(10_000);
    const pace = started();
    run(pace, RestAfterMs + FadeMs + 100);
    expect(pace.asleep).toBe(true);
    pace.input();
    expect(pace.asleep).toBe(false);
    expect(pace.level).toBe(0);
    pace.reset();
    const level = run(pace, FadeMs / 2);
    expect(level).toBeGreaterThan(0.2);
    expect(level).toBeLessThan(0.9);
    expect(run(pace, FadeMs)).toBe(1);
    // And rests again after the active window.
    run(pace, ActiveMs + FadeMs + 100);
    expect(pace.asleep).toBe(true);
  });

  it('wakes in the middle of a fade out: back up from the level it had', () => {
    vi.setSystemTime(10_000);
    const pace = started();
    run(pace, RestAfterMs);
    const mid = run(pace, FadeMs / 2);
    expect(mid).toBeGreaterThan(0);
    expect(mid).toBeLessThan(1);
    pace.input();
    expect(run(pace, 32)).toBeGreaterThan(mid);
  });

  it('a long stall does not skip the intro or the fade (a frame counts at most 100 ms)', () => {
    vi.setSystemTime(10_000);
    const pace = started();
    vi.advanceTimersByTime(RestAfterMs * 2);
    pace.advance(false);
    // Only one clamped frame of the intro has passed.
    vi.advanceTimersByTime(1000);
    expect(pace.advance(false)).toBe(1);
  });

  it('comes back from a pause (hidden, shown) awake with its window again, asleep stays asleep', () => {
    vi.setSystemTime(10_000);
    const pace = started();
    run(pace, RestAfterMs + 5000);
    expect(pace.asleep).toBe(true);
    pace.reset();
    expect(pace.asleep).toBe(true);

    const other = started({ restAfterMs: 1000 });
    run(other, 1000);
    expect(other.asleep).toBe(false);
    vi.advanceTimersByTime(60_000);
    other.reset();
    expect(run(other, 100)).toBe(1);
  });

  it('counts its steps', () => {
    const pace = new AmbientPace();
    pace.stepped();
    pace.stepped();
    expect(pace.steps).toBe(2);
  });
});
