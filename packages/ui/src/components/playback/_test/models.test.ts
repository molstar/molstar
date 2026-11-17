/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { Stepper } from '../stepper-model.js';
import {
  Durationless,
  formatElapsed,
  Timeline,
  type TimelineModel,
  timelineKind,
  timelinePosition,
  timelineShown,
  timelineValueText,
} from '../timeline-model.js';

// The playback models in plain Node (playback §2, ui §0.1): stepping, the three kinds, what the rows say, ticks only
// while playing (fake timers).

beforeEach(() => {
  vi.useFakeTimers();
});
afterEach(() => {
  vi.useRealTimers();
});

const steps = (n: number) =>
  Array.from({ length: n }, (_, i) => ({ title: `Step ${i + 1}`, description: `# Step ${i + 1}\nText ${i + 1}` }));

describe('the stepper', () => {
  it('steps, wraps and goes to a step', () => {
    const s = new Stepper({ steps: steps(3) });
    expect([s.count.value, s.index.value, s.title.value]).toEqual([3, 0, 'Step 1']);
    s.previous();
    expect(s.index.value).toBe(2);
    s.next();
    expect(s.index.value).toBe(0);
    s.goTo(9);
    expect(s.title.value).toBe('Step 3');
    expect(s.description.value).toBe('# Step 3\nText 3');
    s.openLink('#2');
    expect(s.index.value).toBe(1);
  });

  it('plays every interval from two steps on, and ticks only while playing', () => {
    const one = new Stepper({ steps: steps(1) });
    one.togglePlay();
    expect(one.playing.value).toBe(false);
    expect(vi.getTimerCount()).toBe(0);
    const s = new Stepper({ steps: steps(3), intervalMs: 1000, loop: false });
    s.togglePlay();
    expect(vi.getTimerCount()).toBe(1);
    vi.advanceTimersByTime(1000);
    expect(s.index.value).toBe(1);
    vi.advanceTimersByTime(2000);
    // The last step without looping: it stops there.
    expect([s.index.value, s.playing.value]).toEqual([2, false]);
    expect(vi.getTimerCount()).toBe(0);
  });

  it('keeps its place when the steps change, none with none', () => {
    const s = new Stepper({ steps: steps(4), index: 3 });
    s.setSteps(steps(2));
    expect(s.index.value).toBe(1);
    s.setSteps([]);
    expect([s.index.value, s.title.value]).toEqual([-1, undefined]);
  });
});

describe('a timed timeline', () => {
  it('seeks within its length, plays at fps times its speed, wraps when periodic', () => {
    const seen: number[] = [];
    const t = new Timeline({ duration: 2, fps: 10, speeds: [1, 2], onSeek: (at) => seen.push(at) });
    expect(timelineKind(t)).toBe('timed');
    t.seek(5);
    expect(t.time!.value).toBe(2);
    t.seek(0.4);
    expect(timelinePosition(t)).toBe('0.4 / 2.0 s');
    expect(timelineValueText(t)).toBe('0.4 s of 2.0 s');
    t.setSpeed(2);
    t.togglePlay();
    vi.advanceTimersByTime(500);
    expect(t.time!.value).toBeCloseTo(1.4);
    vi.advanceTimersByTime(400);
    // 0.4 + 0.2 × 9 = 2.2: past the end, round to 0.2.
    expect(t.time!.value).toBeCloseTo(0.2);
    t.togglePlay();
    expect(vi.getTimerCount()).toBe(0);
    expect(seen.length).toBeGreaterThan(5);
  });

  it('stops at the end unless periodic, and plays again from the start', () => {
    const t = new Timeline({ duration: 1, fps: 10, periodic: false });
    t.togglePlay();
    vi.advanceTimersByTime(2000);
    expect([t.time!.value, t.playing.value]).toEqual([1, false]);
    t.togglePlay();
    expect(t.time!.value).toBe(0);
    t.dispose();
  });

  it('is not shown without a length', () => {
    const t = new Timeline({ duration: 0 });
    expect(timelineShown(t)).toBe(false);
    t.setLength(3);
    expect(timelineShown(t)).toBe(true);
  });
});

describe('a counted timeline', () => {
  it('seeks whole frames and says them as people count', () => {
    const t = new Timeline({ count: 10, fps: 5 });
    expect(timelineKind(t)).toBe('counted');
    t.seek(2.4);
    expect(t.index!.value).toBe(2);
    expect(timelinePosition(t)).toBe('3 / 10');
    expect(timelineValueText(t)).toBe('frame 3 of 10');
    expect(timelineValueText(new Timeline({ count: 4, unit: 'model' }))).toBe('model 1 of 4');
    t.togglePlay();
    vi.advanceTimersByTime(1000);
    expect(t.index!.value).toBe(7);
    vi.advanceTimersByTime(600);
    // Past the last frame to the first.
    expect(t.index!.value).toBe(0);
    t.dispose();
    expect(timelineShown(new Timeline({ count: 1 }))).toBe(false);
  });
});

describe('a durationless timeline', () => {
  it('runs until paused: elapsed and steps, said as running or paused', () => {
    const d = new Durationless({ stepMs: 100, speeds: [0.5, 1, 2] });
    expect(timelineKind(d)).toBe('durationless');
    expect(timelineShown(d)).toBe(true);
    expect(timelineValueText(d)).toBe('paused at 0 s, step 0');
    d.togglePlay();
    vi.advanceTimersByTime(12_000);
    expect(d.steps!.value).toBe(120);
    expect(timelineValueText(d)).toBe('running, 12 s, step 120');
    expect(timelinePosition(d)).toBe('12 s · step 120');
    d.setSpeed(2);
    vi.advanceTimersByTime(1000);
    expect(d.steps!.value).toBe(140);
    d.togglePlay();
    expect(vi.getTimerCount()).toBe(0);
  });

  it('steps once only while paused, restarts, and says what it has', () => {
    const d = new Durationless({ elapsed: false });
    d.stepOnce!();
    d.stepOnce!();
    expect(timelineValueText(d)).toBe('paused at step 2');
    d.togglePlay();
    d.stepOnce!();
    expect(d.steps!.value).toBe(2);
    d.restart!();
    expect([d.steps!.value, d.playing.value]).toEqual([0, true]);
    d.dispose();
    const spin = new Durationless({ steps: false, restart: false, stepOnce: false });
    expect([spin.restart, spin.stepOnce, spin.steps]).toEqual([undefined, undefined, undefined]);
    spin.togglePlay();
    vi.advanceTimersByTime(12_500);
    expect(timelineValueText(spin)).toBe('running, 12 s');
    spin.dispose();
  });

  it('a model of its own (no reference): only what it has', () => {
    const playing = { value: false } as TimelineModel['playing'];
    const m: TimelineModel = { playing, togglePlay: () => {} };
    expect([timelineKind(m), timelineValueText(m), timelinePosition(m)]).toEqual(['durationless', 'paused', '']);
  });
});

describe('time run', () => {
  it('is whole seconds, then minutes and hours', () => {
    expect([formatElapsed(0), formatElapsed(12.9), formatElapsed(65), formatElapsed(3725)]).toEqual([
      '0 s',
      '12 s',
      '1:05',
      '1:02:05',
    ]);
  });
});
