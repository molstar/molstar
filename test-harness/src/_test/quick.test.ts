/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick case's control scaled to the case: a step of kilobytes gets a control of 3% of it (4 KB at least)
// and is checked at three quarters of that; a step that allocates nothing keeps `leak`'s 88 bytes (the rule, from given
// readings: the unit project's windows are noisier than the perf tier's). One measured case, on a stand-in step.
// @module-tag medium

import { describe, expect, it } from 'vitest';
import { LEAK_BYTES, leakingBytes, measureBytesPerStep, SCALED_MIN, SCALED_SHARE, scaledBytes } from '../quick.js';

const options = { warmSteps: 1_000, steps: 500, windows: 3, maxWarmSteps: 1_000 };
const sink: unknown[] = [null];

/** A step of about `bytes` a step: an array of that size (garbage at the next step), its slots varying with the index. */
const allocating = (bytes: number) => (i: number) => {
  sink[0] = new Array(Math.round((bytes - 16) / 8)).fill(i);
};

const read = (step: (i: number) => void) => measureBytesPerStep(step, options);

describe('scaledBytes', () => {
  const r = (bytesPerStep: number, spread: number) => ({ bytesPerStep, spread });

  it('keeps the 88 bytes for a step that allocates little and steadily, boundaries included', () => {
    expect(scaledBytes(r(0, 0))).toBe(LEAK_BYTES);
    expect(scaledBytes(r(1_024, 32))).toBe(LEAK_BYTES);
    expect(scaledBytes(r(0, 6))).toBe(LEAK_BYTES);
  });

  it('scales to 3% of a step above a kilobyte or spreading above 32 bytes, 4 KB at least', () => {
    expect(scaledBytes(r(1_025, 0))).toBe(SCALED_MIN);
    expect(scaledBytes(r(10, 33))).toBe(SCALED_MIN);
    expect(scaledBytes(r(1_000_000, 50))).toBe(30_000);
  });

  it('is 3% of the reversed pan, whose windows are 0.4 KB apart and whose control read 2.3 KB of 4 KB (2026-11-24)', () => {
    expect(scaledBytes(r(389_722, 436))).toBe(Math.ceil(SCALED_SHARE * 389_722));
    // the readings given: at 11.7 KB the offset of -2.3 KB leaves 9.4 KB, over three quarters (8.8 KB)
    expect(11_692 - 2_374).toBeGreaterThan(0.75 * scaledBytes(r(389_722, 436)));
  });

  it('is eight spreads where that is more', () => {
    expect(scaledBytes(r(1_000_000, 25_000))).toBe(200_000);
  });
});

describe('the control by the case', () => {
  it('a step of kilobytes: the scaled control is read, and one a hundredth of the step is not lost', () => {
    const step = allocating(20_000);
    const r = read(step);
    const bytes = scaledBytes(r);
    expect(bytes).toBeGreaterThanOrEqual(SCALED_MIN);
    const c = read(leakingBytes(step, bytes));
    expect(c.bytesPerStep - r.bytesPerStep).toBeGreaterThan(0.75 * bytes);
  });

  it('a step of kilobytes: a control that allocates nothing fails the scaled check', () => {
    const step = allocating(20_000);
    const r = read(step);
    const bytes = scaledBytes(r);
    // a control that allocated nothing (a measurement blind to a leak) reads under three quarters
    const c = read(step);
    expect(c.bytesPerStep - r.bytesPerStep).toBeLessThan(0.75 * bytes);
  });
});
