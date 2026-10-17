// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// prism-ui's hot paths keep nothing a step (`test:heap`): growth over n and 4n steps is flat, the control caught.

import { describe, expect, it } from 'vitest';
import {
  controlCaught,
  describeScaling,
  heapScaling,
  isFlat,
  isRetainedFlat,
  leaking,
  retainedScaling,
} from '../../_test/perf.js';
import { commandMenuTyping, forUnchanged } from './hot-cases.js';

const gc = (globalThis as { gc?: () => void }).gc!;
const settle = (step: (i: number) => void, n: number) => {
  gc();
  gc();
  for (let i = 0; i < n; i++) step(i);
};

describe('hot path heap guards', () => {
  it('a For given the same rows in a new array', () => {
    const { step } = forUnchanged();
    settle(step, 20_000);
    const r = heapScaling(2_000, step);
    expect(isFlat(r), describeScaling(r, 2_000)).toBe(true);
  });

  it('typing in the command menu', () => {
    // A keystroke allocates about 1.2 MB in jsdom by design (hot.quick.ts): what it keeps is measured after collections.
    const { step } = commandMenuTyping();
    settle(step, 1_000);
    const r = retainedScaling(200, step);
    expect(isRetainedFlat(r, 200), describeScaling(r, 200)).toBe(true);
  });

  it('a keystroke that keeps what it made is caught (positive control)', () => {
    const { step } = commandMenuTyping();
    settle(step, 1_000);
    const kept: unknown[] = [];
    // About 64 KB kept a keystroke: the 4n run keeps 51 MB more than the n run.
    const keeping = (i: number) => {
      step(i);
      kept.push(new Uint8Array(64 << 10));
    };
    expect(isRetainedFlat(retainedScaling(200, keeping), 200)).toBe(false);
  });

  it('the same step allocating an object each is caught (positive control)', () => {
    const { step } = forUnchanged();
    settle(step, 20_000);
    expect(controlCaught(() => heapScaling(2_000, leaking(step)))).toBe(true);
  });
});
