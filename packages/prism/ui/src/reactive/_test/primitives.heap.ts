/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The base primitives keep nothing a step (`test:heap`): heap growth over n and 4n steps is flat; the positive control,
// the same step allocating an object each (`leaking`), is caught. The steps are primitives-cases.ts's.

import { setFlagsFromString } from 'node:v8';
import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../../_test/perf.js';
import {
  bindWrite,
  deriveRecompute,
  deriveSkip,
  frameOffers,
  keyedSelectorMove,
  keyedStateSet,
  stateSet,
  watchRender,
} from './primitives-cases.js';

// The collector on the main thread: a concurrent sweep finishing inside a measured run moved its reading by hundreds of
// kilobytes either way (2026-10-17), far over what a step of these primitives could allocate.
setFlagsFromString('--single-threaded-gc');

const N = 10_000;
const gc = (globalThis as { gc?: () => void }).gc!;

/**
 * Collected first (what earlier cases left would be collected inside the guard and deoptimize the code they share),
 * then warmed by count, so each step runs optimized: under V8's defaults a step can sit 15,000 to 60,000 steps in its
 * middle tier, allocating (packages/sequences/AGENTS.md).
 */
function settle(step: (i: number) => void) {
  gc();
  gc();
  for (let i = 0; i < Warm; i++) step(i);
}

/** Steps before measuring, so each runs optimized: under V8's defaults a step can sit 15,000 to 60,000 steps in its middle tier, allocating (packages/sequences/AGENTS.md). */
const Warm = 100_000;
const cases = {
  stateSet,
  deriveRecompute,
  deriveSkip,
  bindWrite,
  watchRender,
  keyedSelectorMove,
  keyedStateSet,
  frameOffers,
};

describe('primitives heap guards', () => {
  for (const [name, make] of Object.entries(cases)) {
    it(name, () => {
      const { step, sink } = make();
      settle(step);
      const r = heapScaling(N, step);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(sink()).not.toBe(0.5);
    });
  }

  it('the same step allocating an object each is caught (positive control)', () => {
    const { step } = stateSet();
    settle(step);
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });
});
