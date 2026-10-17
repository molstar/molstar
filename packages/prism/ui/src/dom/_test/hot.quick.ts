// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of hot.heap.ts (`test:perf`): prism-ui's hot paths in jsdom, bytes a step against their baseline.

import { expect, it } from 'vitest';
import { allowance, measureBytesPerStep, quickCase } from '../../_test/perf.js';
import { commandMenuTyping, forUnchanged } from './hot-cases.js';

it('a For given the same rows in a new array', ({ task }) => {
  const { step, sink } = forUnchanged();
  quickCase(task, 'hot: For, the same rows anew', step);
  expect(sink()).toBe(200);
});

it('typing in the command menu', ({ task }) => {
  const { step, sink } = commandMenuTyping();
  const { r } = quickCase(task, 'hot: typing in the command menu', step, {
    steps: 200,
    warmSteps: 400,
    control: false,
    why: "jsdom allocates 1.2 MB a step here, its windows ±100 B apart: `leak`'s 88 B is under that; the control below is at the guard's own scale",
  });
  expect(sink()).toBeGreaterThan(0);
  // The positive control at the guard's scale: a step allocating a quarter more reads over the case's allowance (a tenth
  // of it), so a regression that size is caught.
  const extra = Math.ceil(r.bytesPerStep / 4 / 4);
  let kept: number[] = [];
  const c = measureBytesPerStep(
    (i) => {
      step(i);
      kept = new Array<number>(extra).fill(i & 1);
    },
    { steps: 200, warmSteps: 200, maxWarmSteps: 200 },
  );
  expect(kept.length).toBe(extra);
  expect(c.bytesPerStep, c.note).toBeGreaterThan(r.bytesPerStep + allowance(r.bytesPerStep, 4));
});
