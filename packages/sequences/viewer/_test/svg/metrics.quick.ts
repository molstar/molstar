/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of metrics.heap.ts (`test:perf`): `metrics` for a known font string, and `measure`, which allocates
// only the width it returns (V8 boxes it to a caller that does not inline it: 16 bytes, the guard's allowance), the
// guard's steps (metrics-cases.ts), each with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from '../heap.js';
import { measuring, metrics } from './metrics-cases.js';

it('metrics() once a font string is known', ({ task }) => {
  const { step, sum } = metrics();
  quickCase(task, 'metrics: metrics() once a font string is known', step);
  expect(sum()).toBeGreaterThan(0);
});

it('measure() allocates only the width it returns', ({ task }) => {
  const { step, sum } = measuring();
  quickCase(task, 'metrics: measure()', step);
  expect(sum()).toBeGreaterThan(0);
});
