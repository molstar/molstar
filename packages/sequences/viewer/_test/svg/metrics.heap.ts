/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `metricsMeasurer` allocates nothing per call but the number it returns (layout prepares labels through it at edge rate, a
// long list at a time): `metrics` returns a record made once per font string; `measure` over strings of ASCII and of other
// code points makes no string, no array and no record, and V8 boxes the fractional width it returns to a caller that does
// not inline it (16 bytes), which is all the guard allows. The steps are metrics-cases.ts's, shared with the quick twin
// (metrics.quick.ts). The positive control: `metrics` calls allocating an object each (`leaking`), through the same
// measurement, are caught.

import { expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking, SCALED_LIMIT } from '../heap.js';
import { measuring, metrics, N } from './metrics-cases.js';

it('metrics() allocates nothing once a font string is known', () => {
  const { step, sum } = metrics();
  const r = heapScaling(N, step);
  expect(sum()).toBeGreaterThan(0);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
  expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
});

it('measure() allocates only the width it returns', () => {
  const { step, sum } = measuring();
  const r = heapScaling(N, step);
  expect(sum()).toBeGreaterThan(0);
  expect(r.scaled, describeScaling(r, N)).toBeLessThanOrEqual(3 * N * 16 + SCALED_LIMIT);
});
