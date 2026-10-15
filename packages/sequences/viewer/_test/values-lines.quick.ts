/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of values-lines.heap.ts (`test:perf`): the values frames through the conformance kit (values-lines-cases.ts,
// the guard's case): one measurement whose step cycles through the frames the kit's `heap` hook is handed at each zoom
// (kit-frames.ts), with its positive control; the guard's per-frame detail is the deep guard's.

import { expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { cycling, withLinearFrames } from './kit-frames.js';
import { valuesLines } from './values-lines-cases.js';

/** A frame of the line costs 50 to 200 us and the step cycles through a handful of them: windows of 250. */
const FRAME: QuickOptions = { warmSteps: 1_000, steps: 250, windows: 3, minWarmSteps: 1_500, maxWarmSteps: 6_000 };

it('a frame of values colored by their values, at every zoom', { timeout: 300_000 }, async ({ task }) => {
  const checks = await withLinearFrames(valuesLines, (frames) => {
    quickCase(task, 'values-lines: every zoom', cycling(frames), FRAME);
  });
  expect(checks.length).toBeGreaterThan(0);
  expect(checks.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
});
