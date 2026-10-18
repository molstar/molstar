/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of header-hover.heap.ts (`test:perf`): the hover alternating between two track headers and a group's
// header, described and projected, the guard's step (header-hover-cases.ts, the viewer opened afresh), with its
// positive control.

import { expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { headers } from './header-hover-cases.js';

it('headers alternating, each described and projected', async ({ task }) => {
  const h = await headers();
  const before = h.ch.rebuilds.hover;
  // A step takes about 0.5 ms (a hover change, the channels' rebuild): windows of 250 steps.
  quickCase(task, 'header-hover: headers alternating', h.step, {
    warmSteps: 1_000,
    steps: 250,
    windows: 3,
  });
  expect(h.counts).toEqual([40, 10, 5]);
  expect(h.ch.rebuilds.hover - before).toBeGreaterThanOrEqual(5 * 1_000);
  expect(h.sink).toBeGreaterThan(0);
  h.t.ctx.dispose();
});
