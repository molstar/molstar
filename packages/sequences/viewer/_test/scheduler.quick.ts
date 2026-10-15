/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of scheduler.heap.ts (`test:perf`): a frame and a timer with stable callbacks through `ctx.scheduler`
// (scheduler-cases.ts, the guard's case), with its control.

import { expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { looping } from './scheduler-cases.js';

it('a frame loop and a timer with stable callbacks through ctx.scheduler, quick', ({ task }) => {
  const l = looping();
  quickCase(task, 'scheduler: frame and timer', l.step);
  expect(l.frames).toBe(0);
  l.ctx.dispose();
});
