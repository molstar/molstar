/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of recording.heap.ts (`test:perf`): a recorded frame of 4,000 rects and 1,000 text runs, the guard's
// step (recording-cases.ts), with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from '../heap.js';
import { FRAMES, recorded } from './recording-cases.js';

it('recorded frames', ({ task }) => {
  const { step, rects } = recorded();
  quickCase(task, 'recording: a recorded frame', step, FRAMES);
  expect(rects()).toBe(4001);
});
