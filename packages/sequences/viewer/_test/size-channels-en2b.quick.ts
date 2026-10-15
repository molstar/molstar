/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of size-channels-en2b.heap.ts (`test:perf`): its layers drawn through pans, a frame a pan and a whole
// redraw (the guard's case, size-channels-en2b-cases.ts), with its control.

import { describe, it } from 'vitest';
import { quickCase } from './heap.js';
import { en2bFrame } from './size-channels-en2b-cases.js';

/** A frame of several draws (about 0.5 ms) still read 2 to 6 bytes a step in some windows at 5,000 steps. */
const LONG = { warmSteps: 2_000, steps: 1_000, windows: 3, maxWarmSteps: 60_000 };

describe("frames of EN2b's size channels, quick", () => {
  it('pans', async ({ task }) => {
    quickCase(task, 'size-channels-en2b: pans', await en2bFrame(), LONG);
  });
});
