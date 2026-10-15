/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of size-channels.heap.ts (`test:perf`): its layers drawn through pans, a frame a pan and a whole
// redraw (the guard's case, size-channels-cases.ts), with its control.

import { describe, it } from 'vitest';
import { quickCase } from './heap.js';
import { sizeFrame } from './size-channels-cases.js';

/** A frame of several draws (about 0.5 ms) still read 1 to 4 bytes a step in some windows at 5,000 steps. */
const LONG = { warmSteps: 2_000, steps: 1_000, windows: 3, maxWarmSteps: 60_000 };

describe('frames of layers under size themes, quick', () => {
  it('pans', async ({ task }) => {
    quickCase(task, 'size-channels: pans', await sizeFrame(), LONG);
  });
});
