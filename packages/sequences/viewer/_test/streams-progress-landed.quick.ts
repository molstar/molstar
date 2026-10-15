/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-progress-landed.heap.ts (`test:perf`): the overlay's frame while every region of the
// streamed genome fades in (streams-progress-cases.ts, the guard's case), with its positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { progressCase } from './streams-progress-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe("the overlay's fade-in, quick", () => {
  it('every region landed, fading', async ({ task }) => {
    const { g, frame, next, rectsOf, SkinToken } = await progressCase('landed');
    frame();
    expect(rectsOf(g.backend, SkinToken.pending).n).toBeGreaterThan(2);
    expect(next[0]).toBeGreaterThan(0);
    quickCase(task, 'streams-progress-landed: fading', frame, FRAME);
    g.ctx.dispose();
  });
});
