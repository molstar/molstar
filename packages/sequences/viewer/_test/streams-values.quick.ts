/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-values.heap.ts (`test:perf`): a pan's pass and frame of the values draw beside the genes,
// variants and letters at 100 kb (streams-values-cases.ts, the guard's case), with its positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { valuesFrame } from './streams-values-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe("ST1's values draw, quick", () => {
  it('a pass and a frame, 100 kb at DPR 2, with genes, variants and letters', async ({ task }) => {
    const { g, frame } = await valuesFrame();
    frame();
    expect(g.backend.rectItems.n).toBeGreaterThan(500);
    const asked = g.asked();
    quickCase(task, 'streams-values: 100 kb', frame, FRAME);
    expect(g.asked()).toBe(asked);
    g.ctx.dispose();
  });
});
