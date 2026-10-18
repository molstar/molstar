/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-items-folded.heap.ts (`test:perf`): a pan's pass and composer frame over the folded genes
// at TP53 (streams-items-folded-cases.ts, the guard's case), with its positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { foldedFrame } from './streams-items-folded-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe("the folded genes' pans, quick", () => {
  it('a pass and a frame, TP53 folded', async ({ task }) => {
    const g = await foldedFrame();
    g.frame();
    const names = g.backend.textItems.text.slice(0, g.backend.textItems.n);
    expect(names.filter((x) => x === 'TP53').length).toBe(1);
    const asked = g.requests();
    quickCase(task, 'streams-items-folded: TP53', g.frame, FRAME);
    expect(g.requests()).toBe(asked);
    g.dispose();
  });
});
