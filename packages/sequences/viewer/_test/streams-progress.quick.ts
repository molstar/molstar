/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-progress.heap.ts (`test:perf`): the overlay's frame over the streamed genome with every
// region pending and with every region failed (streams-progress-cases.ts, the guard's cases), each with its positive
// control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { progressCase } from './streams-progress-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe("the overlay's skeletons and error bands, quick", () => {
  for (const c of ['pending', 'failed'] as const)
    it(`every region ${c}`, async ({ task }) => {
      const { g, frame, next, rectsOf, SkinToken } = await progressCase(c);
      frame();
      expect(rectsOf(g.backend, c === 'failed' ? SkinToken.conflict : SkinToken.pending).n).toBeGreaterThan(2);
      if (c === 'failed') expect(g.backend.pathItems.n).toBeGreaterThan(10);
      expect(next[0] > 0).toBe(c === 'pending');
      quickCase(task, `streams-progress: ${c}`, frame, FRAME);
      g.ctx.dispose();
    });
});
