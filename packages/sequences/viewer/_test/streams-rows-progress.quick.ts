/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-rows-progress.heap.ts (`test:perf`): the overlay's frame over a lazy alignment's rows with
// every tile pending and with every tile failed (streams-rows-progress-cases.ts, the guard's cases), each with its
// positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { rowsFrame } from './streams-rows-progress-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe('the 2-D skeleton and error bands, quick', () => {
  for (const c of ['pending', 'failed'] as const)
    it(`every tile ${c}`, async ({ task }) => {
      const { t, frame, drawn } = await rowsFrame(c);
      frame();
      expect(drawn()).toBeGreaterThan(0);
      expect(t.composer.animating).toBe(c === 'pending');
      quickCase(task, `streams-rows-progress: ${c}`, frame, FRAME);
      t.ctx.dispose();
    });
});
