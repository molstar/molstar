/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-letters.heap.ts (`test:perf`): a pan's pass and composer frame over the streamed letters
// at 10 kb and at 100 positions (streams-letters-cases.ts, the guard's cases), each with its positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { LETTER_VIEWS, lettersFrame } from './streams-letters-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe("the streamed letters' pans, quick", () => {
  for (const at of LETTER_VIEWS)
    it(`a pass and a frame, ${at}`, async ({ task }) => {
      const g = await lettersFrame(at);
      g.frame();
      if (at === '100 b') expect(g.backend.textItems.n).toBeGreaterThan(50);
      else expect(g.backend.rectItems.n).toBeGreaterThan(500);
      const asked = g.requests();
      quickCase(task, `streams-letters: ${at}`, g.frame, FRAME);
      expect(g.requests()).toBe(asked);
    });
});
