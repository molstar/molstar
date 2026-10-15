/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-miniature.heap.ts (`test:perf`): a pan's pass and composer frame with the `strip`
// navigator's map over the streamed genome at 200 kb (streams-miniature-cases.ts, the guard's case), with its positive
// control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { miniatureFrame } from './streams-miniature-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe('a miniature over streamed layers, quick', () => {
  it("a pass and a frame with the strip's map, 200 kb", async ({ task }) => {
    const g = await miniatureFrame();
    try {
      const p = g.pending();
      expect(p.maps).toBeGreaterThan(0);
      expect(p.n).toBe(0);
      const asked = g.requests();
      quickCase(task, 'streams-miniature: strip map 200 kb', g.frame, FRAME);
      expect(g.requests()).toBe(asked);
    } finally {
      g.dispose();
    }
  });
});
