/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-items.heap.ts (`test:perf`): a pan's pass and composer frame over the streamed genes and
// variants at 100 kb, 1 Mb and 200 Mb (streams-items-cases.ts, the guard's cases), each with its positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { ITEMS_VIEWS, itemsCase } from './streams-items-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe("the streamed items' pans, quick", () => {
  for (const [label, start, stop] of ITEMS_VIEWS)
    it(`a pass and a frame, ${label}`, async ({ task }) => {
      const g = await itemsCase(start, stop);
      g.frame();
      expect(g.backend.rectItems.n).toBeGreaterThan(100);
      if (start === 20_000_000) expect(g.backend.pathItems.n).toBeGreaterThan(0);
      const asked = g.asked();
      quickCase(task, `streams-items: ${label}`, g.frame, FRAME);
      expect(g.asked()).toBe(asked);
      g.dispose();
    });
});
