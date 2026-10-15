/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-port.heap.ts (`test:perf`): a pan's pass, and its pass and composer frame, over the
// streamed genome at the whole view and at 10 kb, and with dynamic nodes held (streams-port-cases.ts, the guard's
// cases), each with its positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { PORT_CASES, portFrame } from './streams-port-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe("the streamed genome's pans, quick", () => {
  for (const { name, at, draw, dynamic } of PORT_CASES)
    it(name, async ({ task }) => {
      const g = await portFrame(at, draw, dynamic);
      if (dynamic) expect(g.manager.nodes.entries.size).toBeGreaterThan(at === 'whole' ? -1 : 0);
      const asked = g.requests();
      quickCase(task, `streams-port: ${name}`, g.frame, FRAME);
      expect(g.requests()).toBe(asked);
      if (!dynamic) expect(g.manager.stats.current.hits).toBeGreaterThan(0);
    });
});
