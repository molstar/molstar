/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams.heap.ts (`test:perf`): the region cache's passes and a draw's reads (streams-cases.ts, the
// guard's cases), each with its positive control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { case9Pass, drawReads, memoryPass } from './streams-cases.js';

describe("the region cache's per-frame work, quick", () => {
  it("a pass over case 9's cached view (remote)", async ({ task }) => {
    const { s, step } = await case9Pass();
    const requests = s.stats.current.requests;
    expect(s.stats.current.hits).toBeGreaterThan(0);
    quickCase(task, 'streams: case 9 pass', step);
    expect(s.stats.current.requests).toBe(requests);
  });

  it('a pass over an alignment in memory zoomed out (local tiles)', async ({ task }) => {
    const { s, step, computed } = await memoryPass();
    expect(computed).toBeGreaterThan(0);
    quickCase(task, 'streams: in-memory pass', step);
    expect(s.stats.current.computed).toBe(computed);
  });

  it("a draw's reads: the lanes in view, a block's levels, its regions and their states", async ({ task }) => {
    const { step } = await drawReads();
    quickCase(task, 'streams: draw reads', step);
  });
});
