/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the stream manager's region cache (LY4; seqspec:viewer#performance), on `heapScaling` (a step is one
// frame's reads): a pass over a cached view with the view moving among held regions, on case 9's 18,000 lazy rows (a
// remote source) and on an alignment in memory (local tiles); and what a draw reads per frame of a layer: its regions
// and their states, the lanes in view and each block's levels. The positive control: case 9's pass allocating an object
// a step too (`leaking`), measured alike, is caught. The cases are streams-cases.ts's, shared with the quick twin
// (streams.quick.ts). Untimed: `heapScaling` runs `N` steps as its own warm-up, after the cases loaded every move's view.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { case9Pass, drawReads, memoryPass } from './streams-cases.js';

const N = 20_000;

describe("the region cache's per-frame work allocates nothing", () => {
  it("a pass over case 9's cached view (remote), the view moving among held regions", async () => {
    const { s, step } = await case9Pass();
    const requests = s.stats.current.requests;
    expect(s.stats.current.hits).toBeGreaterThan(0);
    const r = heapScaling(N, step);
    expect(s.stats.current.requests).toBe(requests);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('a pass over an alignment in memory zoomed out (local tiles held), the view moving', async () => {
    const { s, step, computed } = await memoryPass();
    expect(computed).toBeGreaterThan(0);
    const r = heapScaling(N, step);
    expect(s.stats.current.computed).toBe(computed);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("a draw's reads: the lanes in view, a block's levels, its regions and their states", async () => {
    const { step } = await drawReads();
    const r = heapScaling(N, step);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('the same passes allocating an object each are caught (positive control)', async () => {
    const { step } = await case9Pass();
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });
});
