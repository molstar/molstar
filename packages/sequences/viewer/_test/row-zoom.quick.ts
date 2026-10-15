/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of row-zoom.heap.ts (`test:perf`): the Alt wheel over rows and over a dense matrix, and the vertical
// overview's bottom grip, each on a viewer loaded afresh (row-zoom-cases.ts, the guard's cases), each with its control
// and no layout built while measured.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { altWheel, bottomGrip, matrixWheel } from './row-zoom-cases.js';

/**
 * The grip's cap: even under the quick flags its moves allocate 96 bytes each until about the 9,000th, then nothing
 * (2026-10-14, Node 24, M4); the warm-up runs on until it reads quiet.
 */
const GRIP = { warmSteps: 2_000, maxWarmSteps: 60_000 };

describe('the row zoom, quick', () => {
  it('the Alt wheel over the rows', async ({ task }) => {
    const { step, stats } = await altWheel();
    const builds = stats.builds;
    quickCase(task, 'row-zoom: alt wheel', step);
    expect(stats.builds).toBe(builds);
  });

  it('the Alt wheel over a dense matrix', async ({ task }) => {
    const { step, heights, stats } = await matrixWheel();
    expect(heights).toBe(2);
    const builds = stats.builds;
    quickCase(task, 'row-zoom: matrix wheel', step);
    expect(stats.builds).toBe(builds);
  });

  it("the vertical overview's bottom grip", async ({ task }) => {
    const { step, owner, stats, release } = await bottomGrip();
    expect(owner).toBe('verticalOverview');
    const builds = stats.builds;
    quickCase(task, 'row-zoom: bottom grip', step, GRIP);
    expect(stats.builds).toBe(builds);
    release();
  });
});
