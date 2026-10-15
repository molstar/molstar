/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of overlay-paths.heap.ts (`test:perf`): the vertical overview's scroll, the keyboard cursor's ring,
// a tree clade's hover and the two held drags' pending shapes (the guard's cases, overlay-paths-cases.ts), each with
// its control.

import { describe, it } from 'vitest';
import { quickCase } from './heap.js';
import { clade, overview, plotDrag, ring, selectDrag } from './overlay-paths-cases.js';

/** The hover and drag steps (about 20 microseconds) still read 1 to 24 bytes a step in their windows at 6,000 steps. */
const LONG = { warmSteps: 2_000, steps: 1_000, windows: 3, maxWarmSteps: 60_000 };

describe('overlay paths, quick', () => {
  it("the vertical overview and scrollbar's overlay through a vertical scroll", async ({ task }) => {
    quickCase(task, 'overlay-paths: overview', (await overview()).step);
  });

  it("the keyboard cursor's ring through a pan", async ({ task }) => {
    quickCase(task, 'overlay-paths: ring', (await ring()).step, LONG);
  });

  it("a tree clade's hover, on and off", async ({ task }) => {
    quickCase(task, 'overlay-paths: clade', (await clade()).step, LONG);
  });

  it("a held select drag's pending band on an alignment's rows", async ({ task }) => {
    const c = await selectDrag();
    quickCase(task, 'overlay-paths: select drag', c.step, LONG);
    c.held?.();
  });

  it("a held drag's pending rectangle on a dot plot", async ({ task }) => {
    const c = await plotDrag();
    quickCase(task, 'overlay-paths: dot plot drag', c.step, LONG);
    c.held?.();
  });
});
