/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of layout.heap.ts (`test:perf`): LY1's hot paths, hit testing, layout lookups, fits and labels,
// `layoutRuns` on resize and a wrapped line's lanes (layout-cases.ts, the guard's cases), each with its control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { fits, hitPoints, laneHeights, lookups, runsOnResize } from './layout-cases.js';

describe("LY1's hot paths, quick", () => {
  it('hit testing per pointer move', async ({ task }) => {
    const { step, kinds, nodes } = await hitPoints();
    expect([kinds > 10, nodes > 10]).toEqual([true, true]);
    quickCase(task, 'layout: hit testing', step);
  }, 120_000);

  it('layout lookups per frame, the view moving', async ({ task }) => {
    quickCase(task, 'layout: lookups', await lookups());
  });

  it('fits, numbers and cached labels', async ({ task }) => {
    quickCase(task, 'layout: fits and labels', await fits());
  });

  it('layoutRuns on resize', async ({ task }) => {
    const { step, runs } = await runsOnResize();
    expect(runs).toEqual([1, 35_000]);
    quickCase(task, 'layout: layoutRuns', step);
  });

  it("a wrapped line's lanes", async ({ task }) => {
    const { step, lines, segments } = await laneHeights();
    expect([lines > 200, segments > 1]).toEqual([true, true]);
    quickCase(task, 'layout: lane heights', step);
  });
});
