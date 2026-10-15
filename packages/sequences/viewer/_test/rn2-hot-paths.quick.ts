/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of rn2-hot-paths.heap.ts (`test:perf`): the same cases (rn2-hot-paths-cases.ts). A kit case is one
// measurement whose step cycles through the linear frames the kit's `heap` hook is handed (every zoom, below the row's
// minimum and above; kit-frames.ts), with its positive control; the guard's per-frame detail, and the wrapped geometry,
// are the deep guard's. The pairs' hit and lit pairs measure one step.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { cycling, withLinearFrames } from './kit-frames.js';
import { CASES, caseName, kitRun, litPairsFrame, pairsHit } from './rn2-hot-paths-cases.js';

/** A frame costs 10 to 200 us and the step cycles through a handful of them: windows of 250. */
const FRAME: QuickOptions = { warmSteps: 500, steps: 250, windows: 3, minWarmSteps: 750 };

describe("RN2's hot paths, quick", () => {
  it.for(CASES)('$entry on $doc', { timeout: 300_000 }, async (c, { task }) => {
    const checks = (
      await withLinearFrames(
        (heap) => kitRun(c, heap),
        (frames) => {
          quickCase(task, `rn2-hot-paths: ${caseName(c)}`, cycling(frames), FRAME);
        },
      )
    ).filter((x) => x.check === 'heap');
    expect(checks.length).toBeGreaterThan(0);
    expect(checks.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
  });

  it('pairs: the hit over a grid of fractional points beside one long pair', async ({ task }) => {
    quickCase(task, 'rn2-hot-paths: pairs hit', await pairsHit(), FRAME);
  });

  it("pairs: the overlay's lit pairs through the channels' runs and the end index", async ({ task }) => {
    quickCase(task, 'rn2-hot-paths: pairs lit', await litPairsFrame(), FRAME);
  });
});
