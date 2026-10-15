/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of rn2-draws.heap.ts (`test:perf`): the same cases (rn2-draws-cases.ts). A kit case is one measurement
// whose step cycles through the linear frames the kit's `heap` hook is handed (every zoom, below the row's minimum and
// above, every layer; kit-frames.ts), with its positive control; the guard's per-frame detail, and the wrapped and ring
// geometries, are the deep guard's. The pairs overlay's modes measure one frame each.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { cycling, withLinearFrames } from './kit-frames.js';
import { CASES, kitRun, PAIRS_MODES, pairsOverlay } from './rn2-draws-cases.js';

/** A frame costs 10 to 200 us and the step cycles through up to two dozen of them: windows of 250. */
const FRAME: QuickOptions = { warmSteps: 500, steps: 250, windows: 3, minWarmSteps: 750 };

describe("RN2's feature and logo frames, quick", () => {
  it.for(CASES)('$entry on $doc', { timeout: 300_000 }, async (c, { task }) => {
    const checks = (
      await withLinearFrames(
        (heap) => kitRun(c, heap),
        (frames) => {
          quickCase(task, `rn2-draws: ${c.entry} on ${c.doc}`, cycling(frames), FRAME);
        },
      )
    ).filter((x) => x.check === 'heap');
    expect(checks.length).toBeGreaterThan(0);
    expect(checks.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
  });
});

describe('the pairs overlay, quick', () => {
  it.for(PAIRS_MODES)('arcs %s, many lit ends', async (mode, { task }) => {
    const { frame, paths } = await pairsOverlay(mode);
    expect(paths()).toBeGreaterThan(0);
    quickCase(task, `rn2-draws: pairs overlay: arcs ${mode}`, frame, FRAME);
  });
});
