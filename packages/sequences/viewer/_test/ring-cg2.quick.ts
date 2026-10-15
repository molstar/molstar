/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of ring-cg2.heap.ts (`test:perf`): CG2's ring frames, labels with leaders, arrowheads and chords, over
// 240 and 120,000 positions, each on a viewer loaded afresh (ring-cg2-cases.ts, the guard's case), each with its
// control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { CG2, ringCg2 } from './ring-cg2-cases.js';

describe("CG2's ring frames, quick", () => {
  it.for([240, 120_000])('labels with leaders, arrowheads and chords, %i positions', async (n, { task }) => {
    const { step, be } = await ringCg2(n);
    quickCase(task, `ring-cg2: ${n} positions`, step, CG2);
    expect(be.textItems.n).toBeGreaterThan(0);
    expect(be.pathItems.n).toBeGreaterThan(0);
  });
});
