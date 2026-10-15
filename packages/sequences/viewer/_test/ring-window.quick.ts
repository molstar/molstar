/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of ring-window.heap.ts (`test:perf`): ring frames with the window, the linear section moving and the
// pointer round the rings and over the labels, over 600 and 120,000 positions, each on a viewer loaded afresh
// (ring-window-heap-run.ts, the guards' cases), each with its control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { RING_WINDOW, setup } from './ring-window-heap-run.js';

describe("CG3's ring window and labels, quick", () => {
  it.for([600, 120_000])('ring frames with the window and the labels, %i positions', async (n, { task }) => {
    const s = await setup('frames', n);
    quickCase(task, `ring-window: frames, ${n} positions`, s.step, RING_WINDOW);
    expect(s.sectors()).toBeGreaterThan(0);
    expect(s.labelHits()).toBeGreaterThan(0);
  });
});
