/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of ring.heap.ts (`test:perf`): whole ring frames with polar hovers, over 240 and 120,000 positions,
// each on a viewer loaded afresh (ring-cases.ts, the guard's case), each with its control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { RING, RING_LARGE, ring } from './ring-cases.js';

describe('ring frames, quick', () => {
  it.for([240, 120_000])('whole frames with polar hovers round the rings, %i positions', async (n, { task }) => {
    const { step, be } = await ring(n);
    quickCase(task, `ring: ${n} positions`, step, n > 1_000 ? RING_LARGE : RING);
    expect(be.sectorItems.n).toBeGreaterThan(0);
  });
});
