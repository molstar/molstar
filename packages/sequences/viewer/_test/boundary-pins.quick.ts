/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of boundary-pins.heap.ts (`test:perf`): the same hover sweep over the pinned boundaries
// (boundary-pins-cases.ts), a hit and an overlay frame a step, with its positive control.

import { describe, expect, it } from 'vitest';
import { boundaryPins } from './boundary-pins-cases.js';
import { quickCase } from './heap.js';

describe('pinned boundaries over a hover sweep, quick', () => {
  it('a hit and an overlay frame a move', async ({ task }) => {
    const { frame, boundary } = await boundaryPins();
    quickCase(task, 'boundary-pins: hover sweep', frame);
    // Not vacuous: the pins hovered.
    expect(boundary()).toBeGreaterThan(1_000);
  });
});
