/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of split-units.heap.ts (`test:perf`): the pointer moving over split rows on units of their own, an
// overlay frame a move (the guard's case, split-units-cases.ts), with its control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { unitsHover } from './split-units-cases.js';

describe('hover frames over split rows on their own units, quick', () => {
  it('a hover and an overlay frame a move', async ({ task }) => {
    const { frame, hovers } = await unitsHover();
    quickCase(task, 'split-units: hover frames', frame);
    expect(hovers()).toBeGreaterThan(1_000);
  });
});
