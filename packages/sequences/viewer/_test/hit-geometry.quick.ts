/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of hit-geometry.heap.ts (`test:perf`): the same hover sweeps (hit-geometry-cases.ts), over the sized
// shapes and over an own unit's row in "residue" mode, a hit and an overlay frame a step, each with its positive
// control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { ownRow, SHAPES, sizedShapes } from './hit-geometry-cases.js';

describe('hover frames over sized shapes, quick', () => {
  it('dots, bands, boundaries, pair ends and arcs', async ({ task }) => {
    const { frame, hovers } = await sizedShapes();
    quickCase(task, 'hit-geometry: sized shapes', frame, SHAPES);
    // Not vacuous: the hover was rebuilt.
    expect(hovers()).toBeGreaterThan(1_000);
  });

  it('"residue" mode over an own unit\'s row', async ({ task }) => {
    const { frame } = await ownRow();
    quickCase(task, 'hit-geometry: own row residue', frame);
  });
});
