/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of dot-plot-progressive.heap.ts (`test:perf`): the same case (dot-plot-progressive-cases.ts), the
// stacked lanes' frames past the budget over two pans, with its positive control.

import { describe, expect, it } from 'vitest';
import { DRAW, stacked } from './dot-plot-progressive-cases.js';
import { quickCase } from './heap.js';

describe("the dot plot's progressive draw, quick", () => {
  it('stacked lanes, tiles missing past the budget', async ({ task }) => {
    const g = await stacked();
    g.draw();
    const tiles = g.visits.tiles;
    g.draw();
    expect(g.images()).toBeGreaterThan(0);
    expect(g.t.composer.again).toBe(true);
    quickCase(task, 'dot-plot-progressive: stacked lanes', g.draw, DRAW);
    expect(g.visits.tiles).toBe(tiles);
    g.dispose();
  });
});
