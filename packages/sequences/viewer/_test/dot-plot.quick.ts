/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of dot-plot.heap.ts (`test:perf`): the same cases (dot-plot-cases.ts), the draw over two views, the
// draw zoomed out over two pans and the pointer over its cells, each on a viewer loaded afresh with its positive
// control.

import { describe, expect, it } from 'vitest';
import { DRAW, plot, zoomedOut } from './dot-plot-cases.js';
import { quickCase } from './heap.js';

describe('the dot plot, quick', () => {
  it('its draw over two views', async ({ task }) => {
    const g = await plot();
    g.draw();
    g.draw();
    expect(g.be.pathItems.n).toBeGreaterThan(10);
    expect(g.images()).toBeGreaterThan(0);
    quickCase(task, 'dot-plot: draw', g.draw, DRAW);
    g.dispose();
  });

  it('its draw zoomed out over two pans', async ({ task }) => {
    const g = await zoomedOut();
    g.draw();
    g.draw();
    expect(g.images()).toBeGreaterThan(0);
    quickCase(task, 'dot-plot: draw zoomed out', g.draw, DRAW);
    g.dispose();
  });

  it('the pointer over its cells', async ({ task }) => {
    const g = await plot();
    g.hover();
    expect(g.t.ctx.hover.live.current.pairIndex).toBeGreaterThanOrEqual(0);
    quickCase(task, 'dot-plot: hover', g.hover);
    g.dispose();
  });
});
