/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot's heap guard (the memory rule; seqspec:viewer#performance): its draw (grid tiles at level 0 over two
// views: the fit, about 3 device px a cell, and zoomed in; the repeat's segments), its draw zoomed out past a device
// pixel a cell (2,400 residues at 0.45 px: level 2 or past it, the pyramid's bins, over two pans), and the pointer
// swept over its cells with an overlay frame each (the partner's hit, both positions lit, the rows and the hovered
// cell) allocate nothing once warm (the cases in dot-plot-cases.ts, shared with the quick twin). Each case loads the
// viewer afresh, warmed and settled by count (`frameScaling`). The positive control: the draw's frames allocating an
// object each (`leaking`), through the same warm-up and measurement, are caught.

import { describe, expect, it } from 'vitest';
import { plot, zoomedOut } from './dot-plot-cases.js';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';

/** 1,000 steps: they catch about 22 bytes a step, the trimmed 500 about 43 (heap.ts). */
const N = 1_000;
const guard = (frame: () => void) => {
  const r = frameScaling(frame, N);
  return isFlat(r) ? '' : describeScaling(r, N);
};

describe('the dot plot allocates nothing per frame', () => {
  it('its draw over two views: grid tiles at two levels and the segments', async () => {
    const g = await plot();
    g.draw();
    g.draw();
    // Not vacuous: the tiles and the segments are drawn.
    expect(g.be.pathItems.n).toBeGreaterThan(10);
    expect(g.images()).toBeGreaterThan(0);
    expect(guard(g.draw)).toBe('');
    g.dispose();
  }, 300_000);

  it("its draw zoomed out: past level 1, the pyramid's bins, over two pans", async () => {
    const g = await zoomedOut();
    g.draw();
    g.draw();
    expect(g.images()).toBeGreaterThan(0);
    expect(guard(g.draw)).toBe('');
    g.dispose();
  }, 300_000);

  it("the pointer over its cells: the partner's hit, both positions lit, the rows and the cell outlined", async () => {
    const g = await plot();
    g.hover();
    expect(g.t.ctx.hover.live.current.pairIndex).toBeGreaterThanOrEqual(0);
    expect(guard(g.hover)).toBe('');
    g.dispose();
  }, 300_000);

  it("the draw's frames allocating an object each are caught (positive control)", async () => {
    const g = await plot();
    expect(controlCaught(() => frameScaling(leaking(g.draw), N))).toBe(true);
    g.dispose();
  }, 300_000);
});
