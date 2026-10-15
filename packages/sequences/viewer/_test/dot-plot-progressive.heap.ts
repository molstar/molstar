/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot's progressive draw allocates nothing per frame (the memory rule; seqspec:viewer#performance): stacked
// lanes on a progressive backend, each lane's own grid, its fine tiles never built (the budget spent: none a frame
// once the coarse levels are), so every frame takes the budget, finds tiles missing, draws the coarse level instead,
// asks for another frame (`more`, the composer's `again`) and draws the plot's entries again alone, over two pans (the
// case in dot-plot-progressive-cases.ts, shared with the quick twin). Warmed and settled by count (`frameScaling`). The
// positive control: the same frames allocating an object each (`leaking`), through the same warm-up and measurement,
// are caught.

import { describe, expect, it } from 'vitest';
import { stacked } from './dot-plot-progressive-cases.js';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';

const N = 1_000;
const guard = (frame: () => void) => {
  const r = frameScaling(frame, N);
  return isFlat(r) ? '' : describeScaling(r, N);
};

describe("the dot plot's progressive draw allocates nothing per frame", () => {
  it('stacked lanes, tiles missing past the budget: the coarse level drawn, another frame asked, the entry alone', async () => {
    const g = await stacked();
    g.draw();
    const tiles = g.visits.tiles;
    g.draw();
    // Not vacuous: tiles drawn (the coarse ones), more asked, nothing built.
    expect(g.images()).toBeGreaterThan(0);
    expect(g.t.composer.again).toBe(true);
    expect(guard(g.draw)).toBe('');
    expect(g.visits.tiles).toBe(tiles);
    g.dispose();
  }, 300_000);

  it('the same frames allocating an object each are caught (positive control)', async () => {
    const g = await stacked();
    expect(controlCaught(() => frameScaling(leaking(g.draw), N))).toBe(true);
    g.dispose();
  }, 300_000);
});
