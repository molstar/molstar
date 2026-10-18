/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The row zoom's heap guards (seqspec:viewer#performance; KS2): a row scale alone rescales the rows (`Layouts.rescale`:
// the rows' height and the entries below them), never rebuilding the layout, so the Alt wheel's zoom and the vertical
// overview's bottom grip are flat per event, with no build counted, and so is the Alt wheel over a dense matrix, whose
// track follows the row scale. Before KS2 each event rebuilt the layout (about 2,120 bytes). The positive control: the
// Alt wheel's calls allocating an object each (`leaking`), through the same warm-up (`prime`, by count) and
// measurement, are caught. Each case loads the viewer afresh (`vi.resetModules`; the cases in row-zoom-cases.ts, shared
// with the quick twin, row-zoom.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, HEAP_DEEP, heapScaling, isFlat, leaking } from './heap.js';
import { altWheel, bottomGrip, matrixWheel } from './row-zoom-cases.js';

const N = HEAP_DEEP ? 50_000 : 25_000;

/** As input.heap's: at least 300,000 steps, then until three chunks in a row allocate under 2 KB. */
function prime(step: (i: number) => void) {
  const gc = (globalThis as { gc?: () => void }).gc;
  let k = 0;
  for (let chunk = 0, quiet = 0; chunk < 100 && (quiet < 3 || chunk < 15); chunk++) {
    gc?.();
    const before = process.memoryUsage().heapUsed;
    for (let j = 0; j < 20_000; j++) step(k++);
    quiet = process.memoryUsage().heapUsed - before < 2048 ? quiet + 1 : 0;
  }
}

describe('the row zoom rescales the rows: no rebuild, flat', () => {
  it('the Alt wheel: two hits a call, no build; an object allocated a call is caught (positive control)', async () => {
    const { step, counter, stats } = await altWheel();
    prime(step);
    const [builds, rescales] = [stats.builds, stats.rescales];
    counter.hits = 0;
    const calls = 4000;
    for (let i = 0; i < calls; i++) step(i);
    // The wheel's own hit, and the pointer's again after the rescale: two a call, not three.
    expect(counter.hits / calls).toBeLessThanOrEqual(2);
    expect(counter.hits / calls).toBeGreaterThan(1);
    const r = heapScaling(N, step);
    expect(stats.builds).toBe(builds);
    expect(stats.rescales - rescales).toBeGreaterThan(calls + N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(
      controlCaught(() => {
        const leaky = leaking(step);
        prime(leaky);
        for (let i = 0; i < calls; i++) leaky(i);
        return heapScaling(N, leaky);
      }),
    ).toBe(true);
    expect(counter.sink).toBeGreaterThanOrEqual(0);
  }, 180_000);

  it('the Alt wheel over a dense matrix: its track rescaled each call, no build, flat', async () => {
    const { step, heights, stats, sink } = await matrixWheel();
    // The matrix's track (the first scrolled entry) follows the zoom.
    expect(heights).toBe(2);
    prime(step);
    const [builds, rescales] = [stats.builds, stats.rescales];
    const r = heapScaling(N, step);
    expect(stats.builds).toBe(builds);
    expect(stats.rescales - rescales).toBeGreaterThan(N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(sink()).toBeGreaterThanOrEqual(0);
  }, 180_000);

  it("the vertical overview's bottom grip: rescaled each move, no build, flat", async () => {
    const { step, idle, owner, stats, release } = await bottomGrip();
    expect(idle).toBe(true);
    expect(owner).toBe('verticalOverview');
    prime(step);
    const [builds, rescales] = [stats.builds, stats.rescales];
    const r = heapScaling(N, step);
    expect(stats.builds).toBe(builds);
    expect(stats.rescales - rescales).toBeGreaterThan(N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    release();
  }, 180_000);
});
