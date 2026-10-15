/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The 2-D skeleton allocates nothing per frame (seqspec:sources#levels-and-the-zoom-gate; seqspec:viewer#performance):
// `progressRows`, the overlay over a lazy alignment's rows (3,000 × 400 through the tiles protocol), live, DPR 2, with
// every tile pending (skeletons shimmering, the clock advancing a frame each step) and with every tile failed (error
// bands). The positive control: the pending frames allocating an object each (`leaking`), warmed and measured alike,
// are caught. The cases are streams-rows-progress-cases.ts's, shared with the quick twin
// (streams-rows-progress.quick.ts); each loads the viewer afresh (`vi.resetModules`) and is warmed and settled by count
// (`frameScaling`).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { rowsFrame } from './streams-rows-progress-cases.js';

const N = 1_000;
/**
 * The failed frame is cheap (about 25 microseconds) and tiers up late: with `frameScaling`'s warm-up alone (2,000 frames,
 * then a settle that ends on two windows allocating nothing, about 4,000 in all) the climb to the top tier fell inside
 * the measured run, which the heap growth counts (37 bytes a step, 2026-10-15) and the settle's exact reading of the
 * objects did not see. These frames run this many more by count first.
 */
const TIER_UP = 40_000;
/** The file's warm-up and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  for (let i = 0; i < TIER_UP; i++) frame();
  return frameScaling(frame, N);
};

describe('the 2-D skeleton and error bands allocate nothing per frame', () => {
  for (const c of ['pending', 'failed'] as const)
    it(`every tile ${c}`, async () => {
      const { t, frame, drawn } = await rowsFrame(c);
      frame();
      // Not vacuous: the shades drawn; a pending skeleton shimmers.
      expect(drawn()).toBeGreaterThan(0);
      expect(t.composer.animating).toBe(c === 'pending');
      const r = measure(frame);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      t.ctx.dispose();
    }, 300_000);

  it('the same pending frames allocating an object each are caught (positive control)', async () => {
    const { t, frame } = await rowsFrame('pending');
    expect(controlCaught(() => measure(leaking(frame) as () => void))).toBe(true);
    t.ctx.dispose();
  }, 300_000);
});
