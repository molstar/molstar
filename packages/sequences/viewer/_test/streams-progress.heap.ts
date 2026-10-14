/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3b's heap guards (seqspec:sources#levels-and-the-zoom-gate; the memory rule): the overlay's progressive loading
// allocates nothing per frame. A live composer over the streamed genome (genes, variants, letters, a signal) at 100 kb,
// DPR 2: every region pending (skeletons shimmering, the clock advancing a frame each step, and the surface's next
// frame decided from `animating`, as `SectionSurface` does), and every region failed (error bands, hatches and hints).
// Each case on fresh modules; the fade is streams-progress-landed.heap.ts's (run after these in one isolate, it
// allocated where alone it did not, as ST3a-2's glyph case). A frame allocating an object a step fails the same guard (the
// positive control).

import { describe, expect, it, vi } from 'vitest';
import {
  controlCaught,
  describeScaling,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from './heap.js';

const N = 1_000;
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame();
  // Settled after the timed warm-up (GS-02): a climb to TurboFan can outlast it.
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? '' : describeScaling(r, N);
};

describe("the overlay's skeletons and error bands allocate nothing per frame", () => {
  for (const c of ['pending', 'failed'] as const)
    it(`every region ${c}`, async () => {
      vi.resetModules();
      const { overlayFrame, rectsOf, SkinToken } = await import('./stream-progress-fixtures.js');
      const { g, frame, next } = await overlayFrame(c);
      frame();
      // Not vacuous: the shades drawn (and the hatches over failed regions); a frame asked while one shimmers.
      expect(rectsOf(g.backend, c === 'failed' ? SkinToken.conflict : SkinToken.pending).n).toBeGreaterThan(2);
      if (c === 'failed') expect(g.backend.pathItems.n).toBeGreaterThan(10);
      expect(next[0] > 0).toBe(c === 'pending');
      expect(guard(frame)).toBe('');
      g.ctx.dispose();
    }, 120_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    vi.resetModules();
    const { overlayFrame } = await import('./stream-progress-fixtures.js');
    const { frame } = await overlayFrame('pending');
    expect(controlCaught(() => measure(leaking(frame) as () => void))).toBe(true);
  }, 120_000);
});
