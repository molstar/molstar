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
// The fade is streams-progress-landed.heap.ts's. A frame allocating an object a step fails the same guard (the
// positive control). The cases are streams-progress-cases.ts's, shared with the quick twin (streams-progress.quick.ts);
// each loads the viewer afresh (`vi.resetModules`) and is warmed and settled by count (`frameScaling`).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { progressCase } from './streams-progress-cases.js';

const N = 1_000;

describe("the overlay's skeletons and error bands allocate nothing per frame", () => {
  for (const c of ['pending', 'failed'] as const)
    it(`every region ${c}`, async () => {
      const { g, frame, next, rectsOf, SkinToken } = await progressCase(c);
      frame();
      // Not vacuous: the shades drawn (and the hatches over failed regions); a frame asked while one shimmers.
      expect(rectsOf(g.backend, c === 'failed' ? SkinToken.conflict : SkinToken.pending).n).toBeGreaterThan(2);
      if (c === 'failed') expect(g.backend.pathItems.n).toBeGreaterThan(10);
      expect(next[0] > 0).toBe(c === 'pending');
      const r = frameScaling(frame, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      g.ctx.dispose();
    }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const { g, frame } = await progressCase('pending');
    expect(controlCaught(() => frameScaling(leaking(frame), N))).toBe(true);
    g.ctx.dispose();
  }, 300_000);
});
