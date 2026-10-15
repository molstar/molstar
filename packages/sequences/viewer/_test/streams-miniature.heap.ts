/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A miniature over streamed 1-D layers allocates nothing per pan once its regions are held (GM-04;
// seqspec:sources#the-viewers-side): over the streamed genome with the `strip` navigator, a pan's pass plans the map's
// regions every pass (each 1-D layer's level over its block's whole range) and the frame draws the tracks, the overlay
// and the navigator's map, whose layers read the map's level. Not vacuous: the map is drawn from its regions (no
// placeholder band in it), and nothing is asked while measured. The positive control: the same frame allocating an
// object a step (`leaking`), warmed, settled and measured alike, is caught. The case is streams-miniature-cases.ts's,
// shared with the quick twin (streams-miniature.quick.ts); it loads the viewer afresh (`vi.resetModules`) and is warmed
// and settled by count (`frameScaling`).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { miniatureFrame } from './streams-miniature-cases.js';

/** 1,000 steps, as the other stream guards: they catch about 22 bytes a step. */
const N = 1_000;

describe('a miniature over streamed layers allocates nothing per pan once its regions are held (GM-04)', () => {
  it("a pass and a frame with the strip's map, 200 kb", async () => {
    const g = await miniatureFrame();
    try {
      const p = g.pending();
      expect(p.maps).toBeGreaterThan(0);
      expect(p.n).toBe(0);
      const asked = g.requests();
      const r = frameScaling(g.frame, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(g.requests()).toBe(asked);
    } finally {
      g.dispose();
    }
  }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await miniatureFrame();
    try {
      expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
    } finally {
      g.dispose();
    }
  }, 300_000);
});
