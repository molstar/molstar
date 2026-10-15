/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3a-2's heap guards (seqspec:sources#the-viewers-side; the memory rule), regions held: a pan's pass and composer
// frame allocate nothing at 100 kb (genes as items region by region, variants as density), at 1 Mb (both as density
// from `count` bins) and at 200 Mb (the letters gate's hatch and hint). A frame allocating an object a step fails the same
// guard (the positive control). The glyphs and the placeholders are streams-items-glyphs.heap.ts's; the cases are
// streams-items-cases.ts's, shared with the quick twin (streams-items.quick.ts). Warmed and settled by count
// (`frameScaling`): streamed items drawn in color (a field theme's default reaches them) climb to TurboFan slowly, and
// the settle runs until the frame stops allocating.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { ITEMS_VIEWS, itemsCase } from './streams-items-cases.js';

const N = 1_000;

describe("the streamed items' pans allocate nothing per frame once their regions are held", () => {
  for (const [label, start, stop] of ITEMS_VIEWS)
    it(`a pass and a frame, ${label}`, async () => {
      const g = await itemsCase(start, stop);
      g.frame();
      // Not vacuous: genes and variants drawn (rects), at 200 Mb the hatch too.
      expect(g.backend.rectItems.n).toBeGreaterThan(100);
      if (start === 20_000_000) expect(g.backend.pathItems.n).toBeGreaterThan(0);
      const asked = g.asked();
      const r = frameScaling(g.frame, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(g.asked()).toBe(asked);
      g.dispose();
    }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await itemsCase(100_000_000, 100_100_000);
    expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
    g.dispose();
  }, 300_000);
});
