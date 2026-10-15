/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3b's values guard (seqspec:sources#the-viewers-side; the memory rule): ST1's values draw over held regions beside
// the genes, variants and letters, a pan's pass and frame at 100 kb and DPR 2, allocates nothing. ST3a-2's guards leave
// the values track out: with it they measured about 23 bytes a frame (here 31, the region's values and stride read per
// bin, a load gone polymorphic beside the items' regions; `Layer.bins` now reads them once a region). ST1's
// streams-port.heap.ts guards values alone. A frame allocating an object a step fails the same guard (the positive
// control). The case is streams-values-cases.ts's, shared with the quick twin (streams-values.quick.ts); it loads the
// viewer afresh (`vi.resetModules`) and is warmed and settled by count (`frameScaling`): the genes and variants beside
// the values are drawn in color, which tiers up slowly.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { valuesFrame } from './streams-values-cases.js';

const N = 1_000;

describe("ST1's values draw allocates nothing per frame over held regions", () => {
  it('a pass and a frame, 100 kb at DPR 2, with genes, variants and letters', async () => {
    const { g, frame } = await valuesFrame();
    frame();
    // Not vacuous: the bars and the genes drawn.
    expect(g.backend.rectItems.n).toBeGreaterThan(500);
    const asked = g.asked();
    const r = frameScaling(frame, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(g.asked()).toBe(asked);
    g.ctx.dispose();
  }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const { g, frame } = await valuesFrame();
    expect(controlCaught(() => frameScaling(leaking(frame), N))).toBe(true);
    g.ctx.dispose();
  }, 300_000);
});
