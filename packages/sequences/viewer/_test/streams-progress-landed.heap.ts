/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3b's fade guard (seqspec:sources#levels-and-the-zoom-gate; the memory rule), its own file: every region of the
// streamed genome landed after a second of pending frames, the clock held so the fade never ends; an overlay frame
// allocates nothing. A frame allocating an object a step fails the same guard (the positive control). The case is
// streams-progress-cases.ts's, shared with the quick twin (streams-progress-landed.quick.ts); it loads the viewer afresh
// (`vi.resetModules`) and is warmed and settled by count (`frameScaling`).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { progressCase } from './streams-progress-cases.js';

const N = 1_000;

describe("the overlay's fade-in allocates nothing per frame", () => {
  it('every region landed, fading', async () => {
    const { g, frame, next, rectsOf, SkinToken } = await progressCase('landed');
    frame();
    // Not vacuous: the shade over every region just drawn; a frame asked while they fade.
    expect(rectsOf(g.backend, SkinToken.pending).n).toBeGreaterThan(2);
    expect(next[0]).toBeGreaterThan(0);
    const r = frameScaling(frame, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    g.ctx.dispose();
  }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const { g, frame } = await progressCase('landed');
    expect(controlCaught(() => frameScaling(leaking(frame), N))).toBe(true);
    g.ctx.dispose();
  }, 300_000);
});
