/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The folded genes' heap guard (the lane fold, ST4b follow-up): a pan's pass and composer frame over the demo's GENCODE
// slices at TP53 (about 30 transcripts folded into two lanes) allocate nothing once the regions are held, the fold's
// `fold(fit)` and the folded items' read included. Not vacuous: the folded items are what is drawn (two lanes, TP53
// labeled once), and nothing is asked while measured. A frame allocating an object a step fails the same guard (the
// positive control). The case is streams-items-folded-cases.ts's, shared with the quick twin
// (streams-items-folded.quick.ts); warmed and settled by count (`frameScaling`).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { foldedFrame } from './streams-items-folded-cases.js';

const N = 1_000;

describe("the folded genes' pans allocate nothing per frame once their regions are held", () => {
  it('a pass and a frame, TP53 folded', async () => {
    const g = await foldedFrame();
    g.frame();
    const names = g.backend.textItems.text.slice(0, g.backend.textItems.n);
    expect(names.filter((x) => x === 'TP53').length).toBe(1);
    const asked = g.requests();
    const r = frameScaling(g.frame, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(g.requests()).toBe(asked);
    g.dispose();
  }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await foldedFrame();
    expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
    g.dispose();
  }, 300_000);
});
