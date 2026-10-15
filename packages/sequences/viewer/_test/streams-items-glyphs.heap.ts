/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3a-2's heap guards (seqspec:sources#the-viewers-side; the memory rule): at 5 kb, regions held, a pan's pass and
// composer frame allocate nothing (variants as glyphs region by region, stacked or as points, genes as items, letters as
// cells); with every region pending, a pan's frame draws the placeholder bands and allocates nothing. A frame
// allocating an object a step fails the same guard (the positive control). The cases are streams-items-glyphs-cases.ts's,
// shared with the quick twin (streams-items-glyphs.quick.ts); each loads the viewer afresh and is warmed and settled
// by count (`frameScaling`).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { glyphsCase, pendingBands } from './streams-items-glyphs-cases.js';

const N = 1_000;

describe('the streamed glyphs and placeholders allocate nothing per frame', () => {
  it('a pass and a frame, 5 kb', async () => {
    const g = await glyphsCase('glyphs');
    g.frame();
    expect(g.backend.glyphItems.n).toBeGreaterThan(0);
    const asked = g.asked();
    const r = frameScaling(g.frame, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(g.asked()).toBe(asked);
    g.dispose();
  }, 300_000);

  it('a pass and a frame, 5 kb, the variants as points (unstacked glyphs, region by region)', async () => {
    const g = await glyphsCase('points');
    g.frame();
    expect(g.backend.glyphItems.n).toBeGreaterThan(0);
    const asked = g.asked();
    const r = frameScaling(g.frame, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(g.asked()).toBe(asked);
    g.dispose();
  }, 300_000);

  it('a frame with every region pending (placeholder bands)', async () => {
    const g = await glyphsCase('pending');
    g.frame();
    expect(await pendingBands(g)).toBeGreaterThan(2);
    const r = frameScaling(g.frame, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    g.dispose();
  }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await glyphsCase('glyphs');
    expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
    g.dispose();
  }, 300_000);
});
