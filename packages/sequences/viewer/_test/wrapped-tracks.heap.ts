/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Tracks drawn in full in the text view allocate nothing per frame after warm-up (seqspec:viewer#the-text-view;
// seqspec:viewer#performance): a tall values line and an area, sites as glyphs, domains with their labels and
// secondary structure over every line, with the axes and the labels on, drawn whole each step as the rows scroll down
// and back. The positive control: the same steps allocating an object each (`leaking`), warmed and measured alike, are
// caught. The guard loads the viewer afresh (`vi.resetModules`) and warms by count (`warmFrame`); its case is
// wrapped-tracks-cases.ts's, shared with the quick twin (wrapped-tracks.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { wrappedTracks } from './wrapped-tracks-cases.js';

/** A deep guard's steps: 500 catch about 43 bytes a step, these about 22. */
const STEPS = 1_000;

/** The guard's measurement of a step: warmed and settled by count, counted anew, then `STEPS` and `4 * STEPS`. */
const measure = (d: Awaited<ReturnType<typeof wrappedTracks>>, step: (i: number) => void) => {
  warmFrame(step);
  d.resetCounts();
  return heapScaling(STEPS, step);
};

describe('tracks drawn in full in the text allocate nothing after warm-up', () => {
  it('a frame of a values line, an area, glyphs, labeled ranges and secondary structure, axes and labels on', async () => {
    const d = await wrappedTracks('all');
    const r = measure(d, d.step);
    // The frames drew the lines' paths, the glyphs and the axes' text.
    const { paths, glyphs, texts } = d.state();
    expect(paths).toBeGreaterThan(STEPS);
    expect(glyphs).toBeGreaterThan(0);
    expect(texts).toBeGreaterThan(0);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 300_000);

  it('the same frames allocating an object each are caught (positive control)', async () => {
    const d = await wrappedTracks('all');
    expect(controlCaught(() => measure(d, leaking(d.step)))).toBe(true);
  }, 300_000);
});
