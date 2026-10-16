/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An alignment drawn as text with frames allocates nothing per frame after warm-up (seqspec:viewer#rows-in-text,
// seqspec:viewer#frames, seqspec:viewer#performance): the globins' rows on their `espript` cells, their labels and
// numbers, a box on a row, two frames (one across a line break) with their titles, helices as cylinders with their
// labels and glyph markers, redrawn each step as the text scrolls down and back. The positive control: the same steps
// allocating an object each (`leaking`), warmed and measured alike, are caught. The guard loads the viewer afresh
// (`vi.resetModules`) and warms by count (`warmFrame`); its case is text-rows-cases.ts's, shared with the quick twin
// (text-rows.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { textRows } from './text-rows-cases.js';

/** A deep guard's steps: 500 catch about 43 bytes a step, these about 22. */
const STEPS = 1_000;

/** The guard's measurement of a step: warmed and settled by count, counted anew, then `STEPS` and `4 * STEPS`. */
const measure = (d: Awaited<ReturnType<typeof textRows>>, step: (i: number) => void) => {
  warmFrame(step);
  d.resetCounts();
  return heapScaling(STEPS, step);
};

describe('an alignment as text with frames allocates nothing after warm-up', () => {
  it('a frame of the rows on their cells, a box, frames and their titles, cylinders and markers, scrolled', async () => {
    const d = await textRows('draw');
    const r = measure(d, d.step);
    // The frames drew the outlines and cylinders (paths) and the letters, labels and numbers (texts).
    const { paths, texts } = d.state();
    expect(paths).toBeGreaterThan(STEPS);
    expect(texts).toBeGreaterThan(STEPS);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 300_000);

  it('the same frames allocating an object each are caught (positive control)', async () => {
    const d = await textRows('draw');
    expect(controlCaught(() => measure(d, leaking(d.step)))).toBe(true);
  }, 300_000);
});
