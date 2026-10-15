/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Line styles allocate nothing per frame after warm-up (seqspec:format#themes; seqspec:viewer#performance): a ranges
// layer in every style drawn whole each step, in tracks as the view pans and in the text as the rows scroll (the
// cases in line-styles-cases.ts: dashed and wavy strokes, doubles, boxes, highlights, and the text's styles behind its
// letters). The positive control: the same steps allocating an object each (`leaking`), warmed and measured alike, are
// caught. Each case loads the viewer afresh.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { lineStyles } from './line-styles-cases.js';

/** A deep guard's steps: 500 catch about 43 bytes a step, these about 22. */
const STEPS = 1_000;

describe('line styles allocate nothing after warm-up', () => {
  it('tracks: every style panned left and right, the section redrawn whole each step', async () => {
    const d = await lineStyles('tracks');
    const r = frameScaling(d.step, STEPS);
    // The frame drew the styles' strokes (dashed, dotted, waves) and rectangles.
    expect(d.state().paths).toBeGreaterThan(5);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 300_000);

  it('text: the styled underline lanes and the styles behind the letters, scrolled down and back', async () => {
    const d = await lineStyles('text');
    const r = frameScaling(d.step, STEPS);
    expect(d.state().paths).toBeGreaterThan(10);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 300_000);

  it('the same steps allocating an object each are caught (positive control)', async () => {
    const d = await lineStyles('tracks');
    expect(controlCaught(() => frameScaling(leaking(d.step), STEPS))).toBe(true);
  }, 300_000);
});
