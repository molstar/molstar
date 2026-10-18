/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text view's lanes allocate nothing per frame or per pointer event after warm-up (seqspec:viewer#the-text-view;
// seqspec:viewer#performance): the lanes drawn whole each step as the rows scroll down and back (`draw`), and the
// pointer moving over the strips and the underlines through the real dispatcher with an overlay frame each (`hover`;
// the cases in wrapped-lanes-cases.ts). The positive control: the same steps allocating an object each (`leaking`),
// warmed and measured alike, are caught. Each case loads the viewer afresh; the draw runs first (the hover allocated
// after it only while a timed warm-up ended inside its tier-up: `FRAME_WARM_CAP`, heap.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { wrappedLanes } from './wrapped-lanes-cases.js';

/** A deep guard's steps: 500 catch about 43 bytes a step, these about 22. */
const STEPS = 1_000;

describe("the text view's lanes allocate nothing after warm-up", () => {
  it('a frame of letters, underlines, a values strip with gaps, a color strip and color beside the letters', async () => {
    const d = await wrappedLanes('draw');
    const r = frameScaling(d.step, STEPS);
    // The frame drew the strips' runs and the letters' cells.
    expect(d.state().rects).toBeGreaterThan(500);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 300_000);

  it('the pointer over the strips and the underlines, each step a hit and an overlay frame', async () => {
    const d = await wrappedLanes('hover');
    const r = frameScaling(d.step, STEPS);
    // A position hovered on most steps (the strips' layers', an underline's item, the letters').
    expect(d.state().hovered).toBeGreaterThan(4 * STEPS);
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 300_000);

  it('the same steps allocating an object each are caught (positive control)', async () => {
    const d = await wrappedLanes('hover');
    expect(controlCaught(() => frameScaling(leaking(d.step), STEPS))).toBe(true);
  }, 300_000);
});
