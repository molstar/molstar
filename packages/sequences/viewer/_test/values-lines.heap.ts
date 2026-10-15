/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Values colored by their values through the conformance kit's heap guard (seqspec:format#layer-types;
// seqspec:viewer#the-render-backend): a frame of `linesDoc`'s line in four bands (a discrete palette, its color changing
// halfway between two positions), its area under a continuous palette and its plain line allocates nothing after
// warm-up, at each of the kit's zooms (one below a device px per position: columns' means and bands). The positive
// control goes through the kit's `heap` hook: at the first three of those frames the same warm-up and measurement of the
// frame allocating an object each (`leaking`) must catch it. The warm-up is by count (`warmFrame`). The case is
// values-lines-cases.ts's, shared with the quick twin (values-lines.quick.ts).
//
// The viewer is loaded afresh per case (`vi.resetModules`). The line's tiny frame (0.62 px a position) stays in V8's
// middle tier for tens of thousands of frames under the vitest module runner (draws.heap.ts's header): the warm-up is
// by count and settles it, so the case no longer depends on running alone.

import { expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { valuesLines } from './values-lines-cases.js';

/** The guard's measurement: warmed and settled by count (`warmFrame`), then `FRAME_N` frames and 4× measured. */
const measure = (frame: (i: number) => void) => {
  warmFrame(frame);
  return heapScaling(FRAME_N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, FRAME_N));
};
/**
 * The positive control as the kit's `heap` hook: a frame's check fails unless its leaking frame is caught. Its warm-up
 * never settles (the frame allocates once optimized: `FRAME_WARM_CAP` frames), so it runs on the first three frames
 * the kit hands over (the fit's two draws and the middle zoom's), not on all of them.
 */
const control = (frame: () => void, where: string) => {
  if (Number(where.split('#')[1]) > 3) return;
  if (!controlCaught(() => measure(leaking(frame)))) throw new Error('the leaking frame was not caught');
};

it('a frame of values colored by their values allocates nothing after warm-up, at every zoom', async () => {
  const heap = await valuesLines(guard);
  expect(heap.length).toBeGreaterThan(0);
  expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
}, 300_000);

it('a values frame allocating an object each is caught, at every zoom (positive control)', async () => {
  const heap = await valuesLines(control);
  expect(heap.length).toBeGreaterThan(0);
  expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
}, 300_000);
