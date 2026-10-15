/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// XR1b's heap guard through the conformance kit (seqspec:viewer#extension-points): a frame of `org.example:stars` (a
// star per item through `mark`, colored by a field provider's value) allocates nothing after warm-up, at each of the
// kit's zooms and on rings (CG1b: drawn through the ring map); and RN3's entries, the rows and `seqv:map` (the cases in
// conformance-cases.ts, shared with the quick twin). Each frame is warmed by count and settled (`warmFrame`: a frame
// called once a step climbs to TurboFan over thousands of frames, `seqv:map`'s layer at 2.08 px past 2,000). The
// positive control goes through the kit's `heap` hook: at each of the stars' frames the same warm-up and measurement
// of the frame allocating an object each (`leaking`) must catch it.

import { expect, it } from 'vitest';
import { RN3, rn3Kit, starsKit } from './conformance-cases.js';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking, warmFrame } from './heap.js';

const N = 2_000;
/** The guard's measurement: warmed and settled by count (`warmFrame`), then `N` frames and `4N` measured. */
const measure = (step: (i: number) => void) => {
  let i = 0;
  const frame = () => step(i++);
  warmFrame(frame);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, N));
};
/** The positive control as the kit's `heap` hook: a frame's check fails unless its leaking frame is caught. */
const control = (frame: () => void) => {
  if (!controlCaught(() => measure(leaking(frame)))) throw new Error('the leaking frame was not caught');
};

it('a stars frame allocates nothing after warm-up, at every zoom and on rings', async () => {
  const heap = (await starsKit(guard)).filter((c) => c.check === 'heap');
  expect(heap.filter((c) => c.geometry === 'circular').length).toBeGreaterThan(0);
  expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
}, 300_000);

it('a stars frame allocating an object each is caught, at every zoom and on rings (positive control)', async () => {
  const heap = (await starsKit(control)).filter((c) => c.check === 'heap');
  expect(heap.filter((c) => c.geometry === 'circular').length).toBeGreaterThan(0);
  expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
}, 300_000);

// RN3's entries through the kit with its heap guard: the rows (single rows at every zoom the kit draws, then the dense
// band) and `seqv:map` (its image built once). Hits and legends are conformance.test.ts's.
for (const rn3 of RN3)
  it(`${rn3[0].name} allocates nothing per frame after warm-up, at every zoom the kit draws`, async () => {
    const heap = (await rn3Kit(rn3, guard)).filter((c) => c.check === 'heap');
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
  }, 300_000);
