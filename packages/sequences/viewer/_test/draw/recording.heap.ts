/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The batches and the recording backend allocate nothing once warm (seqspec:viewer#the-render-backend): recorded
// frames, each of 4,000 rects and 1,000 text runs with every other batch kind, an image and a pan's copy (the step is
// recording-cases.ts's, shared with the quick twin, recording.quick.ts). A step is a frame, over runs of 1,000 and
// 4,000 (RB0b: BK1's 210 bytes a frame grew 210 KB in one run of 1,000, under LIMIT). The positive control: the same
// frames allocating an object each (`leaking`), through the same measurement, are caught.

import { expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../heap.js';
import { recorded } from './recording-cases.js';

it('recorded frames allocate nothing after warm-up', () => {
  const { step, rects } = recorded();
  const s = heapScaling(1000, step);
  expect(rects()).toBe(4001);
  expect(isFlat(s), describeScaling(s, 1000)).toBe(true);
  expect(controlCaught(() => heapScaling(1000, leaking(step)))).toBe(true);
});
