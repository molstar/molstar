/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The batches and the recording backend allocate nothing once warm (seqspec:viewer#the-render-backend): recorded
// frames, each of 4,000 rects and 1,000 text runs with every other batch kind, an image and a pan's copy. A step is a
// frame, over runs of 1,000 and 4,000 (RB0b: BK1's 210 bytes a frame grew 210 KB in one run of 1,000, under LIMIT).
// The positive control: the same frames allocating an object each (`leaking`), through the same measurement, are caught.

import { expect, it } from 'vitest';
import { RecordingBackend } from '../../src/testing/index.js';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../heap.js';
import { drawFrame, Frame, PIXELS } from '../support/frame.js';

it('recorded frames allocate nothing after warm-up', () => {
  const r = new RecordingBackend();
  const f = new Frame(4000, 1000);
  const img = r.images.fromPixels(2, 2, PIXELS);
  const step = (i: number) => {
    r.reset();
    drawFrame(r, f, i, img);
  };
  const s = heapScaling(1000, step);
  expect(r.rectItems.n).toBe(4001);
  expect(isFlat(s), describeScaling(s, 1000)).toBe(true);
  expect(controlCaught(() => heapScaling(1000, leaking(step)))).toBe(true);
});
