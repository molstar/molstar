/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fingerprint's frames and hits allocate nothing after warm-up (seqspec:viewer#extension-points): the conformance
// kit's guard (`heapScaling` per frame, at the kit's three zooms, `tiny` too) on case 6, on 20 poses × 8 kinds
// (stripes) and at a short track (compact), on mixed data and the absent property (the tiles beside the rows), and on
// 20,000 positions zoomed out (merged cells), case 6 and the wide one reversed (AX1); and a million hits over every row
// and stripe, compact rows and a block without the property, each row and each kind hit somewhere first (the cases in
// fingerprint-cases.ts, shared with the quick twin). Each case loads the viewer afresh (`vi.resetModules`): a frame
// built in this module would box from the second case on (047, RN3). Each frame is warmed and settled by count
// (`warmFrame`). The positive control, last: the kit's first frame on the 20 poses allocating an object each
// (`leaking`), warmed and measured alike, is caught.

import { describe, expect, it } from 'vitest';
import { controlCaught, leaking, warmFrame } from '../heap.js';
import { fingerprintKit, HIT_SCENES, hitSweep, SCENES } from './fingerprint-cases.js';
import { describeScaling, heapScaling, isFlat } from './heap.js';

const N = 1_000;
/** Warmed and settled by count (`warmFrame`: the wide scene's tiny draw climbed past 2,000 frames). */
const measure = (frame: () => void) => {
  warmFrame(frame);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, N));
};

describe("the fingerprint's frames allocate nothing after warm-up, at every zoom the kit draws", () => {
  for (const name of SCENES)
    it(name, async () => {
      const heap = (await fingerprintKit(name, guard)).filter((c) => c.check === 'heap');
      expect(heap.length).toBeGreaterThanOrEqual(6);
      expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
    }, 300_000);
});

describe('the hit allocates nothing: a million over every row and stripe', () => {
  it.each(HIT_SCENES)(
    '%s',
    async (name) => {
      const { step } = await hitSweep(name);
      // Warmed past the hit's tier-up (compact: about 1.4 bytes a hit for 1.5 million hits, then none), a thousand a
      // frame.
      let w = 0;
      warmFrame(() => {
        for (let q = 0; q < 1_000; q++) step(w++);
      });
      const r = heapScaling(250_000, step);
      expect(isFlat(r), describeScaling(r, 250_000)).toBe(true);
    },
    300_000,
  );
});

describe('the positive control', () => {
  it("the kit's first frame on the 20 poses, allocating an object each, is caught", async () => {
    const caught: boolean[] = [];
    await fingerprintKit('many', (frame) => {
      if (caught.length === 0) caught.push(controlCaught(() => measure(leaking(frame) as () => void)));
    });
    expect(caught).toEqual([true]);
  }, 300_000);
});
