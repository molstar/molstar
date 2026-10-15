/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN3's hot paths allocate nothing per frame after warm-up (seqspec:viewer#extension-points;
// seqspec:viewer#the-render-backend, seqspec:sources#levels-and-the-zoom-gate): the axis view's map synced through a
// pan (the window mapped through the blocks), the navigator's axis and its window with the channels' marks per lane,
// the labels axis, the grips; the rows' tiles painted with their images made (coarse under fine, the overview's pending
// shade). Each case loads the viewer and its frame afresh (rn3-cases.ts, shared with the quick twin rn3.quick.ts). The
// positive control: the rows' tiles of held rows allocating an object a frame (`leaking`), through the same warm-up and
// measurement, are caught. The warm-up is by count (`frameScaling`: the composer's frames compile while they run).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, frameScaling, isFlat, leaking } from './heap.js';
import { RN3_CASES, rowsHeld } from './rn3-cases.js';

const N = FRAME_N;

describe("RN3's hot paths allocate nothing per frame after warm-up", () => {
  for (const c of RN3_CASES)
    it(c.name, async () => {
      const r = frameScaling(await c.open(), N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
    }, 300_000);

  it("the rows' tiles, held rows, an object allocated a frame (the positive control)", async () => {
    const frame = await rowsHeld();
    expect(controlCaught(() => frameScaling(leaking(frame), N))).toBe(true);
  }, 300_000);
});
