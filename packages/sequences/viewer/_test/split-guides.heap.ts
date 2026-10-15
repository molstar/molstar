/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A split track's row guides (seqspec:viewer#render-props) allocate nothing per frame after warm-up: case 13's variants
// grid (22 rows) redrawn whole on a pan at DPR 2, each frame drawing the guides under the rows' data, on the
// recording backend. Not vacuous, counted over the measured frames alone: most frames drew the 21 guides (the grid tone
// at the guide's alpha, 2 device px tall; a pan's strips and the whole frame each draw them). The positive control: the
// same frames allocating an object each (`leaking`), through the same warm-up and measurement, are caught.
//
// The viewer is loaded afresh (`vi.resetModules`) and the frames are warmed by count (`warmFrame`: until they stop
// allocating); the case is split-guides-cases.ts's, shared with the quick twin (split-guides.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { guidesFrame } from './split-guides-cases.js';

const N = FRAME_N;

/** The guard's measurement: warmed by count, `warmed` called, then `N` frames and `4N` measured. */
function measure(frame: (i: number) => void, warmed = () => {}) {
  let i = 0;
  const step = () => frame(i++);
  warmFrame(step);
  warmed();
  return heapScaling(N, step);
}

describe("a split track's row guides allocate nothing after warm-up", () => {
  it("whole frames of case 13's variants grid on a pan, the guides drawn", async () => {
    const { frame, full, guides } = await guidesFrame();
    let full0 = 0;
    let measured = 0;
    const r = measure(
      (i) => {
        frame(i);
        measured++;
      },
      () => {
        full0 = full();
        measured = 0;
      },
    );
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    // Not vacuous, over the measured frames: the frames drew every guide.
    expect(measured).toBeGreaterThanOrEqual(5 * N);
    expect(full() - full0).toBeGreaterThan(0.9 * measured);
    expect(guides()).toBeGreaterThan(0);
  }, 300_000);

  it('the same frames allocating an object each are caught (positive control)', async () => {
    const { frame } = await guidesFrame();
    expect(controlCaught(() => measure(leaking(frame)))).toBe(true);
  }, 300_000);
});
