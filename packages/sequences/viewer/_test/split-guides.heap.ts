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
// Its own file, the viewer loaded afresh (`vi.resetModules`): cases in one file share V8's feedback.

import { describe, expect, it, vi } from "vitest";
import {
  controlCaught,
  describeScaling,
  FRAME_N,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from "./heap.js";

const N = FRAME_N;

/** The guard's measurement: warmed by count and by time, `warmed` called, then `N` frames and `4N` measured. */
function measure(frame: (i: number) => void, warmed = () => {}) {
  const t0 = performance.now();
  let i = 0;
  for (; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  warmed();
  return heapScaling(N, () => frame(i++));
}

/** Case 13 composed at DPR 2, and a frame: a pan and a whole frame, the guides counted. */
async function setup() {
  vi.resetModules();
  const { Dirty, SkinToken } = await import("../src/draw/index.js");
  const { COPY_BACK, composed, recording } = await import("./composer-fixtures.js");
  const { case13 } = await import("../../core/_test/format/cases.js");
  const t = await composed(case13(), { width: 1000, height: 900, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(3, 200, "api");
  t.frame(be, Dirty.all);
  const GRID = SkinToken.grid;
  let guides = 0;
  let full = 0;
  const frame = (i: number) => {
    t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, "pointer");
    t.frame(be, Dirty.all);
    const r = be.rectItems;
    let n = 0;
    for (let j = 0; j < r.n; j++) if (r.color[j] === GRID && r.alpha[j] < 255 && r.h[j] === 2) n++;
    guides += n;
    if (n >= 21) full++;
  };
  return { frame, full: () => full, guides: () => guides };
}

describe("a split track's row guides allocate nothing after warm-up", () => {
  it("whole frames of case 13's variants grid on a pan, the guides drawn", async () => {
    const { frame, full, guides } = await setup();
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
  }, 120_000);

  it("the same frames allocating an object each are caught (positive control)", async () => {
    const { frame } = await setup();
    expect(controlCaught(() => measure(leaking(frame)))).toBe(true);
  }, 120_000);
});
