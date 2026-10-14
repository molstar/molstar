/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The overlay's hover frames allocate nothing after warm-up (CS1b): the pointer moving over 2 × 200,000 positions
// zoomed out with a selection and a highlight on them, and over an alignment's rows with a row's residues selected and
// cells highlighted (its row marks), each move an overlay frame on the recording backend: the channels' hover rebuilt,
// row marks, position marks per lane, the hovered row's band, renderers' overlays and the ring. Each frame is warmed by
// count and by time, then settled (`settleFrame`), as composer.heap.ts's are. The positive control: the hover frames
// over the rows allocating an object each (`leaking`), through the same warm-up, settle and measurement, are caught.
//
// Each case loads the viewer afresh (`vi.resetModules`), as composer.heap.ts does.

import { describe, expect, it, vi } from "vitest";
import {
  controlCaught,
  describeScaling,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from "./heap.ts";

const N = 1_000;
const WARM = FRAME_WARM;
const WARM_MS = FRAME_WARM_MS;
/** The guards' measurement: warmed by count and by time, settled (`settleFrame`), then `N` frames and `4N` measured. */
function measure(frame: (i: number) => void) {
  const t0 = performance.now();
  let i = 0;
  for (; i < WARM || performance.now() - t0 < WARM_MS; i++) frame(i);
  settleFrame(() => frame(i++), FRAME_SETTLE_MS);
  return heapScaling(N, () => frame(i++));
}

describe("overlay frames allocate nothing after warm-up", () => {
  it.each(["scale", "rows"] as const)(
    "hover frames over %s with a selection and a highlight",
    async (kind) => {
      vi.resetModules();
      const { hoverFrames } = await import("./overlay-fixtures.ts");
      const { frame, hovers } = await hoverFrames(kind);
      const r = measure(frame);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      // Not vacuous: the hovered element changed on most moves.
      expect(hovers()).toBeGreaterThan(WARM);
    },
    120_000,
  );

  it("hover frames over rows allocating an object each are caught (positive control)", async () => {
    vi.resetModules();
    const { hoverFrames } = await import("./overlay-fixtures.ts");
    const { frame } = await hoverFrames("rows");
    expect(controlCaught(() => measure(leaking(frame)))).toBe(true);
  }, 120_000);
});
