/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3b's values guard (seqspec:sources#the-viewers-side; the memory rule): ST1's values draw over held regions beside
// the genes, variants and letters, a pan's pass and frame at 100 kb and DPR 2, allocates nothing. ST3a-2's guards leave
// the values track out: with it they measured about 23 bytes a frame (here 31, the region's values and stride read per
// bin, a load
// gone polymorphic beside the items' regions; `Layer.bins` now reads them once a region). Values alone, run first in
// this file, made this case allocate where alone it did not (as ST3a-2's glyph case): ST1's streams-port.heap.ts
// guards values alone. A frame allocating an object a step fails the same guard (the positive control).

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
} from "./heap.js";

const N = 1_000;
/** Six times the shared warm-up: the genes and variants beside the values are drawn in color, which tiers up slowly
 * (streams-items.heap.ts); with the shared warm-up alone this case measured 0 or 42 bytes a step, run to run (ST4b). */
const WARM = 6 * FRAME_WARM;
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame();
  // Settled after the timed warm-up (GS-02): a climb to TurboFan can outlast it.
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};

async function valuesFrame() {
  const { progressGenome } = await import("./stream-progress-fixtures.js");
  const g = await progressGenome(100_000_000, 100_100_000, { values: true, items: true, sequence: true });
  await g.settle();
  let i = 0;
  const frame = () => {
    g.vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, "pointer");
    g.m.pass();
    g.frame();
  };
  return { g, frame };
}

describe("ST1's values draw allocates nothing per frame over held regions", () => {
  it("a pass and a frame, 100 kb at DPR 2, with genes, variants and letters", async () => {
    vi.resetModules();
    const { g, frame } = await valuesFrame();
    frame();
    // Not vacuous: the bars and the genes drawn.
    expect(g.backend.rectItems.n).toBeGreaterThan(500);
    const asked = g.asked();
    expect(guard(frame)).toBe("");
    expect(g.asked()).toBe(asked);
    g.ctx.dispose();
  }, 120_000);

  it("the same frame allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const { frame } = await valuesFrame();
    expect(controlCaught(() => measure(leaking(frame) as () => void))).toBe(true);
  }, 120_000);
});
