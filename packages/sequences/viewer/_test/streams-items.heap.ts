/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3a-2's heap guards (seqspec:sources#the-viewers-side; the memory rule), regions held: a pan's pass and composer
// frame allocate nothing at 100 kb (genes as items region by region, variants as density), at 1 Mb (both as density
// from `count` bins) and at 200 Mb (the letters gate's hatch and hint). A frame allocating an object a step fails the same
// guard (the positive control). The glyphs and the placeholders are streams-items-glyphs.heap.ts's.

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
/**
 * ST3b: streamed items drawn in color (a field theme's default reaches them now) tier up slowly: a frame measured 160,
 * then 16, then 4 bytes over runs of 5,000 frames (a GC between each), and so at e3716faa with a uniform color; with the
 * shared warm-up alone the 5 kb and 200 Mb cases measured 35 to 45 bytes a step. These frames warm six times as long.
 */
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

describe("the streamed items' pans allocate nothing per frame once their regions are held", () => {
  for (const [label, start, stop] of [
    ["100 kb", 100_000_000, 100_100_000],
    ["1 Mb", 100_000_000, 101_000_000],
    ["200 Mb", 20_000_000, 220_000_000],
  ] as const)
    it(`a pass and a frame, ${label}`, async () => {
      vi.resetModules();
      const { itemsFrame } = await import("./stream-items-fixtures.ts");
      const g = await itemsFrame(start, stop);
      g.frame();
      // Not vacuous: genes and variants drawn (rects), at 200 Mb the hatch too.
      expect(g.backend.rectItems.n).toBeGreaterThan(100);
      if (start === 20_000_000) expect(g.backend.pathItems.n).toBeGreaterThan(0);
      const asked = g.asked();
      expect(guard(g.frame)).toBe("");
      expect(g.asked()).toBe(asked);
      g.dispose();
    }, 120_000);

  it("the same frame allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const { itemsFrame } = await import("./stream-items-fixtures.ts");
    const g = await itemsFrame(100_000_000, 100_100_000);
    expect(controlCaught(() => measure(leaking(g.frame) as () => void))).toBe(true);
  }, 120_000);
});
