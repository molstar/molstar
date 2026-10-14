/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3a's heap guards (seqspec:sources#the-viewers-side; the memory rule): over the streamed genome with a sequence
// track and its regions held, a pan's pass and composer frame allocate nothing: at 10 kb (the letters' cells thinner
// than a device px, read raw per device column from `StreamView.letters`) and at 100 positions (the letters drawn); a
// frame allocating an object a step fails the same guard (the positive control). Each case loads the viewer and its frame
// afresh (`vi.resetModules`).

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
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame();
  // Then until a frame stops allocating: at 10 kb 3 s were about 3,000 frames, and the climb to TurboFan ran about
  // 8,000 more at 55 to 270 bytes a frame (V8 13.6, the pinned Node, 2026-10-13).
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};

describe("the streamed letters' pans allocate nothing per frame once their regions are held", () => {
  for (const at of ["10 kb", "100 b"] as const)
    it(`a pass and a frame, ${at}`, async () => {
      vi.resetModules();
      const { genomeFrame } = await import("./stream-fixtures.js");
      const g = await genomeFrame(at, true, false, true);
      // Not vacuous: the letters drawn (100 b), their cells (10 kb).
      g.frame();
      if (at === "100 b") expect(g.backend.textItems.n).toBeGreaterThan(50);
      else expect(g.backend.rectItems.n).toBeGreaterThan(500);
      const asked = g.requests();
      expect(guard(g.frame)).toBe("");
      expect(g.requests()).toBe(asked);
    }, 120_000);

  it("the same frame allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const { genomeFrame } = await import("./stream-fixtures.js");
    const g = await genomeFrame("100 b", true, false, true);
    expect(controlCaught(() => measure(leaking(g.frame) as () => void))).toBe(true);
  }, 120_000);
});
