/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1's heap guards (seqspec:sources#the-viewers-side; the memory rule): over the streamed hg38-sized genome with its
// regions held, a pan's pass, and a pan's pass and composer frame (the values drawn from the regions,
// `StreamView.values` per device column), allocate nothing, at the whole view and at 10 kb; a frame that allocates
// fails the same guard (the positive control). Each case loads the viewer and its frame afresh (`vi.resetModules`), as
// rn3.heap.ts does.

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
const WARM = FRAME_WARM;
const WARM_MS = FRAME_WARM_MS;
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < WARM || performance.now() - t0 < WARM_MS; i++) frame();
  // Then until a frame stops allocating: at 10 kb with the frame drawn, 3 s were about 9,000 frames, and the climb to
  // TurboFan ran about 12,000 more at 8 to 122 bytes a frame (V8 13.6, the pinned Node, 2026-10-13).
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};

describe("the streamed genome's pans allocate nothing per frame once its regions are held", () => {
  for (const at of ["whole", "10 kb"] as const)
    for (const draw of [false, true])
      it(`${draw ? "a pass and a frame" : "a pass"}, ${at}`, async () => {
        vi.resetModules();
        const { genomeFrame } = await import("./stream-fixtures.js");
        const g = await genomeFrame(at, draw);
        const asked = g.requests();
        expect(guard(g.frame)).toBe("");
        expect(g.requests()).toBe(asked);
        expect(g.manager.stats.current.hits).toBeGreaterThan(0);
      }, 120_000);

  // ST1b-a: with dynamic nodes held (the manifest's genes drawn, every data region's peak), the same.
  for (const at of ["whole", "10 kb"] as const)
    it(`a pass and a frame, ${at}, dynamic nodes held`, async () => {
      vi.resetModules();
      const { genomeFrame } = await import("./stream-fixtures.js");
      const g = await genomeFrame(at, true, true);
      expect(g.manager.nodes.entries.size).toBeGreaterThan(at === "whole" ? -1 : 0);
      const asked = g.requests();
      expect(guard(g.frame)).toBe("");
      expect(g.requests()).toBe(asked);
    }, 120_000);

  it("the same frame allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const { genomeFrame } = await import("./stream-fixtures.js");
    const g = await genomeFrame("10 kb", true);
    expect(controlCaught(() => measure(leaking(g.frame) as () => void))).toBe(true);
  }, 120_000);
});
