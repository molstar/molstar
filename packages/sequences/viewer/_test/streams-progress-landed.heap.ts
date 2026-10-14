/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3b's fade guard (seqspec:sources#levels-and-the-zoom-gate; the memory rule), its own file (streams-progress.heap.ts
// says why): every region of the streamed genome landed after a second of pending frames, the clock held so the fade
// never ends; an overlay frame allocates nothing. A frame allocating an object a step fails the same guard (the positive
// control).

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
/** The file's warm-up, settle and measurement, as its guards and its control run them. */
const measure = (frame: () => void) => {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame();
  // Settled after the timed warm-up (GS-02): a climb to TurboFan can outlast it.
  settleFrame(frame, FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};

describe("the overlay's fade-in allocates nothing per frame", () => {
  it("every region landed, fading", async () => {
    vi.resetModules();
    const { overlayFrame, rectsOf, SkinToken } = await import("./stream-progress-fixtures.ts");
    const { g, frame, next } = await overlayFrame("landed");
    frame();
    // Not vacuous: the shade over every region just drawn; a frame asked while they fade.
    expect(rectsOf(g.backend, SkinToken.pending).n).toBeGreaterThan(2);
    expect(next[0]).toBeGreaterThan(0);
    expect(guard(frame)).toBe("");
    g.ctx.dispose();
  }, 120_000);

  it("the same frame allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const { overlayFrame } = await import("./stream-progress-fixtures.ts");
    const { frame } = await overlayFrame("landed");
    expect(controlCaught(() => measure(leaking(frame) as () => void))).toBe(true);
  }, 120_000);
});
