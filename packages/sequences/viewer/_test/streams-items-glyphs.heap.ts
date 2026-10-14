/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3a-2's heap guards (seqspec:sources#the-viewers-side; the memory rule): at 5 kb, regions held, a pan's pass and
// composer frame allocate nothing (variants as glyphs region by region, stacked or as points, genes as items, letters as
// cells); with every region pending, a pan's frame draws the placeholder bands and allocates nothing. A frame
// allocating an object a step fails the same guard (the positive control). Apart from streams-items.heap.ts: run after
// its cases in one isolate, this one allocated.

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

describe("the streamed glyphs and placeholders allocate nothing per frame", () => {
  it("a pass and a frame, 5 kb", async () => {
    vi.resetModules();
    const { itemsFrame } = await import("./stream-items-fixtures.js");
    const g = await itemsFrame(100_040_000, 100_045_000);
    g.frame();
    expect(g.backend.glyphItems.n).toBeGreaterThan(0);
    const asked = g.asked();
    expect(guard(g.frame)).toBe("");
    expect(g.asked()).toBe(asked);
    g.dispose();
  }, 120_000);

  it("a pass and a frame, 5 kb, the variants as points (unstacked glyphs, region by region)", async () => {
    vi.resetModules();
    const { itemsFrame } = await import("./stream-items-fixtures.js");
    const g = await itemsFrame(100_040_000, 100_045_000, false, "points");
    g.frame();
    expect(g.backend.glyphItems.n).toBeGreaterThan(0);
    const asked = g.asked();
    expect(guard(g.frame)).toBe("");
    expect(g.asked()).toBe(asked);
    g.dispose();
  }, 120_000);

  it("a frame with every region pending (placeholder bands)", async () => {
    vi.resetModules();
    const { itemsFrame } = await import("./stream-items-fixtures.js");
    const { SkinToken } = await import("../src/draw/backend.js");
    const g = await itemsFrame(100_000_000, 100_100_000, true);
    g.frame();
    let bands = 0;
    for (let k = 0; k < g.backend.rectItems.n; k++) if (g.backend.rectItems.color[k] === SkinToken.pending) bands++;
    expect(bands).toBeGreaterThan(2);
    expect(guard(g.frame)).toBe("");
    g.dispose();
  }, 120_000);

  it("the same frame allocating an object a step is caught (positive control)", async () => {
    vi.resetModules();
    const { itemsFrame } = await import("./stream-items-fixtures.js");
    const g = await itemsFrame(100_040_000, 100_045_000);
    expect(controlCaught(() => measure(leaking(g.frame) as () => void))).toBe(true);
  }, 120_000);
});
