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
// shade). Each case loads the viewer and its frame (rn3-fixtures.ts) afresh (`vi.resetModules`), as draws.heap.ts does
// (RN1): a frame built here, shared by the cases, read two copies' double fields through one load site and boxed them
// (45 bytes a frame from the second case on). The positive control: the rows' tiles of held rows allocating an object a
// frame (`leaking`), through the same warm-up and measurement, are caught.

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
const WARM = FRAME_WARM;
/** Composer frames compile concurrently while they run: warmed by count and by time (composer.heap.ts). */
const WARM_MS = FRAME_WARM_MS;
/** The guards' measurement: warmed by count and by time, then `N` frames and `4N` measured. */
const measure = (frame: (i: number) => void) => {
  const t0 = performance.now();
  for (let i = 0; i < WARM || performance.now() - t0 < WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? "" : describeScaling(r, N);
};

describe("RN3's hot paths allocate nothing per frame after warm-up", () => {
  for (const stacked of [false, true])
    for (const labels of [false, true])
      it(`the navigator, ${stacked ? "stacked in lanes" : "side by side"}, ${labels ? "labels" : "hairline"}`, async () => {
        vi.resetModules();
        const { navigatorFrame } = await import("./rn3-fixtures.js");
        expect(guard(await navigatorFrame(stacked, labels))).toBe("");
      }, 120_000);
  for (const remote of [false, true])
    it(`the rows' tiles, ${remote ? "from the region cache" : "held rows"}`, async () => {
      vi.resetModules();
      const { rowsFrame } = await import("./rn3-fixtures.js");
      expect(guard(await rowsFrame(remote))).toBe("");
    }, 120_000);
  it("the rows' tiles, held rows, an object allocated a frame (the positive control)", async () => {
    vi.resetModules();
    const { rowsFrame } = await import("./rn3-fixtures.js");
    const frame = await rowsFrame(false);
    expect(controlCaught(() => measure(leaking(frame)))).toBe(true);
  }, 120_000);
  for (const kind of ["single", "dense", "remote"] as const)
    it(`the rows' own draws, ${kind}: cells, gaps, letters, gutters, the tree and the data columns`, async () => {
      vi.resetModules();
      const { rowsDrawFrame } = await import("./rn3-fixtures.js");
      expect(guard(await rowsDrawFrame(kind))).toBe("");
    }, 120_000);
  for (const bare of [false, true])
    it(`a pan's composer frames, ${bare ? "no track in miniature" : "seqv:map and the navigator's map"}: images kept`, async () => {
      vi.resetModules();
      const { mapFrame } = await import("./rn3-fixtures.js");
      const { frame } = await mapFrame(bare);
      expect(guard(frame)).toBe("");
    }, 120_000);
});
