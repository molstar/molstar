/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames of arcs whose partners lie beyond the window allocate nothing after warm-up: case 7 zoomed into barstar, every
// contact's partner off-window in barnase, so each arc is pushed as its pieces in view (cut at the window and the
// band's far edge) with a wedge where it meets the edge; pans move the view each frame (`Arcs.frame` recomputed), and
// a hover on an end lights its arcs on the overlay. The positive control runs the same frames with an object allocated
// each (`leaking`) through the same warm-up and measurement, which must catch it. The case loads the viewer afresh
// (`vi.resetModules`), as the composer's guards do.

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
} from "./heap.ts";

async function setup() {
  vi.resetModules();
  const { lit } = await import("./overlay-fixtures.ts");
  const { recording, COPY_BACK } = await import("./composer-fixtures.ts");
  const { layerOf, record } = await import("./draw-fixtures.ts");
  const { case7 } = await import("../../core/_test/format/cases.ts");
  const t = await lit(case7(), { width: 1000, height: 400, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.zoomTo(1, 34, 86);
  await t.ctx.graph.update();
  t.frame(be, 7);
  // The hover label's strings are one allocation per hovered element by design (IN2's guard): the overlay alone.
  t.ctx.behaviors.get("hoverLabel")?.remove();
  const { k, section } = layerOf(t.ctx, "Contacts");
  expect(section).toBe(t.section);
  const s = t.ctx.graph.assembled.value?.view?.prepared[k].state as { frame: Float64Array };
  const glyphs = Array.from(record(t.ctx, "Contacts").rec.glyphItems.x.subarray(0, 12), (x) => x / 2);
  return { t, be, s, glyphs };
}

/** The guard's measurement: warmed by count and by time, then `heapScaling` over `n` frames and `4n`. */
function run(n: number, frame: (i: number) => void) {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  let i = 0;
  return heapScaling(n, () => frame(i++));
}

describe("frames of arcs whose partners are beyond the window", () => {
  it("pans, a hover on an end and its overlay allocate nothing; an object allocated per frame is caught", async () => {
    const { t, be, s, glyphs } = await setup();
    expect(glyphs.length).toBeGreaterThan(5);
    const into = (await import("./composer-fixtures.ts")).recording();
    const y = s.frame[0];
    const frame = (i: number) => {
      t.vp?.panBy((i & 1) === 0 ? 5 : -5, 0, "pointer");
      t.frame(be, 7);
      t.hover(glyphs[i % glyphs.length], y);
      t.overlay(into);
    };
    const r = run(FRAME_N, frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(controlCaught(() => run(FRAME_N, leaking(frame)))).toBe(true);
  }, 240_000);
});
