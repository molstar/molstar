/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hover frames over split rows on units of their own allocate nothing after warm-up (HR1; a row's own marks): the
// pointer moving over the chains' rows (a unit's residue off its items, an item on them) with a highlight on one unit
// and an accent on the shared sequence, each move an overlay frame on the recording backend: the hit's unit residue,
// the channels' hover rebuilt per own unit, the column split into its rows' bands. The positive control runs the same
// frames allocating an object each (`leaking`) through the same warm-up and measurement, which must catch it. The case
// loads the viewer afresh (`vi.resetModules`).

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
  const { splitUnitsDoc } = await import("./split-units-fixtures.ts");
  const { Dirty } = await import("../src/draw/index.ts");
  const { Inputs, PointerType } = await import("../src/index.ts");
  const t = await lit(splitUnitsDoc({ highlights: true }), { width: 1000, height: 400, dpr: 2 });
  t.ctx.behaviors.get("hoverLabel")?.remove();
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const ys = [0, 1, 2].map((r) => t.yOf(t.section, "t-chains", r));
  const frame = (i: number) => {
    p.x = t.xOf(t.section, 0, 8 + (i % 44));
    p.y = ys[i % 3];
    t.ctx.input.pointer(p);
    t.frame(be, Dirty.overlay);
  };
  const channels = (t.composer as unknown as { channels: { rebuilds: { hover: number } } }).channels;
  return { frame, hovers: () => channels.rebuilds.hover };
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

describe("hover frames over split rows on their own units", () => {
  it("allocate nothing; an object allocated per frame is caught", async () => {
    const { frame, hovers } = await setup();
    const r = run(FRAME_N, frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(hovers()).toBeGreaterThan(FRAME_WARM);
    expect(controlCaught(() => run(FRAME_N, leaking(frame)))).toBe(true);
  }, 240_000);
});
