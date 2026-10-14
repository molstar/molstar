/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hover frames over sized shapes allocate nothing after warm-up (HR2): the pointer moving over dots, bands, boundaries'
// glyphs (the end boundary's too, as ranges, points and stacked variants, the pointer around its edge), pair ends and wide
// arcs, each move a hit (the layer's themes set for its `hit`, the glyphs near the point
// measured), the hover rebuilt and an overlay frame on the recording backend; and over an own unit's row in
// `UI.OwnRowHover` "residue" mode. The positive control runs the first case's frames allocating an object each
// (`leaking`) through the same warm-up and measurement, which must catch it. Each case loads the viewer afresh
// (`vi.resetModules`).

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
  const { hitGeometryDoc, S } = await import("./hit-geometry-fixtures.ts");
  const { Dirty } = await import("../src/draw/index.ts");
  const { Inputs, PointerType } = await import("../src/index.ts");
  const t = await lit(hitGeometryDoc(), { width: 900, height: 600, dpr: 2 });
  t.ctx.behaviors.get("hoverLabel")?.remove();
  // The whole sequence in view (fit), so both ends' boundaries are on screen.
  t.vp?.setView(1, 0, "api");
  const be = recording(COPY_BACK);
  t.frame(be, Dirty.all);
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(S);
  p.type = PointerType.move;
  const ys = ["t-dots", "t-bars", "t-ticks", "t-pairs", "t-arcs"].map((r) => t.yOf(S, r));
  const ends = ["t-ticks", "t-tick-points", "t-tick-variants"].map((r) => t.yOf(S, r));
  const x0 = t.xOf(S, 0, 8);
  // The edge after the last position: the end boundary's glyph, the pointer on either side of it.
  const vp = t.ctx.layout.axis(S);
  const xEnd = vp ? vp.x(0, (vp.shownOf(0, 399) + 1) as never) : 0;
  const frame = (i: number) => {
    if (i % 4 === 3) {
      p.x = xEnd + ((i % 9) - 5) * 1.5;
      p.y = ends[i % 3] + ((i % 5) - 2) * 2;
    } else {
      p.x = x0 + (i % 60) * 2.3;
      p.y = ys[i % 5] + ((i % 7) - 3) * 2;
    }
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

describe("hover frames over sized shapes", () => {
  it("allocate nothing; an object allocated per frame is caught", async () => {
    const { frame, hovers } = await setup();
    const r = run(FRAME_N, frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(hovers()).toBeGreaterThan(FRAME_WARM);
    expect(controlCaught(() => run(FRAME_N, leaking(frame)))).toBe(true);
  }, 240_000);

  it('"residue" mode over an own unit\'s row allocates nothing', async () => {
    vi.resetModules();
    const { lit } = await import("./overlay-fixtures.ts");
    const { recording, COPY_BACK } = await import("./composer-fixtures.ts");
    const { splitUnitsDoc } = await import("./split-units-fixtures.ts");
    const { Dirty } = await import("../src/draw/index.ts");
    const { Inputs, PointerType, SeqViewerConfig } = await import("../src/index.ts");
    const t = await lit(splitUnitsDoc({ highlights: true }), { width: 1000, height: 400, dpr: 2 });
    t.ctx.behaviors.get("hoverLabel")?.remove();
    t.ctx.config.set(SeqViewerConfig.UI.OwnRowHover, "residue");
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
    const r = run(FRAME_N, frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
  }, 240_000);
});
