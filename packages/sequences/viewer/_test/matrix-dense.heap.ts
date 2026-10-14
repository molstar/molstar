/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// KS2's dense matrix allocates nothing and builds no tile per frame or per pointer move after warm-up
// (seqspec:viewer#performance): whole frames of a 3,000-state genotype matrix drawn through its per-device-row tiles,
// at DPR 1.5 and 1, the pointer moving over it (the reduced hit, the hover and the outlined cell), at two zooms (cells
// and two positions a device column); and two stacked blocks of one state at DPR 4/3 (height·dpr not whole), where a
// tile height from each entry's rounded top and bottom alternated between the entries and rebuilt every tile each
// frame. The positive
// control: the same steps allocating an object each (`leaking`), through the same warm-up and measurement, are caught.
// Each case loads the viewer afresh (`vi.resetModules`).

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

/** The guards' measurement: warmed by count and by time, `warmed` called, then `FRAME_N` steps and 4× measured. */
function measure(step: (i: number) => void, warmed = () => {}) {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) step(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => step(settled++), FRAME_SETTLE_MS);
  warmed();
  let i = 0;
  return heapScaling(FRAME_N, () => step(i++));
}

/** The genotype matrix composed, and a step: a pointer move and a whole frame, the zoom switched every 64 steps. */
async function setup(dpr: number, stacked = false) {
  vi.resetModules();
  const { composed, recording } = await import("./composer-fixtures.ts");
  const { genotypeDoc } = await import("./genotype-fixtures.ts");
  const { entryOf } = await import("./draw-fixtures.ts");
  const { Inputs, PointerType } = await import("../src/index.ts");
  const { MatrixVisits } = await import("../src/draw/matrix-tiles.ts");
  const doc = genotypeDoc({ samples: 3000, positions: 2000, stacked });
  const t = await composed(doc, { width: 1000, height: 400, dpr });
  const vp = t.vp;
  if (!vp) throw new Error("no viewport");
  const be = recording();
  t.ctx.behaviors.get("hoverLabel")?.remove();
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex(t.section);
  p.type = PointerType.move;
  const cells = [vp.live.current.cell * 40, vp.live.current.cell];
  const e = entryOf(t.ctx, "G");
  const step = (i: number) => {
    // Every 64 steps the other zoom; the pointer down and across the matrix's rows between.
    if ((i & 63) === 0) vp.setView(cells[(i >> 6) & 1], 0, "api");
    p.x = e.left + 3.5 + ((i * 7) % 600);
    p.y = e.top + 0.25 + ((i * 0.37) % (e.height - 0.5));
    t.ctx.input.pointer(p);
    t.composer.redraw();
    t.frame(be, 7);
  };
  return { step, visits: MatrixVisits, hover: () => t.ctx.hover.live.current.state };
}

async function run(dpr: number, stacked = false) {
  const { step, visits, hover } = await setup(dpr, stacked);
  let images = 0;
  let built = 0;
  const r = measure(step, () => {
    images = visits.images;
    built = visits.tiles;
  });
  return { r, images: visits.images - images, tiles: visits.tiles - built, hover: hover() };
}

describe("KS2's dense matrix frames allocate nothing after warm-up", () => {
  it.each([1.5, 1])(
    "the per-device-row tiles, the reduced hit and the outlined cell, at two zooms, DPR %f",
    async (dpr) => {
      const { r, images, tiles, hover } = await run(dpr);
      expect(images).toBeGreaterThan(0);
      expect(tiles).toBe(0);
      expect(hover).toBeGreaterThanOrEqual(0);
      expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    },
    120_000,
  );

  it("two stacked blocks of one matrix state at DPR 4/3 (height·dpr not whole, the entries' tops apart): no tile built, flat", async () => {
    const { r, images, tiles } = await run(4 / 3, true);
    expect(images).toBeGreaterThan(0);
    expect(tiles).toBe(0);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
  }, 120_000);

  it("the same steps allocating an object each are caught (positive control)", async () => {
    const { step } = await setup(1.5);
    expect(controlCaught(() => measure(leaking(step)))).toBe(true);
  }, 120_000);
});
