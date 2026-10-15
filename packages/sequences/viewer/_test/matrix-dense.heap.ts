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
// frame. The positive control: the same steps allocating an object each (`leaking`), through the same warm-up and
// measurement, are caught. Each case loads the viewer afresh (`vi.resetModules`) and is warmed by count until its step
// stops allocating (`warmFrame`); the cases are matrix-dense-cases.ts's, shared with the quick twin
// (matrix-dense.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { DENSE, setup } from './matrix-dense-cases.js';

/** The guards' measurement: warmed by count, `warmed` called, then `FRAME_N` steps and 4× measured. */
function measure(step: (i: number) => void, warmed = () => {}) {
  let i = 0;
  const frame = () => step(i++);
  warmFrame(frame);
  warmed();
  return heapScaling(FRAME_N, frame);
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
  for (const [dpr, stacked] of DENSE)
    it(stacked
      ? "two stacked blocks of one matrix state at DPR 4/3 (height·dpr not whole, the entries' tops apart): no tile built, flat"
      : `the per-device-row tiles, the reduced hit and the outlined cell, at two zooms, DPR ${dpr}`, async () => {
      const { r, images, tiles, hover } = await run(dpr, stacked);
      expect(images).toBeGreaterThan(0);
      expect(tiles).toBe(0);
      expect(hover).toBeGreaterThanOrEqual(0);
      expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    }, 300_000);

  it('the same steps allocating an object each are caught (positive control)', async () => {
    const { step } = await setup(1.5);
    expect(controlCaught(() => measure(leaking(step)))).toBe(true);
  }, 300_000);
});
