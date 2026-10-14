/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VL-13: a stacked density column in a band thinner than its colors stays in the band: every rect a positive height
// inside it (an SVG export gets no negative or out-of-band rect).

import { expect, it } from "vitest";
import {
  DensityKind,
  type DrawView,
  densityAdd,
  densityDraw,
  densitySlot,
  densityStart,
  RectBatch,
} from "../src/draw/api.js";

it("six colors in a band of three device rows", () => {
  densityStart({ left: 0, right: 10, dpr: 1 } as DrawView);
  for (let s = 0; s < 6; s++) densityAdd(3, densitySlot(0x100000 * (s + 1)), 1);
  const rects = new RectBatch();
  densityDraw(rects, 10, 13, DensityKind.stacked, -1);
  const out = Array.from({ length: rects.n }, (_, i) => [rects.y[i], rects.h[i]]);
  expect(out.length).toBeGreaterThan(0);
  for (const [y, h] of out) {
    expect(h).toBeGreaterThan(0);
    expect(y).toBeGreaterThanOrEqual(10);
    expect(y + h).toBeLessThanOrEqual(13);
  }
  // The band filled once, top to bottom.
  expect(out.reduce((a, [, h]) => a + h, 0)).toBe(3);
});
