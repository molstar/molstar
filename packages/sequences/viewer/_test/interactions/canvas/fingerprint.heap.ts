/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The fingerprint's frames through the Canvas2D backend allocate nothing per frame (seqspec:viewer#the-render-backend),
// in Chromium (`HEAP=1`: heap-browser.ts): 20 poses × 8 kinds (stripes), the same in a short track (compact), mixed
// data and the absent property (tiles beside rows), and 20,000 positions zoomed out (merged cells). The positive
// control: the frames allocating garbage each, sized to the run (`leaking`), warmed and measured alike, are caught.
// Each case loads the viewer afresh (`vi.resetModules`).

import { expect, it, vi } from "vitest";
import { controlCaught, heapScaling, isFlat, leaking } from "../../canvas/heap-browser.js";

const N = 1_000;
const WARM = 2_000;

/** The file's warm-up and measurement of a frame, as its guards and its control run them. */
function measure(step: (i: number) => void) {
  for (let i = 0; i < WARM; i++) step(i);
  return heapScaling(N, step);
}

/**
 * What each case's frame draws, so the guard measures the path it names: stripes on grounds (the group walk), compact
 * rects (no ground), tiles beside rows, tiles alone, and merged thin cells.
 */
const DRAWS = {
  many: (d: Drawn) => d.grounds > 0 && d.rects > d.grounds && d.images === 0,
  compact: (d: Drawn) => d.grounds === 0 && d.rects > 0 && d.images === 0,
  mixed: (d: Drawn) => d.grounds > 0 && d.rects > d.grounds && d.images > 0,
  absent: (d: Drawn) => d.rects === 0 && d.images > 0,
  wide: (d: Drawn) => d.grounds > 0 && d.rects > d.grounds,
};
interface Drawn {
  rects: number;
  grounds: number;
  images: number;
}

for (const name of ["many", "compact", "mixed", "absent", "wide"] as const)
  it(`${name}: frames through Canvas2D are flat`, async () => {
    vi.resetModules();
    const { fingerprintFrames } = await import("./fingerprint-fixtures.js");
    const t = await fingerprintFrames(name);
    try {
      const d = t.drawn();
      expect(DRAWS[name](d), JSON.stringify(d)).toBe(true);
      const r = measure(t.step);
      expect(isFlat(r), JSON.stringify(r)).toBe(true);
    } finally {
      t.dispose();
    }
  }, 120_000);

it("the same frames allocating garbage each, sized to the run, are caught (positive control)", async () => {
  vi.resetModules();
  const { fingerprintFrames } = await import("./fingerprint-fixtures.js");
  const t = await fingerprintFrames("many");
  try {
    expect(controlCaught(() => measure(leaking(t.step, N)))).toBe(true);
  } finally {
    t.dispose();
  }
}, 120_000);
