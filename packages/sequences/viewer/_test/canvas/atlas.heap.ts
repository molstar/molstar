/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The glyph atlas's path allocates nothing per frame (seqspec:viewer#the-render-backend; CS2 letters), in Chromium: a
// lazy alignment's lettered rows (18 px, cell 14, DPR 2: about 2,000 letters a frame, each a cell hit) drawn whole, and
// panned through the copy (the band copied, its strip's letters drawn) by a pointer and by a program (origin `api`,
// committed at once: its record is built only when read). Each case loads the viewer afresh (`vi.resetModules`). The
// positive controls: the full frames allocating garbage each, sized to the run (`leaking`), warmed and measured alike,
// are caught; so is an object injected every 16th letter in the atlas's draw (about 6 KB a frame).

import { expect, it, vi } from "vitest";
import { controlCaught, heapScaling, isFlat, leaking } from "./heap-browser.js";

const N = { full: 200, pan: 1_000, api: 1_000 };
const WARM_MS = 3_000;

/** The steps the last `measure` ran after its warm-up. */
let measured = 0;
/** The file's warm-up and measurement of a step over `n` steps, as its guards and its controls run them. */
function measure(step: (i: number) => void, n: number, warmed = () => {}) {
  const t0 = performance.now();
  for (let i = 0; i < 200 || performance.now() - t0 < WARM_MS; i++) step(i);
  warmed();
  let i = 0;
  const r = heapScaling(n, () => step(i++));
  measured = i;
  return r;
}

for (const kind of ["full", "pan", "api"] as const)
  it(`${kind} frames of lettered rows through the atlas are flat`, async () => {
    vi.resetModules();
    const { atlasFrames } = await import("./atlas-fixtures.js");
    const t = await atlasFrames(kind);
    try {
      let made = 0;
      let copies = 0;
      const r = measure(t.step, N[kind], () => {
        made = t.shared.atlas?.made ?? 0;
        copies = t.composer.stats.copies;
      });
      expect(isFlat(r), JSON.stringify(r)).toBe(true);
      // Every letter a hit after the warm-up; a pan copies the band.
      expect(made).toBeGreaterThan(20);
      expect(t.shared.atlas?.made).toBe(made);
      if (kind !== "full") expect(t.composer.stats.copies - copies).toBe(measured);
    } finally {
      t.dispose();
    }
  }, 120_000);

it("the same full frames allocating garbage each, sized to the run, are caught (positive control)", async () => {
  vi.resetModules();
  const { atlasFrames } = await import("./atlas-fixtures.js");
  const t = await atlasFrames("full");
  try {
    expect(controlCaught(() => measure(leaking(t.step, N.full), N.full))).toBe(true);
  } finally {
    t.dispose();
  }
}, 120_000);

it("the guard sees an allocation injected in the atlas's draw (positive control)", async () => {
  vi.resetModules();
  const { atlasFrames } = await import("./atlas-fixtures.js");
  const t = await atlasFrames("full");
  try {
    const atlas = t.shared.atlas as NonNullable<typeof t.shared.atlas>;
    const draw = atlas.draw.bind(atlas);
    let n = 0;
    let kept: unknown = null;
    atlas.draw = (g, b, i) => {
      if ((++n & 15) === 0) kept = { n, i, b };
      return draw(g, b, i);
    };
    expect(controlCaught(() => measure(t.step, N.full))).toBe(true);
    expect(kept).not.toBeNull();
  } finally {
    t.dispose();
  }
}, 120_000);
