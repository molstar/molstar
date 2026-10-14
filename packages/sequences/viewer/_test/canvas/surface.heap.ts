/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The surface's events and frames allocate nothing after warm-up (CS2): pointer moves (hover, the overlay), pans (the
// viewport, the tracks through the copy, the overlay) and both at once, each event ORing flags and each step running
// the one frame they asked for, through the composer and the Canvas2D backend into the canvases.
// Each case loads the viewer afresh (`vi.resetModules`), as the composer's guards do. A browser guard resolves about
// 0.5 MB over the 3 × N extra steps (heap-browser.ts): the surface's own paths are clean to a few bytes a step (hover:
// none measured); a frame of all three passes measured 10 to 70 bytes a step in Chromium where each pass alone measured
// none (CS2 hot paths: the composer's, routed). The positive control: hover frames allocating garbage each, sized to
// the run (`leaking`), warmed and measured alike, are caught.

import { expect, it, vi } from 'vitest';
import { controlCaught, heapScaling, isFlat, leaking } from './heap-browser.js';

/** Steps per measured run: pans draw more, so fewer fit before a scavenge. */
const N = { hover: 2_000, pan: 1_000, both: 1_000 };
const WARM_MS = 3_000;

/** The guard's measurement: warmed by count and by time, `warmed` called, then `n` and `4n` steps; and their count. */
function measure(step: (i: number) => void, n: number, warmed = () => {}) {
  const t0 = performance.now();
  for (let i = 0; i < 1_000 || performance.now() - t0 < WARM_MS; i++) step(i);
  warmed();
  let i = 0;
  const r = heapScaling(n, () => step(i++));
  return { r, steps: i };
}

for (const kind of ['hover', 'pan', 'both'] as const)
  it(`${kind} frames through a mounted surface are flat`, async () => {
    vi.resetModules();
    const { surfaceFrames } = await import('./surface-fixtures.js');
    const t = await surfaceFrames(kind);
    try {
      let frames = 0;
      const { r, steps } = measure(t.step, N[kind], () => {
        frames = t.surface.stats.frames;
      });
      expect(isFlat(r), JSON.stringify(r)).toBe(true);
      // Every pan asked for one frame and drew it; a pointer move only when the hover changed.
      const drawn = t.surface.stats.frames - frames;
      if (kind === 'hover') expect(drawn).toBeGreaterThan(steps / 2);
      else expect(drawn).toBe(steps);
    } finally {
      t.dispose();
    }
  }, 120_000);

it('hover frames allocating garbage sized to the run each are caught (positive control)', async () => {
  vi.resetModules();
  const { surfaceFrames } = await import('./surface-fixtures.js');
  const t = await surfaceFrames('hover');
  try {
    expect(controlCaught(() => measure(leaking(t.step, N.hover), N.hover).r)).toBe(true);
  } finally {
    t.dispose();
  }
}, 120_000);
