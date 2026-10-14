/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards on the surface's frame path over rows (CS2; seqspec:viewer#the-render-backend): a pan and scroll frame
// over 12,000 rows (one entry per row, dense rows in one call) with the tree, data columns and lanes, and over lazy
// rows with their tiles cached (the rows-tiles paint path), and a hover frame over the tree, allocate nothing after
// warm-up (browser, `heapScaling`; each case loads the viewer afresh, as the composer's guards do). The positive
// control: the lazy rows' frames allocating garbage each, sized to the run (`leaking`), warmed and measured alike, are
// caught; it runs on the lazy rows' 2,000 frames (150 over 12,000 rows, which a longer run collects in).

import { expect, it, vi } from 'vitest';
import { controlCaught, heapScaling, isFlat, leaking } from './heap-browser.js';

const WARM_MS = 3_000;
/** Steps per measured run: longer runs over 12,000 rows were collected mid-run (a collection hides allocation, heap-browser.ts). */
const N = 150;
/** Lazy rows' steps: a frame of cached tiles is cheaper, so the run is longer. */
const LAZY_N = 2_000;

type Frames = { step: (i: number) => void; dispose: () => void };

/** The guard's measurement: warmed by count and by time, then `n` and `4n` frames of `step`. */
function measure(step: (i: number) => void, n = N) {
  const t0 = performance.now();
  for (let i = 0; i < 1_000 || performance.now() - t0 < WARM_MS; i++) step(i);
  let i = 0;
  return heapScaling(n, () => step(i++));
}
/** A fresh viewer's frames (`vi.resetModules`), measured flat. */
async function guard(make: () => Promise<Frames>, n = N) {
  vi.resetModules();
  const t = await make();
  try {
    const r = measure(t.step, n);
    expect(isFlat(r), JSON.stringify(r)).toBe(true);
  } finally {
    t.dispose();
  }
}
const lazy = async () => {
  const { lazyFrames } = await import('./rows-fixtures.js');
  return lazyFrames();
};

for (const scale of [1, 0.12])
  it(`a pan and scroll frame over 12,000 rows at row scale ${scale} is flat`, async () => {
    await guard(async () => {
      const { familyFrames } = await import('./rows-fixtures.js');
      return familyFrames(scale);
    });
  }, 120_000);

it('a pan and scroll frame over lazy rows with their tiles cached is flat', async () => {
  await guard(lazy, LAZY_N);
}, 120_000);

it("the same lazy rows' frames allocating garbage sized to the run each are caught (positive control)", async () => {
  vi.resetModules();
  const t = await lazy();
  try {
    expect(controlCaught(() => measure(leaking(t.step, LAZY_N), LAZY_N))).toBe(true);
  } finally {
    t.dispose();
  }
}, 120_000);

it('a hover frame over the tree (clades lit, their rows banded) is flat', async () => {
  vi.resetModules();
  const { treeHover } = await import('./rows-fixtures.js');
  const t = await treeHover();
  try {
    for (let i = 0; i < 800; i++) t.step(i);
    expect(t.nodes()).toBeGreaterThan(100);
    const t0 = performance.now();
    for (let i = 0; i < 1_000 || performance.now() - t0 < WARM_MS; i++) t.step(i);
    let i = 0;
    const r = heapScaling(2_000, () => t.step(i++));
    expect(isFlat(r), JSON.stringify(r)).toBe(true);
  } finally {
    t.dispose();
  }
}, 120_000);
