/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of composer.heap.ts (`test:perf`): the same cases (composer-cases.ts), a frame a step, each on a
// viewer loaded afresh with its positive control: the copy's pans at three zooms, the navigators' pans, the wrapped
// scroll under each group label mode, the range labels' pans at two zooms, and the frames over 400,000 positions zoomed
// out (`FRAME` counts; the navigators `NAVIGATOR`, over 400,000 positions `ZOOMED_OUT_COUNTS`).

import { describe, expect, it } from 'vitest';
import {
  COPY_CELLS,
  copyPans,
  FRAME,
  LABEL_CELLS,
  LABEL_MODES,
  NAVIGATOR,
  NAVIGATORS,
  navigatorPans,
  rangeLabels,
  wrappedScroll,
  ZOOMED_OUT_CELLS,
  ZOOMED_OUT_COUNTS,
  zoomedOut,
} from './composer-cases.js';
import { quickCase } from './heap.js';

describe('composer frames, quick', () => {
  it.for(COPY_CELLS)('pans and scrolls through the copy at cell %f', async (cell, { task }) => {
    const { frame, copies } = await copyPans(cell);
    const c0 = copies();
    quickCase(task, `composer: copy pans at ${cell}`, frame, FRAME);
    expect(copies() - c0).toBeGreaterThan(1_000);
  });

  it.for(NAVIGATORS)("every pass through a pan, the %s navigator's window lit", async (nav, { task }) => {
    const { frame } = await navigatorPans(nav);
    quickCase(task, `composer: ${nav} navigator pans`, frame, NAVIGATOR);
  });

  it.for(LABEL_MODES)("a wrapped section's scroll frames with %s group labels", async (mode, { task }) => {
    const { frame } = await wrappedScroll(mode);
    quickCase(task, `composer: wrapped scroll, ${mode} group labels`, frame, FRAME);
  });

  it.for(LABEL_CELLS)('pans with range labels at cell %f', async (cell, { task }) => {
    const { frame, copies } = await rangeLabels(cell);
    const c0 = copies();
    quickCase(task, `composer: range labels at ${cell}`, frame, FRAME);
    expect(copies() - c0).toBeGreaterThan(1_000);
  });

  it.for(ZOOMED_OUT_CELLS)('whole frames and pans over 400,000 positions at cell %f', async (cell, { task }) => {
    const { frame } = await zoomedOut(cell);
    quickCase(task, `composer: 400,000 positions at ${cell}`, frame, ZOOMED_OUT_COUNTS[cell]);
  });
});
