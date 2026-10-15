/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer's frames allocate nothing after warm-up (CS1a; seqspec:viewer#the-render-backend): pans and scrolls
// through the copy and the strips' draws at three zooms (every type), range labels sticking at the view's edges (two
// zooms), a wrapped section's scroll frames (letters, emphasis, each group label mode, the gutter's numbers with
// `none`), and whole frames and pans over 400,000 positions zoomed out (cells, values' columns, the ruler), on a
// recording backend that reports `Copy` and `Back` (the ops logged, not applied); the cases in composer-cases.ts,
// shared with the quick twin. The overlay's guards (hover frames with a selection, a highlight and row loci) are
// CS1b's. Each case is warmed and settled by count (`frameScaling`), then measured over `FRAME_N` frames and `4N`. The
// positive control runs the first case's pans allocating an object a frame (`leaking`) through the same warm-up,
// settle and measurement, which must catch it.
//
// Each case loads the viewer afresh (`vi.resetModules`), as RN1's draw guards do: in one process V8 can leave a draw in
// its middle tier after other documents deoptimized it, where vitest's module runner boxes doubles (RN1).

import { describe, expect, it } from 'vitest';
import {
  COPY_CELLS,
  copyPans,
  LABEL_CELLS,
  LABEL_MODES,
  NAVIGATORS,
  navigatorPans,
  rangeLabels,
  wrappedScroll,
  ZOOMED_OUT_CELLS,
  zoomedOut,
} from './composer-cases.js';
import {
  controlCaught,
  describeScaling,
  FRAME_N,
  FRAME_WARM,
  frameScaling,
  HEAP_DEEP,
  isFlat,
  leaking,
} from './heap.js';

const N = FRAME_N;
const guard = (frame: (i: number) => void) => {
  const r = frameScaling(frame, N);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
};

describe('composer frames allocate nothing after warm-up', () => {
  it.each(COPY_CELLS)(
    'pans and scrolls through the copy at cell %f (every type)',
    async (cell) => {
      const { frame, copies } = await copyPans(cell);
      const c0 = copies();
      guard(frame);
      expect(copies() - c0).toBeGreaterThan(2 * FRAME_WARM);
    },
    300_000,
  );

  it('the same pans allocating an object a frame are caught (positive control)', async () => {
    const { frame } = await copyPans(13);
    expect(controlCaught(() => frameScaling(leaking(frame), N))).toBe(true);
  }, 300_000);

  it.each(NAVIGATORS)(
    "every pass through a pan, the %s navigator's window lit",
    async (nav) => {
      const { frame } = await navigatorPans(nav);
      guard(frame);
    },
    300_000,
  );

  it.each(LABEL_MODES)(
    "a wrapped section's scroll frames with %s group labels",
    async (mode) => {
      const { frame } = await wrappedScroll(mode);
      guard(frame);
    },
    300_000,
  );

  it.each(LABEL_CELLS)(
    "pans with range labels centered, sticking at the view's edges and cut, at cell %f (seqspec:sources#levels-and-the-zoom-gate)",
    async (cell) => {
      const { frame, copies } = await rangeLabels(cell);
      const c0 = copies();
      guard(frame);
      expect(copies() - c0).toBeGreaterThan(2 * FRAME_WARM);
    },
    300_000,
  );

  it.each(HEAP_DEEP ? [0.0019, 0.03, 0.4] : ZOOMED_OUT_CELLS)(
    'whole frames and pans over 400,000 positions at cell %f',
    async (cell) => {
      const { frame } = await zoomedOut(cell);
      guard(frame);
    },
    300_000,
  );
});
