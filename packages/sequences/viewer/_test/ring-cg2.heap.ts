/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CG2's ring draws allocate nothing per frame after warm-up (seqspec:viewer#the-render-backend): whole frames of a
// circular section with items' labels outside the ring (their leaders and texts from the kept placement), arrowheads by
// a direction field and pairs as chords, the pointer moving round the rings (polar hits and the hover). The positive
// control: the same steps allocating an object each (`leaking`), through the same warm-up, settle and measurement, are
// caught. Each case loads the viewer afresh (`vi.resetModules`; the case in ring-cg2-cases.ts, shared with the quick
// twin, ring-cg2.quick.ts) and is warmed and settled by count (`cheapScaling`, heap.ts: a frame of tens of microseconds, whose climb under V8's defaults has quiet spans of thousands of steps).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, cheapScaling, isFlat, leaking } from './heap.js';
import { ringCg2 } from './ring-cg2-cases.js';

describe("CG2's ring frames allocate nothing after warm-up", () => {
  it.each([240, 120_000])(
    'labels with leaders, arrowheads and chords, polar hovers round the rings, %i positions',
    async (n) => {
      const { step, be } = await ringCg2(n);
      const r = cheapScaling(step);
      expect(be.textItems.n).toBeGreaterThan(0);
      expect(be.pathItems.n).toBeGreaterThan(0);
      expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    },
    300_000,
  );

  it('the same steps allocating an object each are caught (positive control)', async () => {
    const { step } = await ringCg2(240);
    expect(controlCaught(() => cheapScaling(leaking(step)))).toBe(true);
  }, 300_000);
});
