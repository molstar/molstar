/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Rings allocate nothing per frame after warm-up (CG1a; seqspec:viewer#the-render-backend): whole frames of a circular
// section (main and overlay passes through the polar map: every type's batches mapped to sectors, paths, glyphs and
// text) with the pointer moving round the rings (polar hits and the hover, the hover label's strings left out as
// overlay-fixtures.ts does), over 240 positions and over 120,000 (reductions per angular column). CG1b: with a
// selection and a highlight across the origin (marks as sectors), the keyboard cursor's ring (its rect mapped each
// frame), the hover as a radial line or the hovered item's sector and outline, and the ruler's minor ticks. The
// positive control: the same steps allocating an object each (`leaking`), through the same warm-up, settle and
// measurement, are caught.
//
// Each case loads the viewer afresh (`vi.resetModules`; the case in ring-cases.ts, shared with the quick twin,
// ring.quick.ts) and is warmed and settled by count (`cheapScaling`, heap.ts: a frame of tens of microseconds, whose climb under V8's defaults has quiet spans of thousands of steps): over 120,000 positions the climb to
// TurboFan runs thousands of frames at 20 to 210 bytes a frame, which the settle waits out.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, cheapScaling, isFlat, leaking } from './heap.js';
import { ring } from './ring-cases.js';

describe('ring frames allocate nothing after warm-up', () => {
  it.each([240, 120_000])(
    'whole frames with polar hovers round the rings, %i positions',
    async (n) => {
      const { step, be } = await ring(n);
      const r = cheapScaling(step);
      expect(be.sectorItems.n).toBeGreaterThan(0);
      expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    },
    300_000,
  );

  it('the same steps allocating an object each are caught (positive control)', async () => {
    const { step } = await ring(240);
    expect(controlCaught(() => cheapScaling(leaking(step)))).toBe(true);
  }, 300_000);
});
