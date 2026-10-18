/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The overlay's hover frames allocate nothing after warm-up (CS1b): the pointer moving over 2 × 200,000 positions
// zoomed out with a selection and a highlight on them, and over an alignment's rows with a row's residues selected and
// cells highlighted (its row marks), each move an overlay frame on the recording backend: the channels' hover rebuilt,
// row marks, position marks per lane, the hovered row's band, renderers' overlays and the ring. Each frame is warmed
// and settled by count (`cheapScaling`, heap.ts), as composer.heap.ts's are. The positive control: the hover
// frames over the rows allocating an object each (`leaking`), through the same warm-up, settle and measurement, are
// caught.
//
// Each case loads the viewer afresh (`vi.resetModules`), as composer.heap.ts does; the cases are overlay-cases.ts's,
// shared with the quick twin (overlay.quick.ts).

import { describe, expect, it } from 'vitest';
import { cheapScaling, controlCaught, describeScaling, FRAME_WARM, isFlat, leaking } from './heap.js';
import { hovering, KINDS } from './overlay-cases.js';

const N = 1_000;

describe('overlay frames allocate nothing after warm-up', () => {
  it.each(KINDS)(
    'hover frames over %s with a selection and a highlight',
    async (kind) => {
      const { frame, hovers } = await hovering(kind);
      const r = cheapScaling(frame, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      // Not vacuous: the hovered element changed on most moves.
      expect(hovers()).toBeGreaterThan(FRAME_WARM);
    },
    300_000,
  );

  it('hover frames over rows allocating an object each are caught (positive control)', async () => {
    const { frame } = await hovering('rows');
    expect(controlCaught(() => cheapScaling(leaking(frame), N))).toBe(true);
  }, 300_000);
});
