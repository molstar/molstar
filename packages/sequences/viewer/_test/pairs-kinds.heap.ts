/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames over cells that hold several kinds of contact allocate nothing after warm-up: the split cells' glyphs and "+n"
// markers through pans, and the overlay lighting the hovered glyph's kind (a hit and a hover per frame). The positive
// control runs the same frames allocating an object each (`leaking`) through the same warm-up and measurement, which
// must catch it; both are warmed and settled by count (`cheapScaling`, heap.ts). The case loads the viewer
// afresh (`vi.resetModules`), as the composer's guards do, and is pairs-kinds-cases.ts's, shared with the quick twin
// (pairs-kinds.quick.ts).

import { describe, expect, it } from 'vitest';
import { cheapScaling, controlCaught, describeScaling, FRAME_N, isFlat, leaking } from './heap.js';
import { kindsFrame } from './pairs-kinds-cases.js';

describe('frames over cells of several kinds', () => {
  it('pans, a hover on each glyph and its overlay allocate nothing; an object allocated per frame is caught', async () => {
    const { frame, multi, inSection } = await kindsFrame();
    expect(multi).toBe(true);
    expect(inSection).toBe(true);
    const r = cheapScaling(frame, FRAME_N);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(controlCaught(() => cheapScaling(leaking(frame), FRAME_N))).toBe(true);
  }, 300_000);
});
