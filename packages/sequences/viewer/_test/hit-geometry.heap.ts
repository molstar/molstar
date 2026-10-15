/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hover frames over sized shapes allocate nothing after warm-up (HR2): the pointer moving over dots, bands, boundaries'
// glyphs (the end boundary's too, as ranges, points and stacked variants, the pointer around its edge), pair ends and
// wide arcs, each move a hit (the layer's themes set for its `hit`, the glyphs near the point measured), the hover
// rebuilt and an overlay frame on the recording backend; and over an own unit's row in `UI.OwnRowHover` "residue" mode
// (the cases in hit-geometry-cases.ts, shared with the quick twin). Warmed and settled by count (`frameScaling`). The
// positive control runs the first case's frames allocating an object each (`leaking`) through the same warm-up and
// measurement, which must catch it. Each case loads the viewer afresh.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, FRAME_WARM, frameScaling, isFlat, leaking } from './heap.js';
import { ownRow, sizedShapes } from './hit-geometry-cases.js';

describe('hover frames over sized shapes', () => {
  it('allocate nothing; an object allocated per frame is caught', async () => {
    const { frame, hovers } = await sizedShapes();
    const r = frameScaling(frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(hovers()).toBeGreaterThan(FRAME_WARM);
    expect(controlCaught(() => frameScaling(leaking(frame)))).toBe(true);
  }, 300_000);

  it('"residue" mode over an own unit\'s row allocates nothing', async () => {
    const { frame } = await ownRow();
    const r = frameScaling(frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
  }, 300_000);
});
