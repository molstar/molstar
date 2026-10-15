/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Carets allocate nothing per hover frame after warm-up (seqspec:interaction#channels, seqspec:viewer#performance):
// the pointer sweeping the design's residues, its inserts among them, each move a hover whose projection onto the wild
// type's section finds the insert's next pair and lights a caret; then the wild type's 29–31, which the design deletes
// (the forward leg's insert: a caret at the design's junction), and a point just past the wild type's block (the end
// reach, its boundary marks asked); an overlay frame each, with a selected insert drawn as a caret every frame (the
// case in insert-carets-cases.ts, shared with the quick twin). Warmed and settled by count (`frameScaling`).
// Positive controls: most frames drew a hover's caret beside the selection's, and the same frames allocating an
// object each (`leaking`), through the same warm-up and measurement, are caught.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, frameScaling, isFlat, leaking } from './heap.js';
import { insertCarets } from './insert-carets-cases.js';

describe('carets over a hover sweep', () => {
  it('allocate nothing; most frames draw a caret; an object allocated per frame is caught', async () => {
    const { frame, share } = await insertCarets();
    const r = frameScaling(frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    // 8 of the 12 points light a caret (the inserts, the deleted stretch).
    expect(share()).toBeGreaterThan(0.5);
    expect(controlCaught(() => frameScaling(leaking(frame)))).toBe(true);
  }, 300_000);
});
