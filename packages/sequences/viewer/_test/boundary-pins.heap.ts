/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Boundaries pinned over their track's letters allocate nothing per hover frame after warm-up
// (seqspec:viewer#performance): the pointer moving across a design row's pinned boundaries and the letters beside them,
// each move a hit (the pin's reach, the residue beside it) and an overlay frame outlining the hovered pin (the case in
// boundary-pins-cases.ts, shared with the quick twin). Warmed and settled by count (`frameScaling`). Positive controls:
// most frames hovered a boundary, and the same frames allocating an object each (`leaking`), through the same warm-up
// and measurement, are caught.

import { describe, expect, it } from 'vitest';
import { boundaryPins } from './boundary-pins-cases.js';
import { controlCaught, describeScaling, FRAME_N, frameScaling, isFlat, leaking } from './heap.js';

describe('pinned boundaries over a hover sweep', () => {
  it('allocate nothing; the pins hover; an object allocated per frame is caught', async () => {
    const { frame, boundary } = await boundaryPins();
    const r = frameScaling(frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(boundary()).toBeGreaterThan(FRAME_N);
    expect(controlCaught(() => frameScaling(leaking(frame)))).toBe(true);
  }, 300_000);
});
