/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guard for the label layer (NG-12), on RB0b's stepped `heapScaling`: `placeLabel` at a panned view each step
// (every zoom frame; the style writes, strings by nature, aside), its width, cut and height fractional (LA-13); the
// step is labels-heap-fixture.ts's `placing`, shared with the quick twin (labels-place.quick.ts). A file of its own.
// The positive control is the template's: the same step with a leak.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../../core/_test/model/heap.js';
import { N, placing } from './labels-heap-fixture.js';

describe('the label layer allocates nothing per frame', () => {
  it('placeLabel at a panned view each step', () => {
    const { step, modes } = placing();
    const r = heapScaling(N, step);
    expect([modes.full > 0, modes.cut > 0, modes.none > 0]).toEqual([true, true, true]);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });
});
