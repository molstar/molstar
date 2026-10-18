/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames of layers drawn through size themes allocate nothing after warm-up: pairs stroked at a width theme's (sorted
// by color then width) with an end radius, ranges at a thickness, points at a radius, values at a thickness, through
// pans. The positive control runs the same frames allocating an object each (`leaking`) through the same warm-up and
// measurement, which must catch it. Both are warmed and settled by count (`frameScaling`). The case loads the viewer
// afresh (`vi.resetModules`), as the composer's guards do, and is size-channels-cases.ts's, shared with the quick twin
// (size-channels.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, frameScaling, isFlat, leaking } from './heap.js';
import { sizeFrame } from './size-channels-cases.js';

describe('frames of layers under size themes', () => {
  it('pans allocate nothing; an object allocated per frame is caught', async () => {
    const frame = await sizeFrame();
    const r = frameScaling(frame, FRAME_N);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(controlCaught(() => frameScaling(leaking(frame), FRAME_N))).toBe(true);
  }, 300_000);
});
