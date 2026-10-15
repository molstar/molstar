/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// EN2b's channels in frames allocate nothing after warm-up: pairs at a height theme's apex with their labels, variants
// and secondary structure at a thickness, through pans. The positive control runs the same frames allocating an object
// each (`leaking`) through the same warm-up and measurement, which must catch it. Both are warmed and settled by count
// (`frameScaling`). The case loads the viewer afresh (`vi.resetModules`), as the composer's guards do, and is
// size-channels-en2b-cases.ts's, shared with the quick twin (size-channels-en2b.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, frameScaling, isFlat, leaking } from './heap.js';
import { en2bFrame } from './size-channels-en2b-cases.js';

describe("frames of EN2b's size channels", () => {
  it('pans allocate nothing; an object allocated per frame is caught', async () => {
    const frame = await en2bFrame();
    const r = frameScaling(frame, FRAME_N);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(controlCaught(() => frameScaling(leaking(frame), FRAME_N))).toBe(true);
  }, 300_000);
});
