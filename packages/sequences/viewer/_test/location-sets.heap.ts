/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LS1's frames allocate nothing after warm-up: a section's `highlights` (a `band` set of 120 items in the main pass, a
// `highlight` set on the overlay) and a stepped bookmark's accent, through pans that copy and draw strips and the
// overlay, and a band set on a wrapped section (painted as a highlight there, LS2). The positive control: the pans
// allocating an object a frame (`leaking`), warmed and measured alike, are caught; a channel set on every frame (a
// rebuild per frame) is not flat either. Each case loads the viewer afresh (`vi.resetModules`), as the composer's
// guards do, and is warmed and settled by count (`frameScaling`, and `cheapScaling` for the wrapped band's cheaper
// frame); the cases are location-sets-cases.ts's, shared with the quick twin (location-sets.quick.ts; the rebuild per
// frame has none: it allocates by design).

import { describe, expect, it } from 'vitest';
import { cheapScaling, controlCaught, describeScaling, FRAME_N, frameScaling, isFlat, leaking } from './heap.js';
import { bandPans, pans, rebuilding, wrappedBand } from './location-sets-cases.js';

describe("location sets' frames", () => {
  it("pans with a band set, a highlight set and a bookmark's accent allocate nothing", async () => {
    const r = frameScaling(await bandPans(), FRAME_N);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
  }, 300_000);

  it('a band set in a wrapped section, painted as a highlight on the overlay (LS2), allocates nothing per frame', async () => {
    const r = cheapScaling(await wrappedBand(), FRAME_N);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
  }, 300_000);

  it('the same pans allocating an object a frame are caught (positive control)', async () => {
    const pan = await pans();
    expect(controlCaught(() => frameScaling(leaking(pan), FRAME_N))).toBe(true);
  }, 300_000);

  it('a channel set on every frame (a rebuild per frame) is not flat', async () => {
    const r = frameScaling(await rebuilding(), 2_000);
    expect(isFlat(r), describeScaling(r, 2_000)).toBe(false);
  }, 300_000);
});
