/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Frames of arcs whose partners lie beyond the window allocate nothing after warm-up: case 7 zoomed into barstar, every
// contact's partner off-window in barnase, so each arc is pushed as its pieces in view (cut at the window and the
// band's far edge) with a wedge where it meets the edge; pans move the view each frame (`Arcs.frame` recomputed), and
// a hover on an end lights its arcs on the overlay (the case in arcs-off-window-cases.ts, shared with the quick twin).
// Warmed and settled by count (`frameScaling`). The positive control runs the same frames with an object allocated
// each (`leaking`) through the same warm-up and measurement, which must catch it. The case loads the viewer afresh.

import { describe, expect, it } from 'vitest';
import { arcsOffWindow } from './arcs-off-window-cases.js';
import { controlCaught, describeScaling, FRAME_N, frameScaling, isFlat, leaking } from './heap.js';

describe('frames of arcs whose partners are beyond the window', () => {
  it('pans, a hover on an end and its overlay allocate nothing; an object allocated per frame is caught', async () => {
    const { frame, glyphs } = await arcsOffWindow();
    expect(glyphs.length).toBeGreaterThan(5);
    const r = frameScaling(frame);
    expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
    expect(controlCaught(() => frameScaling(leaking(frame)))).toBe(true);
  }, 300_000);
});
