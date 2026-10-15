/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ideogram navigator's heap guard (ST5; the memory rule): its strip, redrawn on every frame here (the nav pass is
// drawn by hand, as a layout or skin change would), and the window moved over it allocate nothing once warm. Over
// hg38's committed cytobands on the streamed genome's first three chromosomes (about 200 bands; the case in
// ideogram-cases.ts, shared with the quick twin). Warmed and settled by count (`frameScaling`: the strip's frame climbs
// to TurboFan slowly, as the streamed items' do, and the settle runs until it stops allocating). The positive control:
// the same frames allocating an object each (`leaking`), through the same warm-up and measurement, are caught.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { ideogramFrame } from './ideogram-cases.js';

const N = 1_000;
const guard = (frame: () => void) => {
  const r = frameScaling(frame, N);
  return isFlat(r) ? '' : describeScaling(r, N);
};

describe('the ideogram navigator allocates nothing per frame', () => {
  it('its strip and the window over it, three chromosomes', async () => {
    const g = await ideogramFrame();
    g.frame();
    // Not vacuous: the bands, the centromeres' triangles and the window.
    expect(g.backend.rectItems.n).toBeGreaterThan(80);
    expect(g.backend.pathItems.nv).toBeGreaterThanOrEqual(24);
    expect(guard(g.frame)).toBe('');
    g.dispose();
  }, 300_000);

  it('the same frames allocating an object each are caught (positive control)', async () => {
    const g = await ideogramFrame();
    expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
    g.dispose();
  }, 300_000);
});
