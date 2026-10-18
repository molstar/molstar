/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The overlay paths no other guard covers allocate nothing per frame after warm-up (seqspec:viewer#performance; VL-15):
// the vertical overview and scrollbar's overlay through a vertical scroll (`VerticalNav`, its row bars lit), the
// keyboard cursor's ring through a pan (`cursorRing`), a tree clade's hover (`cladeHover`), a held select drag's
// pending band on an alignment's rows and a dot plot's pending rectangle (`pending`, the moves dispatched as the engine
// sends them), each on the recording backend. Each case loads the viewer afresh (`vi.resetModules`) and is warmed and
// settled by count (`cheapScaling`, heap.ts). The positive control runs the vertical overview's step through
// the same warm-up, settle and measurement with `leaking` (the guard template). The cases are overlay-paths-cases.ts's,
// shared with the quick twin (overlay-paths.quick.ts).

import { describe, expect, it } from 'vitest';
import { cheapScaling, controlCaught, describeScaling, FRAME_N, isFlat, leaking } from './heap.js';
import { clade, type PathCase, overview, plotDrag, ring, selectDrag } from './overlay-paths-cases.js';

/** The guard: the case's step is flat over `FRAME_N` steps and `4 * FRAME_N`, after the warm-up by count. */
function guard(c: PathCase) {
  const r = cheapScaling(c.step, FRAME_N);
  expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
  c.held?.();
}

describe('overlay paths allocate nothing per frame after warm-up', () => {
  it("the vertical overview and scrollbar's overlay through a vertical scroll", async () => {
    guard(await overview());
  }, 300_000);

  it("the keyboard cursor's ring through a pan", async () => {
    guard(await ring());
  }, 300_000);

  it("a tree clade's hover, on and off", async () => {
    guard(await clade());
  }, 300_000);

  it("a held select drag's pending band on an alignment's rows", async () => {
    guard(await selectDrag());
  }, 300_000);

  it("a held drag's pending rectangle on a dot plot", async () => {
    guard(await plotDrag());
  }, 300_000);

  it("the positive control: the overview's step with a leak a step is caught", async () => {
    const { step } = await overview();
    expect(controlCaught(() => cheapScaling(leaking(step), FRAME_N))).toBe(true);
  }, 300_000);
});
