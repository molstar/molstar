/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CG3's ring window drag allocates nothing per pointer event after warm-up (seqspec:interaction#linked-sections): a press
// on the window's body, then moves along it, each setting the linked linear section through the link and drawing the
// ring's overlay. Its own file since a timed warm-up ended inside the drag's climb to TurboFan once
// ring-window.heap.ts's frame cases had run (the "shared feedback" it was blamed on); warmed and settled by count now
// (`frameScaling`). The case is ring-window-heap-run.ts', shared with the quick twin (ring-window-drag.quick.ts). The
// positive control: the same moves allocating an object each (`leaking`), through the same warm-up, settle and
// measurement, are caught.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, isFlat, leaking } from './heap.js';
import { measure, run, setup } from './ring-window-heap-run.js';

/** A deep guard's steps: the trimmed run's 500 catch about 43 bytes a step, these about 22 (the drag once cost 32). */
const STEPS = 1_000;

describe("CG3's ring window drag allocates nothing after warm-up", () => {
  it('a drag of the window, each move setting the linear section', async () => {
    const { r, owner } = await run('drag', 600, STEPS);
    expect(owner).toBe('ringWindow');
    expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
  }, 300_000);

  it('the same moves allocating an object each are caught (positive control)', async () => {
    const s = await setup('drag', 600);
    expect(controlCaught(() => measure(leaking(s.step), STEPS))).toBe(true);
    expect(s.owner()).toBe('ringWindow');
  }, 300_000);
});
