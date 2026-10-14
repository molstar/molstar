/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CG3's ring window and labels allocate nothing per frame or pointer event after warm-up
// (seqspec:interaction#linked-sections): whole ring frames with the window on the linked linear section, that section
// panning and zooming each step, the pointer moving round the rings and over the labels (polar hits, the labels' hits,
// the cursor's window part). The window's drag is ring-window-drag.heap.ts (in one file, the drag allocated after these
// cases had run: their shared feedback). The positive control: the same steps allocating an object each (`leaking`),
// through the same warm-up, settle and measurement, are caught.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, isFlat, leaking } from './heap.js';
import { measure, run, setup } from './ring-window-heap-run.js';

/** A deep guard's steps: the trimmed run's 500 catch about 43 bytes a step, these about 22. */
const STEPS = 1_000;

describe("CG3's ring window and labels allocate nothing after warm-up", () => {
  it.each([600, 120_000])(
    'ring frames with the window, the linear section moving, the pointer round the rings and labels, %i positions',
    async (n) => {
      const { r, sectors, labelHits } = await run('frames', n, STEPS);
      expect(sectors).toBeGreaterThan(0);
      expect(labelHits).toBeGreaterThan(0);
      expect(isFlat(r), describeScaling(r, STEPS)).toBe(true);
    },
    120_000,
  );

  it('the same steps allocating an object each are caught (positive control)', async () => {
    const { step } = await setup('frames', 600);
    expect(controlCaught(() => measure(leaking(step), STEPS))).toBe(true);
  }, 120_000);
});
