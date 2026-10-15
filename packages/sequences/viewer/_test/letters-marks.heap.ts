/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A frame of a text-only letters row's marks (dots, then bars), another polymer's (`of`) and the plain row's under an
// empty `rules` theme, allocates nothing after warm-up (seqspec:viewer#the-render-backend); each case loads the viewer
// afresh and is warmed by count (`cheapScaling`, heap.ts: the frame is a few microseconds). The positive
// control: the design row's dots allocating an object a frame (`leaking`), through the same warm-up and measurement,
// are caught. The cases are letters-marks-cases.ts's, shared with the quick twin (letters-marks.quick.ts).

import { describe, expect, it } from 'vitest';
import { cheapScaling, controlCaught, describeScaling, FRAME_N, isFlat, leaking } from './heap.js';
import { LABELS, MARKS, marksFrame } from './letters-marks-cases.js';

describe("a text-only letters row's marks allocate nothing", () => {
  // The plain row (`blockMarks`): its 48 bytes a frame under a short warm-up were the draw's climb to TurboFan.
  for (const label of LABELS)
    for (const [cell, dpr, what] of MARKS)
      it(`${label}: ${what}`, async () => {
        const c = await marksFrame(label, cell, dpr);
        const r = cheapScaling(c.frame, FRAME_N);
        c.dispose();
        // Not vacuous: the frame draws the marks.
        expect(c.marked).toBeGreaterThan(0);
        expect(isFlat(r), describeScaling(r, FRAME_N)).toBe(true);
      }, 300_000);

  it('Design: dots, an object allocated a frame (the positive control)', async () => {
    const c = await marksFrame('Design', 5, 2);
    const caught = controlCaught(() => cheapScaling(leaking(c.frame), FRAME_N));
    c.dispose();
    expect(caught).toBe(true);
  }, 300_000);
});
