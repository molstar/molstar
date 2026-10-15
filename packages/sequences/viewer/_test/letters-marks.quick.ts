/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of letters-marks.heap.ts (`test:perf`): the design row's and the plain row's letters layer drawn as
// dots and as bars (the guard's cases, letters-marks-cases.ts), each with its control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { LABELS, MARKS, marksFrame } from './letters-marks-cases.js';

describe("a text-only letters row's marks, quick", () => {
  for (const label of LABELS)
    for (const [cell, dpr, what] of MARKS)
      it(`${label}: ${what}`, async ({ task }) => {
        const c = await marksFrame(label, cell, dpr);
        expect(c.marked).toBeGreaterThan(0);
        quickCase(task, `letters-marks: ${label}: ${what}`, c.frame);
        c.dispose();
      });
});
