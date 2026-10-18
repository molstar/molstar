/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of split-guides.heap.ts (`test:perf`): whole frames of case 13's variants grid on a pan at DPR 2, the
// guides drawn (the guard's case, split-guides-cases.ts), with its control. Not vacuous: most frames, counted over the
// whole run, drew all 21 guides.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { guidesFrame } from './split-guides-cases.js';

describe("a split track's row guides, quick", () => {
  it("whole frames of case 13's variants grid on a pan, the guides drawn", async ({ task }) => {
    const { frame, full, guides } = await guidesFrame();
    let calls = 0;
    quickCase(task, 'split-guides: variants grid on a pan', (i) => {
      frame(i);
      calls++;
    });
    expect(full()).toBeGreaterThan(0.9 * calls);
    expect(guides()).toBeGreaterThan(0);
  });
});
