/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of matrix-dense.heap.ts (`test:perf`): the dense matrix's pointer move and whole frame at DPR 1.5 and
// 1, and two stacked blocks of one state at DPR 4/3 (the guard's cases, matrix-dense-cases.ts), each with its control.
// Not vacuous: tiles are drawn as images and none is built once warmed (counted from step 5,000, inside the warm-up).

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { DENSE, setup } from './matrix-dense-cases.js';

/**
 * A step of 0.2 ms whose pointer and zoom alternate, in 2,000-step windows. At a fractional DPR it allocates 16 bytes a
 * frame once optimized (the passes' `dpr` argument boxed: packages/sequences/todo.md), so it never reads quiet: a cap.
 */
const STEP = { warmSteps: 2_000, steps: 2_000, windows: 3, maxWarmSteps: 50_000 };

describe("KS2's dense matrix frames, quick", () => {
  for (const [dpr, stacked] of DENSE)
    it(`${stacked ? 'two stacked blocks' : 'the per-device-row tiles'} at DPR ${dpr}`, async ({ task }) => {
      const { step, visits } = await setup(dpr, stacked);
      let images = 0;
      let built = 0;
      quickCase(
        task,
        `matrix-dense: ${stacked ? 'stacked' : 'tiles'}, DPR ${dpr}`,
        (i) => {
          step(i);
          if (i === 5_000) {
            images = visits.images;
            built = visits.tiles;
          }
        },
        STEP,
      );
      expect(visits.images - images).toBeGreaterThan(0);
      expect(visits.tiles - built).toBe(0);
    });
});
