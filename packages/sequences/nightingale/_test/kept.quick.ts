/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of kept.heap.ts (`test:perf`): a host's highlight change projected and kept under a section's forty
// highlight runs, the guard's case (kept-cases.ts), with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from '../../core/_test/model/quick.js';
import { kept } from './kept-cases.js';

it("a section's highlights kept under a host's", async ({ task }) => {
  const c = await kept();
  quickCase(task, "kept: a section's highlights under a host's", c.step);
  expect(c.lit()).toBeGreaterThan(40 * 1_000);
});
