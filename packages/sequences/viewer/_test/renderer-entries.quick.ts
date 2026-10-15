/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of renderer-entries.heap.ts (`test:perf`): a hit test a move over the kernel entries' own hits
// (renderer-entries-cases.ts, the guard's case), with its control.

import { expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { entryHits } from './renderer-entries-cases.js';

it("hit tests per move over the kernel entries' hits, quick", async ({ task }) => {
  const { step, items, states } = await entryHits();
  expect([items > 20, states > 20]).toEqual([true, true]);
  quickCase(task, 'renderer-entries: hits', step);
});
