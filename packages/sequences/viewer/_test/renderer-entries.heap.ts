/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guard for LY3b's hits (seqspec:viewer#performance): a hit test per pointer move over the kernel entries' own
// hits (letters `of` another polymer, ranges in lanes, pairs' arcs and ends over hundreds of pairs, a matrix's state
// bands), on `heapScaling` (a step one move; the case in renderer-entries-cases.ts, shared with the quick twin,
// renderer-entries.quick.ts). The positive control: the same moves allocating an object each (`leaking`), through the
// same measurement, are caught.

import { expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { entryHits } from './renderer-entries-cases.js';

const N = 50_000;

it("hit tests per move over the kernel entries' hits allocate nothing", async () => {
  const { step, items, states } = await entryHits();
  // Not vacuous: items and matrix states are hit.
  expect([items > 20, states > 20]).toEqual([true, true]);
  const r = heapScaling(N, step);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
  expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
});
