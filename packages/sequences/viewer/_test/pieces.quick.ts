/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of pieces.heap.ts (`test:perf`): the conformance kit on the guard's reversed and ranged pieces
// (pieces-cases.ts, without the cases the trim dropped), each frame the kit draws measured in its `heap` hook, by the
// order the kit draws them in (`#1`, `#2`, ...: a layer, a zoom, the draw and below the row's minimum, a geometry),
// each with its control; and the overlay's hover frames over the scale and the rows on reversed pieces.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { CASES, KIT, reversedHover, run } from './pieces-cases.js';

describe("AX1: reversed and ranged pieces' frames, quick", () => {
  for (const [name, doc] of CASES)
    it(`${name} on ${doc}, reversed`, async ({ task }) => {
      let k = 0;
      const heap = (
        await run(name, doc, (frame) => {
          quickCase(task, `pieces: ${name} on ${doc}, reversed #${++k}`, frame, KIT);
        })
      ).filter((c) => c.check === 'heap');
      expect(heap.length).toBeGreaterThan(0);
      expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
    });
});

describe("AX1: the overlay's hover frames on reversed pieces, quick", () => {
  it.for(['scale', 'rows'] as const)('hover frames over %s, reversed', async (kind, { task }) => {
    const { frame, hovers } = await reversedHover(kind);
    quickCase(task, `pieces: hover frames over ${kind}, reversed`, frame);
    expect(hovers()).toBeGreaterThan(1_000);
  });
});
