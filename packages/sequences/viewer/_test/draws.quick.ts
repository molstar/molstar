/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of draws.heap.ts (`test:perf`): the same kit runs (draws-cases.ts), a renderer on a document each,
// one quick case per geometry with its positive control: the geometry's first full draw, its first layer at the fit
// (`DRAW` counts; ids by kit-cases.ts `frameIds`, `0` that draw). Every frame of the deep guard (130: every layer,
// zoom and `tiny`) took 93 s as quick cases: the closer zooms, `tiny` and the further layers are the deep guard's.

import { describe, expect, it } from 'vitest';
import { CASES, DRAW, drawsKit } from './draws-cases.js';
import { quickCase } from './heap.js';
import { frameIds } from './kit-cases.js';

describe("RN1's frames, quick", () => {
  for (const [name, docs] of Object.entries(CASES))
    for (const doc of docs)
      it(`${name} on ${doc}`, async ({ task }) => {
        const id = frameIds(`draws: ${name} on ${doc} at the fit`);
        const heap = (await drawsKit(name, doc, (frame, at) => quickCase(task, id(at), frame, DRAW), true)).filter(
          (c) => c.check === 'heap',
        );
        expect(heap.length).toBeGreaterThan(0);
        expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
      });
});
