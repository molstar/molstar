/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of letter-maps.heap.ts (`test:perf`): the conformance kit on the guard's cases (letter-maps-cases.ts,
// the reverse complement on the same reduced set), each frame the kit draws measured in its `heap` hook, by the order
// the kit draws them in (`#1`, `#2`, ...: a layer, a zoom, the draw and below the row's minimum, a geometry), each
// with its control. The positive control case of the guard is the control of every frame here.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { CASES, KIT, REVERSED, run } from './letter-maps-cases.js';

describe('LM1: frames drawn through a letter map, quick', () => {
  for (const reversed of [false, true])
    for (const [name, doc] of CASES)
      if (!reversed || REVERSED.has(`${name} on ${doc}`)) {
        const what = `${name} on ${doc}, ${reversed ? 'as a reverse complement' : 'complemented'}`;
        it(what, async ({ task }) => {
          let k = 0;
          const heap = (
            await run(name, doc, reversed, (frame) => {
              quickCase(task, `letter-maps: ${what} #${++k}`, frame, KIT);
            })
          ).filter((c) => c.check === 'heap');
          expect(heap.length).toBeGreaterThan(0);
          expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
        });
      }
});
