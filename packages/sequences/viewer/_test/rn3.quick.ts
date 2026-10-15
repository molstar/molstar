/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of rn3.heap.ts (`test:perf`): the same cases (rn3-cases.ts), each with its positive control.

import { describe, it } from 'vitest';
import { quickCase } from './heap.js';
import { RN3_CASES } from './rn3-cases.js';

describe("RN3's hot paths, quick", () => {
  it.for(RN3_CASES)('$name', async (c, { task }) => {
    quickCase(task, `rn3: ${c.id}`, await c.open());
  });
});
