/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the model's hot paths (M1, M2a, M2b; MD2's coordinate maps; MD3's indexes):
// runs of 250,000 and a million steps, the longer growing at most 64 KB beyond the shorter (vitest.heap.config.ts).
// Untimed: `heapScaling` warms each run by count (JIT, lazy state). The cases are hot-paths-cases.ts's, shared with the
// quick twin (hot-paths.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { FIXTURES, PATHS, sums } from './hot-paths-cases.js';

describe("heap guards: a step's growth does not scale with the run", () => {
  for (const [name, actual, expected] of FIXTURES) {
    it(name, () => {
      expect(actual()).toEqual(expected);
    });
  }

  for (const [name, n, step] of PATHS) {
    it(name, () => {
      const r = heapScaling(n, step);
      expect(isFlat(r), `${name}: ${describeScaling(r, n)}`).toBe(true);
      expect(Number.isFinite(sums())).toBe(true);
    });
  }

  // CM-05: the template's control, one guarded step with a leak, caught by the same measurement.
  it('the same step with a leak is caught (the positive control)', () => {
    const [, n, step] = PATHS[0];
    expect(controlCaught(() => heapScaling(n, leaking(step)))).toBe(true);
  });
});
