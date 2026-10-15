/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of hot-paths.heap.ts (`test:perf`): the model's hot paths (hot-paths-cases.ts, the guard's cases), each
// step with its positive control. `Loci.forEachLocation` is a walk of a million cells, a step of milliseconds, with the
// counts of one; the others are a few lookups.

import { describe, expect, it } from 'vitest';
import { FIXTURES, PATHS, sums } from './hot-paths-cases.js';
import { quickCase } from './quick.js';

/** A step of many milliseconds: few steps, in small windows. */
const WALK = { warmSteps: 20, steps: 5, windows: 3, maxWarmSteps: 200 };

describe("the model's hot paths, quick", () => {
  it('runs over the fixtures the guards describe', () => {
    for (const [name, actual, expected] of FIXTURES) expect(actual(), name).toEqual(expected);
  });

  it.for(PATHS.map(([name, n, step]) => ({ name, n, step })))('$name', async ({ name, n, step }, { task }) => {
    quickCase(task, `hot-paths: ${name}`, step, n < 1_000 ? WALK : {});
    expect(Number.isFinite(sums())).toBe(true);
  });
});
