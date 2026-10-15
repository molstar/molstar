/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of alignment.heap.ts (`test:perf`): the build of a 20k × 300 alignment and the check of its 6M residues
// (alignment-cases.ts, the guard's cases), each a step of many milliseconds measured in bytes a run, each with its
// positive control. The build's bytes a run are its own objects (the residues' arrays are off the JS heap): none a
// residue.

import { describe, expect, it } from 'vitest';
import { buildCase, checkCase } from './alignment-cases.js';
import { quickCase } from './quick.js';

/** A step of many milliseconds: few steps, in small windows. */
const RUN = { warmSteps: 20, steps: 5, windows: 3, maxWarmSteps: 200 };

describe('the alignment build, quick', () => {
  it('builds a 20k × 300 alignment', async ({ task }) => {
    const { step, issues, alphabet } = buildCase();
    quickCase(task, 'alignment: build 20k x 300', step, RUN);
    expect(issues).toEqual([]);
    expect(alphabet()).toBe('protein');
  });

  it('checks 6M residues', async ({ task }) => {
    const { step, residues, found } = checkCase();
    quickCase(task, 'alignment: check 6M residues', step, RUN);
    expect(residues.length).toBeGreaterThan(5_000_000);
    expect(found()).toEqual(['protein', -1, -1]);
  });
});
