/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of molecule.heap.ts (`test:perf`): the build of a 3M-residue DNA string and the detection of its
// alphabet (molecule-cases.ts, the guard's cases), each a step of milliseconds measured in bytes a run, each with its
// positive control.

import { describe, expect, it } from 'vitest';
import { buildCase, detectCase, SIZE } from './molecule-cases.js';
import { quickCase } from './quick.js';

/** A step of many milliseconds: few steps, in small windows. */
const RUN = { warmSteps: 20, steps: 5, windows: 3, maxWarmSteps: 200 };

describe('the molecule build, quick', () => {
  it('builds a 3M-residue DNA string', async ({ task }) => {
    const { step, issues, length } = buildCase();
    quickCase(task, 'molecule: build 3M residues', step, {
      ...RUN,
      control: false,
      why: "a build's bytes vary by tens a run (1,936 to 2,006): leak's 88 cannot be read beside them",
    });
    expect(issues).toEqual([]);
    expect(length()).toBe(SIZE);
  });

  it('detects the alphabet of 3M residues', async ({ task }) => {
    const { step, alphabet } = detectCase();
    quickCase(task, 'molecule: detect the alphabet of 3M residues', step, RUN);
    expect(alphabet()).toBe('dna');
  });
});
