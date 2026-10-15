/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The alignment build at 20k rows (MD1's Opus review): its residues are detected and checked by one counting pass
// over typed arrays, with no string or boxed value per residue. A per-residue allocation over 6M residues would be
// tens of MB; the build's per-alignment objects are a few kB. The build is one run (`heapGrowth`); the check is a
// repeated step (`heapScaling`). Untimed, as every core guard: both measurements warm by count. The cases are
// alignment-cases.ts's, shared with the quick twin (alignment.quick.ts).

import { describe, expect, it } from 'vitest';
import { buildCase, checkCase, LENGTH, ROWS } from './alignment-cases.js';
import { controlCaught, describeScaling, heapGrowth, heapScaling, isFlat, LIMIT, leak } from './heap.js';

describe('alignment build', () => {
  it('builds a 20k × 300 alignment, detecting its alphabet, without allocating per residue', () => {
    const { step, issues, alphabet } = buildCase();
    const growth = heapGrowth(step);
    expect(issues).toEqual([]);
    expect(alphabet()).toBe('protein');
    expect(growth).toBeLessThan(LIMIT);
    // The template's control: the fault it guards against, an object a residue, is caught by the same measurement.
    const cells = ROWS * LENGTH;
    expect(
      controlCaught(() =>
        heapGrowth(() => {
          step();
          for (let k = 0; k < cells; k++) leak(k);
        }),
      ),
    ).toBe(true);
  });

  it('checks 6M residues without allocating', () => {
    const { step, residues, found } = checkCase();
    const r = heapScaling(4, step);
    expect(residues.length).toBeGreaterThan(5_000_000);
    expect(found()).toEqual(['protein', -1, -1]);
    expect(isFlat(r), describeScaling(r, 4)).toBe(true);
    // The control: an object a residue in the same step is caught (a step is 6M residues: one leak a step is not
    // what this guard is for).
    const perResidue = () => {
      step();
      for (let k = 0; k < residues.length; k++) leak(k);
    };
    expect(controlCaught(() => heapScaling(4, perResidue))).toBe(true);
  });
});
