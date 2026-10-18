/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The load path at 3M residues (MD1; M2's optional check): a long DNA string builds without a short-lived string or
// array per residue. The old build upper-cased the string (3M characters) and looked each residue up through
// `String.fromCharCode`; the 128-entry tables read the codes straight into the molecule's array. The build is one run
// (`heapGrowth`); detection is a repeated step (`heapScaling`). Untimed, as every core guard: both measurements warm
// by count. The cases are molecule-cases.ts's, shared with the quick twin (molecule.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapGrowth, heapScaling, isFlat, LIMIT, leak } from './heap.js';
import { buildCase, detectCase, SIZE } from './molecule-cases.js';

describe('molecule build', () => {
  it('builds a 3M-residue DNA string without allocating beyond its own array', () => {
    const { step, issues, length } = buildCase();
    const growth = heapGrowth(step);
    expect(issues).toEqual([]);
    expect(length()).toBe(SIZE);
    expect(growth).toBeLessThan(LIMIT);
    // The template's control: the old build's fault, a short-lived object a residue, is caught by the same measurement.
    const perResidue = () => {
      step();
      for (let k = 0; k < SIZE; k++) leak(k);
    };
    expect(controlCaught(() => heapGrowth(perResidue))).toBe(true);
  });

  it('detects the alphabet of 3M residues without allocating', () => {
    const { step, alphabet } = detectCase();
    const r = heapScaling(4, step);
    expect(alphabet()).toBe('dna');
    expect(isFlat(r), describeScaling(r, 4)).toBe(true);
    // The control: an object a residue in the same step is caught.
    const perResidue = () => {
      step();
      for (let k = 0; k < SIZE; k++) leak(k);
    };
    expect(controlCaught(() => heapScaling(4, perResidue))).toBe(true);
  });
});
