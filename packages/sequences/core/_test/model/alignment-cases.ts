/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The alignment build's cases, for its heap guard and its quick twin (alignment.heap.ts, alignment.quick.ts), at 20k rows
// of 300 (MD1's Opus review): `build` is one build of the alignment, its residues detected and checked by one counting
// pass over typed arrays, with no string or boxed value per residue; `check` is that pass alone over the 6M residues
// the build keeps. Built when a case asks.

import { Alignments, type Issue, Molecules, type ResidueCheck, StringTable } from '../../src/model/index.js';
import { alignment } from './support/generators.js';

export const ROWS = 20_000;
export const LENGTH = 300;

let fixture: ReturnType<typeof alignment> | undefined;
const input = () => {
  fixture ??= alignment({ rows: ROWS, length: LENGTH, seed: 5 });
  return { key: 'big', flavor: 'aligned' as const, gapped: [...fixture.rows], labels: [...fixture.labels] };
};

/** One build of the alignment: the step, and what it found. */
export function buildCase() {
  const data = input();
  const strings = new StringTable();
  const issues: Issue[] = [];
  let alphabet: string | null = null;
  const step = () => {
    issues.length = 0;
    alphabet = Alignments.build(data, 0, strings, issues, '/alignments/0')?.alphabet ?? null;
  };
  return { step, issues, alphabet: () => alphabet };
}

/** The residue check over the built alignment's 6M residues: the step, and what it found. */
export function checkCase() {
  const strings = new StringTable();
  const a = Alignments.build(input(), 0, strings, [], '');
  const residues = a?.rows.residues ?? new Uint8Array(0);
  const check: ResidueCheck = { at: -1, invalid: -1 };
  let alphabet = '';
  const step = () => {
    alphabet = Molecules.checkResidues(residues, undefined, check);
  };
  return { step, residues, found: () => [alphabet, check.at, check.invalid] };
}
