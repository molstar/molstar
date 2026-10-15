/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The load path's cases at 3M residues, for its heap guard and its quick twin (molecule.heap.ts, molecule.quick.ts)
// (MD1; M2's optional check): `build` is one build of a long DNA string, without a short-lived string or array per
// residue (the old build upper-cased the string and looked each residue up through `String.fromCharCode`; the
// 128-entry tables read the codes straight into the molecule's array); `detect` is the alphabet's detection over the
// codes. Built when a case asks.

import { type Issue, Molecules, StringTable } from '../../src/model/index.js';
import { prng } from './support/generators.js';

export const SIZE = 3_000_000;

let text: string | undefined;
/** Mixed case, so that upper-casing has work to do. */
const dna = () => {
  if (text === undefined) {
    const rand = prng(11);
    text = Array.from({ length: SIZE }, () => 'ACGTacgt'[(rand() * 8) | 0]).join('');
  }
  return text;
};

/** One build of the molecule: the step, and what it built. */
export function buildCase() {
  const residues = dna();
  const strings = new StringTable();
  const issues: Issue[] = [];
  let length = 0;
  const step = () => {
    length = Molecules.build({ key: 'chr', residues }, 0, strings, issues, '/molecules/0')?.length ?? 0;
  };
  return { step, issues, length: () => length };
}

/** The alphabet's detection over the 3M codes: the step, and what it found. */
export function detectCase() {
  const code = Uint8Array.from(dna().toUpperCase(), (c) => c.charCodeAt(0));
  const found = { at: -1 };
  let alphabet = '';
  const step = () => {
    alphabet = Molecules.detectAlphabet(code, found);
  };
  return { step, alphabet: () => alphabet };
}
