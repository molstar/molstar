/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The standard chemical component of each letter by alphabet (A is ALA in protein, DA in DNA, A in RNA), and its
// inverse: what a position of a sequence is when nothing says otherwise. The inverse half of the one-letter codes
// (`@molstar/model`'s `getProteinOneLetterCode` and friends), which also read force-field and modified names; these
// tables hold the standard components only, so a letter has one component and a component one letter.

/** The alphabets of a sequence's letters. */
export type SequenceAlphabet = 'protein' | 'dna' | 'rna';

/**
 * The standard component (a CCD code) of each uppercase letter, by alphabet. The unknown letters have the CCD's unknown
 * components: `X` is UNK in protein, `N` is DN in DNA and N in RNA. Ambiguity codes (B, Z, J; R, Y and the like) and the
 * gap have none.
 */
export const StandardComponents: { readonly [A in SequenceAlphabet]: Readonly<Record<string, string>> } = {
  protein: {
    A: 'ALA',
    R: 'ARG',
    N: 'ASN',
    D: 'ASP',
    C: 'CYS',
    Q: 'GLN',
    E: 'GLU',
    G: 'GLY',
    H: 'HIS',
    I: 'ILE',
    L: 'LEU',
    K: 'LYS',
    M: 'MET',
    F: 'PHE',
    P: 'PRO',
    S: 'SER',
    T: 'THR',
    W: 'TRP',
    Y: 'TYR',
    V: 'VAL',
    U: 'SEC',
    O: 'PYL',
    X: 'UNK',
  },
  dna: { A: 'DA', C: 'DC', G: 'DG', T: 'DT', N: 'DN' },
  rna: { A: 'A', C: 'C', G: 'G', U: 'U', N: 'N' },
};

const Letters: { readonly [A in SequenceAlphabet]: ReadonlyMap<string, string> } = {
  protein: inverse(StandardComponents.protein),
  dna: inverse(StandardComponents.dna),
  rna: inverse(StandardComponents.rna),
};

function inverse(table: Readonly<Record<string, string>>): ReadonlyMap<string, string> {
  return new Map(Object.entries(table).map(([letter, component]) => [component, letter]));
}

/** The standard component of a letter (either case) in an alphabet; undefined for a letter that has none. */
export function standardComponent(alphabet: SequenceAlphabet, letter: string): string | undefined {
  const table = StandardComponents[alphabet];
  const upper = letter.toUpperCase();
  return Object.hasOwn(table, upper) ? table[upper] : undefined;
}

/** The letter whose standard component a code is, in an alphabet; undefined for a code that is no standard one there. */
export function standardLetter(alphabet: SequenceAlphabet, component: string): string | undefined {
  return Letters[alphabet].get(component.toUpperCase());
}
