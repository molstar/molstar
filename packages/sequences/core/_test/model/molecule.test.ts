/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Molecules (seqspec:model#molecules). The old stack tested them through the set (set.test: "molecules"); here
// `Molecules.build` is called directly, at the path the set gives it. MD1's new cases: M2's input and an empty polymer.

import { describe, expect, it } from 'vitest';
import { type Issue, type MoleculeInput, MoleculeKind, Molecules, StringTable } from '../../src/model/index.js';

const ascii = (s: string) => Uint8Array.from(s, (c) => c.charCodeAt(0));
const codes = (issues: Issue[]) => issues.map((i) => `${i.severity} ${i.code} ${i.path}`);

function build(input: MoleculeInput) {
  const strings = new StringTable();
  const issues: Issue[] = [];
  const molecule = Molecules.build(input, 0, strings, issues, '/molecules/0');
  return { molecule, issues, strings };
}
const residuesOf = (input: MoleculeInput) => {
  const { molecule } = build(input);
  return molecule && String.fromCharCode(...molecule.code);
};

describe('alphabets', () => {
  it('detects them; HHHHHH is protein', () => {
    const detect = (s: string) => Molecules.detectAlphabet(ascii(s));
    expect(detect('HHHHHH')).toBe('protein');
    expect(detect('ACGTTGCAN')).toBe('dna');
    expect(detect('ACGUUGCA')).toBe('rna');
    expect(detect('MKTAYIAKQR')).toBe('protein');
    expect(detect('ACGTACGTEE')).toBe('protein');
    // 90% are ACGTUN but E is no nucleotide code, so this is protein (it was DNA, then rejected on the E).
    expect(detect('ACGTACGTAE')).toBe('protein');
  });

  it('counts an ambiguity code as a letter that is not ACGTUN (peptides of ambiguity letters)', () => {
    const detect = (s: string) => Molecules.detectAlphabet(ascii(s));
    // His tags and peptides of ambiguity letters (H, M, K, R, S, …) plus a few ACGT: with the ambiguity codes left out
    // of the count, the few ACGT were 100% and these loaded as DNA, silently.
    for (const peptide of ['MGSSHHHHHHSSG', 'HHHHHHA', 'MKHDSWRYVA', 'ACGTRYKMNNACGT']) {
      expect(detect(peptide), peptide).toBe('protein');
    }
    // A nucleotide sequence with a few ambiguity codes is still one: 90% exactly is enough.
    expect(detect('ACGTACGTAR')).toBe('dna');
    expect(detect('ACGTACGTRY')).toBe('protein');
    expect(detect('ACGTACGTACGTACGTACGRYN')).toBe('dna');
  });

  it('reads a peptide that looks like nucleotides as protein, with a warning at the letter', () => {
    const { molecule, issues } = build({ key: 'tag', residues: 'AAAAAAAAAGE' });
    expect(molecule?.alphabet).toBe('protein');
    expect(codes(issues)).toEqual(['warning model.alphabet.detected /molecules/0/residues']);
    expect(issues[0].message).toContain('"E" at 10');
    // 80% is plainly protein, and a nucleotide sequence is plainly one: neither warns.
    expect(codes(build({ key: 'p', residues: 'ACGTACGTEE' }).issues)).toEqual([]);
    expect(build({ key: 'n', residues: 'ACGTACGTAC' }).molecule?.alphabet).toBe('dna');
    // A declared alphabet is never second-guessed: the letter is an error, as before.
    expect(codes(build({ key: 'd', residues: 'AAAAAAAAAGE', alphabet: 'dna' }).issues)).toEqual([
      'error model.residues.invalid /molecules/0/residues',
    ]);
  });

  it('tells the letters of each alphabet from the 128-entry tables, and nothing outside ASCII', () => {
    const letters = {
      protein: 'ACDEFGHIKLMNPQRSTVWYBZJUOX*',
      dna: 'ACGTNRYSWKMBDHV',
      rna: 'ACGUNRYSWKMBDHV',
      generic: '',
    } as const;
    for (const [alphabet, set] of Object.entries(letters) as [keyof typeof letters, string][]) {
      for (let c = -1; c < 300; c++) {
        const generic = c > 32 && c < 127 && c !== 45 && c !== 46;
        expect(Molecules.isLetter(alphabet, c), `${alphabet} ${c}`).toBe(
          alphabet === 'generic' ? generic : c >= 0 && set.includes(String.fromCharCode(c)),
        );
      }
    }
  });
});

describe('residues', () => {
  it('stores them upper case and checks them against the alphabet', () => {
    const rna = build({ key: 'a', residues: 'acgu' }).molecule;
    expect(rna && [rna.alphabet, String.fromCharCode(...rna.code), rna.kind]).toEqual([
      'rna',
      'ACGU',
      MoleculeKind.polymer,
    ]);
    // A gap character is no residue, in any alphabet, declared or detected.
    expect(codes(build({ key: 'b', residues: 'AC-GT' }).issues)).toEqual([
      'error model.residues.invalid /molecules/0/residues',
    ]);
    expect(codes(build({ key: 'b', residues: 'AC.GT', alphabet: 'generic' }).issues)).toEqual([
      'error model.residues.invalid /molecules/0/residues',
    ]);
    expect(residuesOf({ key: 'c', residues: 'xyz1', alphabet: 'generic' })).toBe('XYZ1');
    expect(codes(build({ key: 'd', residues: 'ACGUT' }).issues)).toEqual([
      'error model.residues.invalid /molecules/0/residues',
    ]);
  });

  it('upper-cases a to z and nothing else (ASCII rules only)', () => {
    const every = Array.from({ length: 94 }, (_, i) => String.fromCharCode(33 + i))
      .filter((c) => c !== '-' && c !== '.')
      .join('');
    expect(residuesOf({ key: 'g', residues: every, alphabet: 'generic' })).toBe(
      every.replace(/[a-z]/g, (c) => c.toUpperCase()),
    );
  });

  it('rejects a character outside ASCII with its index, never truncating or folding it', () => {
    // 'Ł' (U+0141) was stored as 'A' and 'ß' upper-cased to "SS": `MKŁß` loaded as `MKASS`.
    const { molecule, issues } = build({ key: 'pdf', residues: 'MKŁß' });
    expect(molecule).toBeNull();
    expect(codes(issues)).toEqual(['error model.residues.invalid /molecules/0/residues']);
    expect(issues[0].message).toContain('U+141 at 2');
    // The index is a UTF-16 index; every code unit above 127 fails, wherever it sits and whatever it truncates to.
    for (let c = 128; c < 0x800; c++) {
      const r = build({ key: 'x', residues: `ACG${String.fromCharCode(c)}T`, alphabet: 'generic' });
      expect(codes(r.issues), `U+${c.toString(16)}`).toEqual(['error model.residues.invalid /molecules/0/residues']);
      expect(r.issues[0].message).toContain(' at 3 ');
    }
    expect(codes(build({ key: 'e', residues: 'A😀' }).issues)).toEqual([
      'error model.residues.invalid /molecules/0/residues',
    ]);
    // A surrogate pair is named by its code point, not its first half (U+D83D), at that half's index.
    expect(build({ key: 'e', residues: 'A😀' }).issues[0].message).toContain('U+1F600 at 1 ');
    expect(build({ key: 'e', residues: 'A\ud83dC' }).issues[0].message).toContain('U+D83D at 1 ');
  });

  it('rejects a polymer with no residues as its content', () => {
    expect(codes(build({ key: 'p', residues: '' }).issues)).toEqual([
      'error model.molecule.content /molecules/0/residues',
    ]);
    expect(codes(build({ key: 'p' }).issues)).toEqual(['error model.molecule.content /molecules/0/residues']);
  });
});

describe('ligands', () => {
  it('builds from heavy-atom elements, rejecting hydrogens, unknown symbols and residues', () => {
    const { molecule, issues } = build({
      key: 'hem',
      kind: 'ligand',
      elements: ['Fe', 'N', 'N', 'c', 'O'],
      chemCompIds: ['HEM'],
    });
    expect(issues).toEqual([]);
    expect([molecule?.kind, molecule?.alphabet, molecule?.length]).toEqual([MoleculeKind.ligand, null, 5]);
    expect(Array.from(molecule?.element ?? [], Molecules.elementSymbol)).toEqual(['Fe', 'N', 'N', 'C', 'O']);
    expect(String.fromCharCode(...(molecule?.code ?? []))).toBe('XXXXX');
    const bad: MoleculeInput[] = [
      { key: 'h', kind: 'ligand', elements: ['C', 'H'] },
      { key: 'q', kind: 'ligand', elements: ['Qq'] },
      { key: 'r', kind: 'ligand', residues: 'C' },
    ];
    expect(bad.flatMap((b) => codes(build(b).issues))).toEqual([
      'error model.ligand.element /molecules/0/elements/1',
      'error model.ligand.element /molecules/0/elements/0',
      'error model.molecule.content /molecules/0/elements',
    ]);
  });
});
