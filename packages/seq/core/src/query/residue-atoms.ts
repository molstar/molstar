/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The heavy atoms of the standard components (_spec/README.md "Queries", Properties; the author, 2026-11-20, on a
// polymer's atoms: "Residues, as a structure viewer reads them"): where a residue's atoms are not held, an atom field
// selects the residues whose standard component has that atom (CA on every amino acid, P on nucleotides), and a
// structure viewer given the residues picks the atoms. The names are the CCD's (PDB's) for the twenty amino acids,
// SEC and the nucleotides; the terminal OXT and the 5' phosphate's absence are not modeled (a residue has what its
// chain's interior one has).

const Backbone = ['N', 'CA', 'C', 'O'];
const Sugar = ['P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"];
const Purine = ['N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4'];
const Pyrimidine = ['N1', 'C2', 'O2', 'N3', 'C4', 'C5', 'C6'];

const SideChains: Readonly<Record<string, string>> = {
  ALA: 'CB',
  ARG: 'CB CG CD NE CZ NH1 NH2',
  ASN: 'CB CG OD1 ND2',
  ASP: 'CB CG OD1 OD2',
  CYS: 'CB SG',
  GLN: 'CB CG CD OE1 NE2',
  GLU: 'CB CG CD OE1 OE2',
  GLY: '',
  HIS: 'CB CG ND1 CD2 CE1 NE2',
  ILE: 'CB CG1 CG2 CD1',
  LEU: 'CB CG CD1 CD2',
  LYS: 'CB CG CD CE NZ',
  MET: 'CB CG SD CE',
  PHE: 'CB CG CD1 CD2 CE1 CE2 CZ',
  PRO: 'CB CG CD',
  SER: 'CB OG',
  THR: 'CB OG1 CG2',
  TRP: 'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
  TYR: 'CB CG CD1 CD2 CE1 CE2 CZ OH',
  VAL: 'CB CG1 CG2',
  SEC: 'CB SE',
  PYL: 'CB CG CD CE NZ C2 O2 CA2 CB2 CG2 CD2 CE2 N2',
  UNK: '',
};

const Bases: Readonly<Record<string, readonly string[]>> = {
  A: [...Purine, 'N6'],
  G: [...Purine, 'O6', 'N2'],
  C: [...Pyrimidine, 'N4'],
  U: [...Pyrimidine, 'O4'],
  T: [...Pyrimidine, 'O4', 'C7'],
  N: [],
};

const table = new Map<string, readonly string[]>();
for (const [code, side] of Object.entries(SideChains))
  table.set(code, Object.freeze([...Backbone, ...(side ? side.split(' ') : [])]));
for (const [base, atoms] of Object.entries(Bases)) {
  table.set(base, Object.freeze([...Sugar, "O2'", ...atoms]));
  if (base !== 'U') table.set(`D${base}`, Object.freeze([...Sugar, ...atoms]));
}

/** A standard component's heavy atom names, or undefined for a code the table lacks. */
export function standardAtoms(code: string): readonly string[] | undefined {
  return table.get(code);
}

const elementSets = new Map<string, readonly string[] | undefined>();

/** A standard component's heavy atoms' element symbols, each once; made once a code. */
export function standardElements(code: string): readonly string[] | undefined {
  if (elementSets.has(code)) return elementSets.get(code);
  const names = table.get(code);
  const out = names && Object.freeze([...new Set(names.map(atomElement))]);
  elementSets.set(code, out);
  return out;
}

/** An atom name's element symbol (`SE` is selenium, `SG` sulfur: the CCD's names of the standard components). */
export function atomElement(name: string): string {
  return name === 'SE' ? 'SE' : name[0];
}
