/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

/**
 * The canonical residues, whose bond orders the built-in bonds always give (the core order table, the backbone C=O
 * and phosphate P=O rules): `auto` leaves them alone (_spec §1.4). Not Mol*'s `AminoAcidNames`, which holds modified
 * residues (PTR, SEP, PCA, the D-forms…) whose side-chain orders the core table does not give.
 */
export const KnownComponents: ReadonlySet<string> = new Set([
  // the 20 standard amino acids
  'ALA',
  'ARG',
  'ASN',
  'ASP',
  'CYS',
  'GLN',
  'GLU',
  'GLY',
  'HIS',
  'ILE',
  'LEU',
  'LYS',
  'MET',
  'PHE',
  'PRO',
  'SER',
  'THR',
  'TRP',
  'TYR',
  'VAL',
  // the standard nucleotides
  'A',
  'C',
  'G',
  'U',
  'DA',
  'DC',
  'DG',
  'DT',
]);
