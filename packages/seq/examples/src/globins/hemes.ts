/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// 4HHB's heme contacts (derived at recording, `scripts/record.ts globins` checks them against the entry, the fixtures
// test against `fixtures/structures/4HHB.cif.gz`): per pair of units, a heme's (by its chain, `auth_asym_id`) and a
// chain's, the chain's residues with a heavy atom within 4 Å of the heme's, 1-based `label_seq_id` (equal to `auth`
// in 4HHB; UniProt's number is one more, its initiator Met counted). Each heme touches its own chain only; C and D
// are A's and B's copies, their contacts not quite the same (C lacks V62, D lacks L31 and L88 and has L91).

/** The four chains: each a molecule's copy (α: A and C, β: B and D), its heme's asym id and author number. */
export const HemoglobinChains = [
  { chain: 'A', of: 'hba', heme: 'E', hemeAuth: '142' },
  { chain: 'B', of: 'hbb', heme: 'G', hemeAuth: '148' },
  { chain: 'C', of: 'hba', heme: 'H', hemeAuth: '142' },
  { chain: 'D', of: 'hbb', heme: 'J', hemeAuth: '148' },
] as const;

/**
 * The heme's heavy atoms (EC-64) as 4HHB's first heme has them (`atom_site` of asym E, in order; every heme of the
 * entry the same 43): names, then elements (`scripts/record.ts globins` checks both against the entry).
 */
export const HemeAtoms: { readonly name: readonly string[]; readonly element: readonly string[] } = heavyAtoms([
  ...['CHA', 'CHB', 'CHC', 'CHD'],
  ...['C1A', 'C2A', 'C3A', 'C4A', 'CMA', 'CAA', 'CBA', 'CGA', 'O1A', 'O2A'],
  ...['C1B', 'C2B', 'C3B', 'C4B', 'CMB', 'CAB', 'CBB'],
  ...['C1C', 'C2C', 'C3C', 'C4C', 'CMC', 'CAC', 'CBC'],
  ...['C1D', 'C2D', 'C3D', 'C4D', 'CMD', 'CAD', 'CBD', 'CGD', 'O1D', 'O2D'],
  ...['NA', 'NB', 'NC', 'ND', 'FE'],
]);

/** Names and their elements: a name's first letter, iron's `FE`. */
function heavyAtoms(name: string[]) {
  return { name, element: name.map((n) => (n === 'FE' ? 'FE' : n[0])) };
}

/** Each heme's contacts per pair of units: the heme's chain, the chain it touches, the residues touched there. */
export const HemeUnitContacts: readonly { heme: string; chain: string; residues: readonly number[] }[] = [
  { heme: 'A', chain: 'A', residues: [42, 43, 45, 46, 58, 61, 62, 83, 86, 87, 91, 93, 97, 98, 101, 136] },
  { heme: 'B', chain: 'B', residues: [31, 41, 42, 63, 66, 67, 70, 88, 92, 96, 98, 102, 103, 106, 141] },
  { heme: 'C', chain: 'C', residues: [42, 43, 45, 46, 58, 61, 83, 86, 87, 91, 93, 97, 98, 101, 136] },
  { heme: 'D', chain: 'D', residues: [41, 42, 63, 66, 67, 70, 91, 92, 96, 98, 102, 103, 106, 141] },
];

/** Each molecule's heme pocket (its first chain's contacts: A's for α, B's for β). */
export const HemeContacts = {
  hba: HemeUnitContacts[0].residues,
  hbb: HemeUnitContacts[1].residues,
} as const;

/** The distal and proximal histidines (E7 and F8), 4HHB's numbers. */
export const HemeHistidines = { hba: { distal: 58, proximal: 87 }, hbb: { distal: 63, proximal: 92 } } as const;
