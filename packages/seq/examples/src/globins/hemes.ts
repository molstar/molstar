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

/** The four chains: each a molecule's copy (α: A and C, β: B and D), its heme's asym id. */
export const HemoglobinChains = [
  { chain: 'A', of: 'hba', heme: 'E' },
  { chain: 'B', of: 'hbb', heme: 'G' },
  { chain: 'C', of: 'hba', heme: 'H' },
  { chain: 'D', of: 'hbb', heme: 'J' },
] as const;

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
