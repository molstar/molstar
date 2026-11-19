/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// 4HHB's heme contacts (derived at recording, `scripts/record.ts globins` checks them against the entry): the residues
// of each heme's own chain with a heavy atom within 4 Å of the heme's, 1-based `label_seq_id` (equal to `auth` in
// 4HHB; UniProt's number is one more, its initiator Met counted). Chains A and B: C and D are their copies (C lacks
// V62 and D L31 and L88, has L91; the view shows the molecules' contacts as A's and B's).

export const HemeContacts = {
  hba: [42, 43, 45, 46, 58, 61, 62, 83, 86, 87, 91, 93, 97, 98, 101, 136],
  hbb: [31, 41, 42, 63, 66, 67, 70, 88, 92, 96, 98, 102, 103, 106, 141],
} as const;

/** The distal and proximal histidines (E7 and F8), 4HHB's numbers. */
export const HemeHistidines = { hba: { distal: 58, proximal: 87 }, hbb: { distal: 63, proximal: 92 } } as const;
