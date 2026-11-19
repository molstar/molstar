/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The shape of what the record script derives from the Complex's structures (`derived-data.ts`): 1FIN's chains and
// copies, the contacts across each copy, and the ligands with their contacting residues. The chains' letters are the
// files' (`1FIN.seq`, the UniProt entries), read where the views read them.

export interface Derived {
  /** The author number of the cyclin A construct's first residue (P20248's numbering). */
  readonly cyclinFirst: number;
  /** 1FIN's chains: label asym id, the entity's system, the copy of the complex. */
  readonly chains: readonly { readonly id: string; readonly system: 'cdk2' | 'cyclin'; readonly copy: number }[];
  /** Each contact across a copy: the copy, CDK2's residue, cyclin A's (0-based, on the construct). */
  readonly contacts: readonly (readonly [copy: number, cdk2: number, cyclin: number])[];
  /** Each ligand: 1FIN's ATPs first, then one of each other entry; its contacts CDK2's residue and its heavy atom. */
  readonly ligands: readonly {
    readonly comp: string;
    readonly entry: string;
    readonly chain: string;
    readonly resolution: number;
    readonly atoms: readonly string[];
    readonly contacts: readonly (readonly [residue: number, atom: number])[];
  }[];
}
