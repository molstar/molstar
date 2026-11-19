/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The shape of what the record script derives from the Complex's structures (`derived-data.ts`): 1FIN's chains and
// copies, the contacts across each copy, and the ligands with their contacting residues.

export interface Derived {
  /**
   * The chains' letters: 1FIN's CDK2 entity and cyclin A construct (as `1FIN.seq` has them), UniProt's P24941 and P20248
   * (as the entries have them): held for the blocks, whose lengths a layout needs before a source opens.
   */
  readonly letters: {
    readonly cdk2: string;
    readonly construct: string;
    readonly P24941: string;
    readonly P20248: string;
  };
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
