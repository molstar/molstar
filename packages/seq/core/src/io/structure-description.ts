/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The host-neutral description of a structure (_spec/README.md "Seq io"): what every adapter of a structure fills (a
// file's reader, a host holding its model), named in mmCIF's words, held as columns. Rows reference rows of other
// tables by index; `structureNodes` turns it into the nodes molecules are read as.

/** An entity's type (`entity.type`). */
export type EntityType = 'polymer' | 'non-polymer' | 'branched' | 'water';

/** `entity` (and `entity_poly`): a row an entity. */
export interface EntityTable {
  /** `entity.id`: what an entity's node is keyed and named by. */
  readonly id: readonly string[];
  readonly type: readonly EntityType[];
  /** `entity_poly.type` (`polypeptide(L)`, `polyribonucleotide`, `polydeoxyribonucleotide`); a polymer's only. */
  readonly polymerType?: readonly (string | undefined)[];
  /** `entity.pdbx_description`. */
  readonly name?: readonly (string | undefined)[];
  /** `entity_poly.pdbx_seq_one_letter_code_can`: a polymer's canonical letters; absent, from its components. */
  readonly letters?: readonly (string | undefined)[];
}

/**
 * The components a position of an entity (`entity_poly_seq`, `pdbx_entity_branch_list`, `pdbx_entity_nonpoly`): a row
 * a position, several at one position its alternatives (`hetero`), the first giving the letter.
 */
export interface ComponentTable {
  /** The entity's row. */
  readonly entity: Int32Array;
  /** The position in the entity, from 0. */
  readonly position: Int32Array;
  /** The component's code (`MSE`, `NAG`). */
  readonly component: readonly string[];
}

/** An entity's identifiers (a CCD code, a PRD id, WURCS, LINUCS): a row each. */
export interface IdentifierTable {
  readonly entity: Int32Array;
  /** An open lowercase word (`ccd`, `prd`, `wurcs`). */
  readonly namespace: readonly string[];
  readonly value: readonly string[];
}

/** The chain instances of the assembly described (`struct_asym` under its operators): a row a unit. */
export interface ChainTable {
  readonly label_asym_id: readonly string[];
  readonly auth_asym_id?: readonly (string | undefined)[];
  /** The entity's row. */
  readonly entity: Int32Array;
  /** The instance's operator (`1_555`, `ASM-X0-1`); absent: the deposited model's chain. */
  readonly operator?: readonly (string | undefined)[];
}

/** The residues a chain instance has (`pdbx_poly_seq_scheme`, `pdbx_branch_scheme`, `pdbx_nonpoly_scheme`). */
export interface ResidueTable {
  /** The chain's row. */
  readonly chain: Int32Array;
  /** The position in its entity, from 0. */
  readonly position: Int32Array;
  readonly auth_seq_id?: Int32Array;
  /** `pdbx_PDB_ins_code` (`''` or absent: none). */
  readonly ins_code?: readonly (string | undefined)[];
  /** 1 where the residue was modeled; absent: every listed residue was, an unlisted one was not. */
  readonly observed?: Uint8Array;
}

/** Links between residues of chain instances (`struct_conn`, a branch's links): a row a bond. */
export interface LinkTable {
  readonly a_chain: Int32Array;
  readonly a_position: Int32Array;
  readonly b_chain: Int32Array;
  readonly b_position: Int32Array;
  /** `covalent`, `disulfide`, `glycosidic`, `metal`, any word. */
  readonly type?: readonly (string | undefined)[];
  readonly a_atom?: readonly (string | null)[];
  readonly b_atom?: readonly (string | null)[];
}

/**
 * Links between residues of one molecule (a branch's `pdbx_entity_branch_link`): a row a bond, read as the molecule's
 * `links` child (molecules.md, Links). A polymer's are not read: a sequence has no links child.
 */
export interface EntityLinkTable {
  /** The entity's row. */
  readonly entity: Int32Array;
  /** Positions in the entity, from 0. */
  readonly a_position: Int32Array;
  readonly b_position: Int32Array;
  /** `glycosidic`, any word. */
  readonly type?: readonly (string | undefined)[];
  readonly a_atom?: readonly (string | null)[];
  readonly b_atom?: readonly (string | null)[];
  /** The bond's order (`sing`, `doub`). */
  readonly order?: readonly (string | null)[];
}

/** The components described (`chem_comp`): a row each. */
export interface CompoundTable {
  readonly id: readonly string[];
  readonly name?: readonly (string | undefined)[];
  /** `chem_comp.type` (`L-PEPTIDE LINKING`, `NON-POLYMER`). */
  readonly type?: readonly (string | undefined)[];
  /** `chem_comp.formula` (`C3 H2 O4 -2`, `Zn 2`): what tells an ion (one atom) from a ligand. */
  readonly formula?: readonly (string | undefined)[];
  /** `chem_comp.mon_nstd_parent_comp_id`, split. */
  readonly parents?: readonly (readonly string[] | undefined)[];
}

/** A molecule's heavy atoms (never hydrogens): a row an atom of a non-polymer or branched entity. */
export interface AtomTable {
  readonly entity: Int32Array;
  /** The position of the atom's residue in its entity, non-decreasing within one. */
  readonly residue: Int32Array;
  readonly name: readonly string[];
  readonly element: readonly string[];
}

/** A structure as its adapters describe it. */
export interface StructureDescription {
  /** The entry's name: its complex's `name`. */
  readonly name?: string;
  readonly entities: EntityTable;
  readonly components: ComponentTable;
  readonly identifiers?: IdentifierTable;
  readonly chains: ChainTable;
  readonly residues: ResidueTable;
  readonly links?: LinkTable;
  /** Links within a molecule (a glycan's own bonds). */
  readonly entityLinks?: EntityLinkTable;
  readonly compounds?: CompoundTable;
  readonly atoms?: AtomTable;
  /** Named values a residue row (`b_factor`, `occupancy`, `plddt`; NaN: none): the channels a `values` node reads. */
  readonly values?: Readonly<Record<string, Float32Array>>;
}
