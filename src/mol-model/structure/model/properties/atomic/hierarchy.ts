/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from 'mol-data/db'
import { Segmentation } from 'mol-data/int'
import { mmCIF_Schema as mmCIF } from 'mol-io/reader/cif/schema/mmcif'
import { ElementSymbol } from '../../types'
import { Element } from '../../../structure'

export const AtomsSchema = {
    type_symbol: Column.Schema.Aliased<ElementSymbol>(mmCIF.atom_site.type_symbol),
    label_atom_id: mmCIF.atom_site.label_atom_id,
    auth_atom_id: mmCIF.atom_site.auth_atom_id,
    label_alt_id: mmCIF.atom_site.label_alt_id,
    pdbx_formal_charge: mmCIF.atom_site.pdbx_formal_charge
    // id, occupancy and B_iso_or_equiv are part of conformation
};

export type AtomsSchema = typeof AtomsSchema
export interface Atoms extends Table<AtomsSchema> { }

export const ResiduesSchema = {
    group_PDB: mmCIF.atom_site.group_PDB,
    label_comp_id: mmCIF.atom_site.label_comp_id,
    auth_comp_id: mmCIF.atom_site.auth_comp_id,
    label_seq_id: mmCIF.atom_site.label_seq_id,
    auth_seq_id: mmCIF.atom_site.auth_seq_id,
    pdbx_PDB_ins_code: mmCIF.atom_site.pdbx_PDB_ins_code
};
export type ResiduesSchema = typeof ResiduesSchema
export interface Residues extends Table<ResiduesSchema> { }

export const ChainsSchema = {
    label_asym_id: mmCIF.atom_site.label_asym_id,
    auth_asym_id: mmCIF.atom_site.auth_asym_id,
    label_entity_id: mmCIF.atom_site.label_entity_id
}
export type ChainsSchema = typeof ChainsSchema
export interface Chains extends Table<ChainsSchema> { }

export interface AtomicData {
    atoms: Atoms,
    residues: Residues,
    chains: Chains
}

export interface AtomicSegments {
    residueSegments: Segmentation<Element>,
    chainSegments: Segmentation<Element>
    // TODO: include entity segments?
}

export interface AtomicKeys {
    // TODO: since Atoms must be sorted now, get rid of keys
    // TODO: include (lazily computed) "entity/chain/residue" indices?

    // assign a key to each residue index.
    residueKey: ArrayLike<number>,
    // assign a key to each chain index
    chainKey: ArrayLike<number>,
    // assigne a key to each chain index
    // also index to the Entities table.
    entityKey: ArrayLike<number>,

    findChainKey(entityId: string, label_asym_id: string): number,

    /** Unique number for each of the residue. Also the index of the 1st occurence of this residue. */
    findResidueKey(entityId: string, label_asym_id: string, label_comp_id: string, auth_seq_id: number, pdbx_PDB_ins_code: string): number
}

type _Hierarchy = AtomicData & AtomicSegments & AtomicKeys
export interface AtomicHierarchy extends _Hierarchy { }