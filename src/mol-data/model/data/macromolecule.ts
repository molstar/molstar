/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from '../../../mol-base/collections/column'
import Table from '../../../mol-base/collections/table'
import { Schema as mmCIF } from '../../../mol-io/reader/cif/schema/mmcif'

const _esCache = Object.create(null);
export interface ElementSymbol extends String { '@type': 'element-symbol' }
export function ElementSymbol(s: string): ElementSymbol {
    return _esCache[s] || (_esCache[s] = s.toUpperCase());
}

export const AtomsSchema = {
    type_symbol: Column.Type.aliased<ElementSymbol>(mmCIF.atom_site.type_symbol),
    label_atom_id: mmCIF.atom_site.label_atom_id,
    auth_atom_id: mmCIF.atom_site.auth_atom_id,
    label_alt_id: mmCIF.atom_site.label_alt_id,
    pdbx_formal_charge: mmCIF.atom_site.pdbx_formal_charge,
    occupancy: mmCIF.atom_site.occupancy,
    B_iso_or_equiv: mmCIF.atom_site.B_iso_or_equiv
};

export interface Atoms extends Table<typeof AtomsSchema> { }

export const ResiduesSchema = {
    group_PDB: mmCIF.atom_site.group_PDB,
    label_comp_id: mmCIF.atom_site.label_comp_id,
    auth_comp_id: mmCIF.atom_site.auth_comp_id,
    label_seq_id: mmCIF.atom_site.label_seq_id,
    auth_seq_id: mmCIF.atom_site.auth_seq_id,
    pdbx_PDB_ins_code: mmCIF.atom_site.pdbx_PDB_ins_code,

    key: Column.Type.int
};

export interface Residues extends Table<typeof AtomsSchema> { }

export const ChainsSchema = {
    label_asym_id: mmCIF.atom_site.label_asym_id,
    auth_asym_id: mmCIF.atom_site.auth_asym_id,
    auth_comp_id: mmCIF.atom_site.auth_comp_id,
    label_entity_id: mmCIF.atom_site.label_entity_id,
    pdbx_PDB_model_num: mmCIF.atom_site.pdbx_PDB_model_num,

    key: Column.Type.int,
    entityIndex: Column.Type.int
}

export interface Chains extends Table<typeof ChainsSchema> { }

export const EntitySchema = mmCIF['entity']
export interface Entities extends Table<typeof EntitySchema> { }

export interface Macromolecule {
    atoms: Atoms,
    residues: Residues,
    chains: Chains,
    entities: Entities
}

export default Macromolecule