/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from '../../../mol-base/collections/column'
import { Shape as mmCIF } from '../../../mol-io/reader/cif/schema/mmcif'

export interface ElementSymbol extends String { '@type': 'element-symbol' }
export function ElementSymbol(s: string): ElementSymbol {
    // TODO: optimize?
    return s.toUpperCase() as any;
}

type Key = { key: Column<number> }

type _Atoms = Pick<mmCIF['atom_site'],
    | 'type_symbol'
    | 'label_atom_id'
    | 'auth_atom_id'
    | 'label_alt_id'
    | 'pdbx_formal_charge'
    | 'occupancy'
    | 'B_iso_or_equiv'>
    & Key
export interface Atoms extends _Atoms {
    source_row: Column<number>
}

type _Residues = Pick<mmCIF['atom_site'],
    | 'group_PDB'
    | 'label_comp_id'
    | 'auth_comp_id'
    | 'label_seq_id'
    | 'auth_seq_id'
    | 'pdbx_PDB_ins_code'>
    & Key
export interface Residues extends _Residues { }

type _Chains = Pick<mmCIF['atom_site'],
    | 'label_asym_id'
    | 'auth_asym_id'
    | 'auth_comp_id'
    | 'label_entity_id'
    | 'pdbx_PDB_model_num'>
    & Key
export interface Chains extends _Chains {
    enityDataIndex: Column<number>
}

type _EntityData = mmCIF['entity']
export interface EntityData extends _EntityData { }

export interface Macromolecule {
    atoms: Atoms,
    residues: Residues,
    chains: Chains,
    entityData: EntityData
}

export default Macromolecule