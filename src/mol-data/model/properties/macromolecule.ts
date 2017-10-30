/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from '../../../mol-base/collections/column'
import { Shape as mmCIF } from '../../../mol-io/reader/cif/schema/mmcif'

export type Table<Data> = { [E in keyof Data]: Column<Data[E]> }

export interface ElementSymbol extends String { '@type': 'element-symbol' }
export function ElementSymbol(s: string): ElementSymbol {
    // TODO: optimize?
    return s.toUpperCase() as any;
}

export interface Atoms extends Table<{
    // unique number for each atom
    key: number,

    id: number,
    type_symbol: ElementSymbol,
    label_atom_id: string,
    auth_atom_id: string,
    label_alt_id: string,
    auth_alt_id: string,
    pdbx_formal_charge: string,
    occupancy: number,
    B_iso_or_equiv: number
}> { }

export interface Residues extends Table<{
    // unique number for each residue
    key: number,

    group_PDB: string,

    label_comp_id: string,
    auth_comp_id: string,
    label_seq_id: number,
    auth_seq_id: number,
    pdbx_PDB_ins_code: string
}> { }

export interface Chains extends Table<{
    // unique number for each chain
    key: number,

    label_asym_id: string,
    auth_asym_id: string
}> { }

export interface Entities extends Table<{
    // unique number for each entity
    // row index to the EntityData table
    key: number,
    label_entity_id: string,
    pdbx_PDB_model_num: number
}> { }

type _EntityData = mmCIF['entity']
export interface EntityData extends _EntityData { }

export interface Macromolecule {
    atoms: Atoms,
    residues: Residues,
    chains: Chains,
    entities: Entities
}

export default Macromolecule