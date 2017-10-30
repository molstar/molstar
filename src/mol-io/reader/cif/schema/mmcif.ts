/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Types, TypedFrame } from '../schema'

const str = Types.str;
const int = Types.int;
const float = Types.float;

const entry = {
    id: str
}

type EntityType = 'polymer' | 'non-polymer' | 'water'

const entity = {
    id: str,
    type: Types.aliased<EntityType>(str),
    src_method: str,
    pdbx_description: str,
    formula_weight: float,
    pdbx_number_of_molecules: int,
    details: str,
    pdbx_mutation: str,
    pdbx_fragment: str,
    pdbx_ec: str
}

const exptl = {
    entry_id: str,
    method: str
}

const cell = {
    entry_id: str,
    length_a: float,
    length_b: float,
    length_c: float,
    angle_alpha: float,
    angle_beta: float,
    angle_gamma: float,
    Z_PDB: int,
    pdbx_unique_axis: str
}

const symmetry = {
    entry_id: str,
    'space_group_name_H-M': str,
    'pdbx_full_space_group_name_H': str,
    cell_setting: str,
    Int_Tables_number: int,
    space_group_name_Hall: str
}

const struct_conf = {
    conf_type_id: str,
    id: str,
    pdbx_PDB_helix_id: int,
    beg_label_comp_id: str,
    beg_label_asym_id: str,
    beg_label_seq_id: int,
    pdbx_beg_PDB_ins_code: str,
    end_label_comp_id: str,
    end_label_asym_id: str,
    end_label_seq_id: int,
    pdbx_end_PDB_ins_code: str,
    beg_auth_comp_id: str,
    beg_auth_asym_id: str,
    beg_auth_seq_id: int,
    end_auth_comp_id: str,
    end_auth_asym_id: str,
    end_auth_seq_id: int,
    pdbx_PDB_helix_class: int,
    details: str,
    pdbx_PDB_helix_length: int
}

const struct_sheet_range = {
    sheet_id: str,
    id: int,
    beg_label_comp_id: str,
    beg_label_asym_id: str,
    beg_label_seq_id: int,
    pdbx_beg_PDB_ins_code: str,
    end_label_comp_id: str,
    end_label_asym_id: str,
    end_label_seq_id: int,
    pdbx_end_PDB_ins_code: str,
    beg_auth_comp_id: str,
    beg_auth_asym_id: str,
    beg_auth_seq_id: int,
    end_auth_comp_id: str,
    end_auth_asym_id: str,
    end_auth_seq_id: int
}

type StructConnTypeId =
    | 'covale'
    | 'covale_base'
    | 'covale_phosphate'
    | 'covale_sugar'
    | 'disulf'
    | 'hydrog'
    | 'metalc'
    | 'mismat'
    | 'modres'
    | 'saltbr'

type BondValueOrder =
    | 'SING'
    | 'DOUB'
    | 'TRIP'
    | 'QUAD'

const struct_conn = {
    id: str,
    conn_type_id: Types.aliased<StructConnTypeId>(str),
    pdbx_PDB_id: str,
    ptnr1_label_asym_id: str,
    ptnr1_label_comp_id: str,
    ptnr1_label_seq_id: int,
    ptnr1_label_atom_id: str,
    pdbx_ptnr1_label_alt_id: str,
    pdbx_ptnr1_PDB_ins_code: str,
    pdbx_ptnr1_standard_comp_id: str,
    ptnr1_symmetry: str,
    ptnr2_label_asym_id: str,
    ptnr2_label_comp_id: str,
    ptnr2_label_seq_id: int,
    ptnr2_label_atom_id: str,
    pdbx_ptnr2_label_alt_id: str,
    pdbx_ptnr2_PDB_ins_code: str,
    ptnr1_auth_asym_id: str,
    ptnr1_auth_comp_id: str,
    ptnr1_auth_seq_id: int,
    ptnr2_auth_asym_id: str,
    ptnr2_auth_comp_id: str,
    ptnr2_auth_seq_id: int,
    ptnr2_symmetry: str,
    pdbx_ptnr3_label_atom_id: str,
    pdbx_ptnr3_label_seq_id: int,
    pdbx_ptnr3_label_comp_id: str,
    pdbx_ptnr3_label_asym_id: str,
    pdbx_ptnr3_label_alt_id: str,
    pdbx_ptnr3_PDB_ins_code: str,
    details: str,
    pdbx_dist_value: float,
    pdbx_value_order: Types.aliased<BondValueOrder>(str)
}

const struct_conn_type = {
    id: Types.aliased<StructConnTypeId>(str),
    criteria: str,
    reference: str
}

const chem_comp_bond = {
    comp_id: str,
    pdbx_stereo_config: str,
    pdbx_ordinal: int,
    pdbx_aromatic_flag: Types.aliased<'Y' | 'N'>(str),
    atom_id_1: str,
    atom_id_2: str,
    value_order: Types.aliased<BondValueOrder>(str)
}

const pdbx_struct_assembly = {
    id: str,
    details: str,
    method_details: str,
    oligomeric_details: str,
    oligomeric_count: int
}

const pdbx_struct_assembly_gen = {
    assembly_id: str,
    oper_expression: str,
    asym_id_list: str
}

const pdbx_struct_oper_list = {
    id: str,
    type: str,
    name: str,
    symmetry_operation: str,
    matrix: Types.matrix(3, 3),
    vector: Types.vector(3)
}

const pdbx_struct_mod_residue = {
    id: int,
    label_asym_id: str,
    label_seq_id: int,
    label_comp_id: str,
    auth_asym_id: str,
    auth_seq_id: int,
    auth_comp_id: str,
    PDB_ins_code: str,
    parent_comp_id: str,
    details: str
}

const atom_site = {
    group_PDB: str,
    id: int,
    type_symbol: str,
    label_atom_id: str,
    label_alt_id: str,
    label_comp_id: str,
    label_asym_id: str,
    label_entity_id: str,
    label_seq_id: int,
    pdbx_PDB_ins_code: str,
    pdbx_formal_charge: str,
    Cartn_x: float,
    Cartn_y: float,
    Cartn_z: float,
    occupancy: float,
    B_iso_or_equiv: float,
    auth_atom_id: str,
    auth_comp_id: str,
    auth_asym_id: str,
    auth_seq_id: int,
    pdbx_PDB_model_num: int
}

export const Schema = {
    entry,
    entity,
    exptl,
    cell,
    symmetry,
    struct_conf,
    struct_sheet_range,
    struct_conn,
    struct_conn_type,
    chem_comp_bond,
    pdbx_struct_assembly,
    pdbx_struct_assembly_gen,
    pdbx_struct_oper_list,
    pdbx_struct_mod_residue,
    atom_site
};

export interface Frame extends TypedFrame<typeof Schema> { }