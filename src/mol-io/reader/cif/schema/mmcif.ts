/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Field, TypedFrame } from '../schema'

const str = Field.str();
const int = Field.int();
const float = Field.float();

const entry = {
    id: str
}

const entity = {
    id: str,
    type: str as Field.Schema<'polymer' | 'non-polymer' | 'water'>,
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
    space_group_name_HM: Field.str({ alias: 'space_group_name_H-M' }),
    pdbx_full_space_group_name_HM: Field.str({ alias: 'pdbx_full_space_group_name_H-M' }),
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
    conn_type_id: str as Field.Schema<StructConnTypeId>,
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
    pdbx_value_order: str as Field.Schema<BondValueOrder>
}

const struct_conn_type = {
    id: str as Field.Schema<StructConnTypeId>,
    criteria: str,
    reference: str
}

const chem_comp_bond = {
    comp_id: str,
    pdbx_stereo_config: str,
    pdbx_ordinal: int,
    pdbx_aromatic_flag: str as Field.Schema<'Y' | 'N'>,
    atom_id_1: str,
    atom_id_2: str,
    value_order: str as Field.Schema<BondValueOrder>
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
    matrix: Field.matrix(3, 3),
    vector: Field.vector(3)
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

const mmCIF = {
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
type mmCIF = TypedFrame<typeof mmCIF>
export default mmCIF;