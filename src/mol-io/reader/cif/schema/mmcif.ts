/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated 'mmCIF' schema file
 *
 * @author mol-star package (src/apps/schema-generator/generate)
 */

import { Database, Column } from 'mol-data/db'

import Schema = Column.Schema

const str = Schema.str;
const int = Schema.int;
const float = Schema.float;
const coord = Schema.coord;

const Aliased = Schema.Aliased;
const Matrix = Schema.Matrix;
const Vector = Schema.Vector;
const List = Schema.List;

export const mmCIF_Schema = {
    atom_site: {
        auth_asym_id: str,
        auth_atom_id: str,
        auth_comp_id: str,
        auth_seq_id: int,
        B_iso_or_equiv: float,
        Cartn_x: coord,
        Cartn_y: coord,
        Cartn_z: coord,
        group_PDB: Aliased<'ATOM' | 'HETATM'>(str),
        id: int,
        label_alt_id: str,
        label_asym_id: str,
        label_atom_id: str,
        label_comp_id: str,
        label_entity_id: str,
        label_seq_id: int,
        occupancy: float,
        type_symbol: str,
        pdbx_PDB_ins_code: str,
        pdbx_PDB_model_num: int,
        pdbx_formal_charge: int,
    },
    cell: {
        angle_alpha: float,
        angle_beta: float,
        angle_gamma: float,
        entry_id: str,
        length_a: float,
        length_b: float,
        length_c: float,
        Z_PDB: int,
        pdbx_unique_axis: str,
    },
    chem_comp: {
        formula: str,
        formula_weight: float,
        id: str,
        mon_nstd_flag: Aliased<'no' | 'n' | 'yes' | 'y'>(str),
        name: str,
        type: Aliased<'D-peptide linking' | 'L-peptide linking' | 'D-peptide NH3 amino terminus' | 'L-peptide NH3 amino terminus' | 'D-peptide COOH carboxy terminus' | 'L-peptide COOH carboxy terminus' | 'DNA linking' | 'RNA linking' | 'L-RNA linking' | 'L-DNA linking' | 'DNA OH 5 prime terminus' | 'RNA OH 5 prime terminus' | 'DNA OH 3 prime terminus' | 'RNA OH 3 prime terminus' | 'D-saccharide 1,4 and 1,4 linking' | 'L-saccharide 1,4 and 1,4 linking' | 'D-saccharide 1,4 and 1,6 linking' | 'L-saccharide 1,4 and 1,6 linking' | 'L-saccharide' | 'D-saccharide' | 'saccharide' | 'non-polymer' | 'peptide linking' | 'peptide-like' | 'L-gamma-peptide, C-delta linking' | 'D-gamma-peptide, C-delta linking' | 'L-beta-peptide, C-gamma linking' | 'D-beta-peptide, C-gamma linking' | 'other'>(str),
        pdbx_synonyms: str,
    },
    chem_comp_bond: {
        atom_id_1: str,
        atom_id_2: str,
        comp_id: str,
        value_order: Aliased<'sing' | 'doub' | 'trip' | 'quad' | 'arom' | 'poly' | 'delo' | 'pi'>(str),
        pdbx_ordinal: int,
        pdbx_stereo_config: Aliased<'E' | 'Z' | 'N'>(str),
        pdbx_aromatic_flag: Aliased<'Y' | 'N'>(str),
    },
    entity: {
        details: str,
        formula_weight: float,
        id: str,
        src_method: Aliased<'nat' | 'man' | 'syn'>(str),
        type: Aliased<'polymer' | 'non-polymer' | 'macrolide' | 'water'>(str),
        pdbx_description: str,
        pdbx_number_of_molecules: float,
        pdbx_mutation: str,
        pdbx_fragment: str,
        pdbx_ec: List(',', x => x),
    },
    entry: {
        id: str,
    },
    exptl: {
        entry_id: str,
        method: Aliased<'X-RAY DIFFRACTION' | 'NEUTRON DIFFRACTION' | 'FIBER DIFFRACTION' | 'ELECTRON CRYSTALLOGRAPHY' | 'ELECTRON MICROSCOPY' | 'SOLUTION NMR' | 'SOLID-STATE NMR' | 'SOLUTION SCATTERING' | 'POWDER DIFFRACTION' | 'INFRARED SPECTROSCOPY' | 'EPR' | 'FLUORESCENCE TRANSFER' | 'THEORETICAL MODEL'>(str),
    },
    struct: {
        entry_id: str,
        title: str,
    },
    struct_conf: {
        beg_label_asym_id: str,
        beg_label_comp_id: str,
        beg_label_seq_id: int,
        beg_auth_asym_id: str,
        beg_auth_comp_id: str,
        beg_auth_seq_id: int,
        conf_type_id: Aliased<'HELX_P' | 'HELX_OT_P' | 'HELX_RH_P' | 'HELX_RH_OT_P' | 'HELX_RH_AL_P' | 'HELX_RH_GA_P' | 'HELX_RH_OM_P' | 'HELX_RH_PI_P' | 'HELX_RH_27_P' | 'HELX_RH_3T_P' | 'HELX_RH_PP_P' | 'HELX_LH_P' | 'HELX_LH_OT_P' | 'HELX_LH_AL_P' | 'HELX_LH_GA_P' | 'HELX_LH_OM_P' | 'HELX_LH_PI_P' | 'HELX_LH_27_P' | 'HELX_LH_3T_P' | 'HELX_LH_PP_P' | 'HELX_N' | 'HELX_OT_N' | 'HELX_RH_N' | 'HELX_RH_OT_N' | 'HELX_RH_A_N' | 'HELX_RH_B_N' | 'HELX_RH_Z_N' | 'HELX_LH_N' | 'HELX_LH_OT_N' | 'HELX_LH_A_N' | 'HELX_LH_B_N' | 'HELX_LH_Z_N' | 'TURN_P' | 'TURN_OT_P' | 'TURN_TY1_P' | 'TURN_TY1P_P' | 'TURN_TY2_P' | 'TURN_TY2P_P' | 'TURN_TY3_P' | 'TURN_TY3P_P' | 'STRN'>(str),
        details: str,
        end_label_asym_id: str,
        end_label_comp_id: str,
        end_label_seq_id: int,
        end_auth_asym_id: str,
        end_auth_comp_id: str,
        end_auth_seq_id: int,
        id: str,
        pdbx_beg_PDB_ins_code: str,
        pdbx_end_PDB_ins_code: str,
        pdbx_PDB_helix_class: str,
        pdbx_PDB_helix_length: int,
        pdbx_PDB_helix_id: str,
    },
    struct_conn: {
        conn_type_id: Aliased<'covale' | 'disulf' | 'hydrog' | 'metalc' | 'mismat' | 'saltbr' | 'modres' | 'covale_base' | 'covale_sugar' | 'covale_phosphate'>(str),
        details: str,
        id: str,
        ptnr1_label_asym_id: str,
        ptnr1_label_atom_id: str,
        ptnr1_label_comp_id: str,
        ptnr1_label_seq_id: int,
        ptnr1_auth_asym_id: str,
        ptnr1_auth_comp_id: str,
        ptnr1_auth_seq_id: int,
        ptnr1_symmetry: str,
        ptnr2_label_asym_id: str,
        ptnr2_label_atom_id: str,
        ptnr2_label_seq_id: int,
        ptnr2_auth_asym_id: str,
        ptnr2_auth_comp_id: str,
        ptnr2_auth_seq_id: int,
        ptnr2_symmetry: str,
        pdbx_ptnr1_PDB_ins_code: str,
        pdbx_ptnr1_label_alt_id: str,
        pdbx_ptnr1_standard_comp_id: str,
        pdbx_ptnr2_PDB_ins_code: str,
        pdbx_ptnr2_label_alt_id: str,
        pdbx_ptnr3_PDB_ins_code: str,
        pdbx_ptnr3_label_alt_id: str,
        pdbx_ptnr3_label_asym_id: str,
        pdbx_ptnr3_label_atom_id: str,
        pdbx_ptnr3_label_comp_id: str,
        pdbx_ptnr3_label_seq_id: int,
        pdbx_PDB_id: str,
        pdbx_dist_value: float,
        pdbx_value_order: Aliased<'sing' | 'doub' | 'trip' | 'quad'>(str),
    },
    struct_conn_type: {
        criteria: str,
        id: Aliased<'covale' | 'disulf' | 'hydrog' | 'metalc' | 'mismat' | 'saltbr' | 'modres' | 'covale_base' | 'covale_sugar' | 'covale_phosphate'>(str),
        reference: str,
    },
    struct_keywords: {
        entry_id: str,
        text: List(',', x => x),
        pdbx_keywords: str,
    },
    struct_sheet_range: {
        beg_label_asym_id: str,
        beg_label_comp_id: str,
        beg_label_seq_id: int,
        end_label_asym_id: str,
        end_label_comp_id: str,
        end_label_seq_id: int,
        beg_auth_asym_id: str,
        beg_auth_comp_id: str,
        beg_auth_seq_id: int,
        end_auth_asym_id: str,
        end_auth_comp_id: str,
        end_auth_seq_id: int,
        id: str,
        sheet_id: str,
        pdbx_beg_PDB_ins_code: str,
        pdbx_end_PDB_ins_code: str,
    },
    symmetry: {
        entry_id: str,
        cell_setting: Aliased<'triclinic' | 'monoclinic' | 'orthorhombic' | 'tetragonal' | 'rhombohedral' | 'trigonal' | 'hexagonal' | 'cubic'>(str),
        Int_Tables_number: int,
        space_group_name_Hall: str,
        'space_group_name_H-M': str,
    },
    pdbx_struct_assembly: {
        method_details: str,
        oligomeric_details: str,
        oligomeric_count: int,
        details: str,
        id: str,
    },
    pdbx_struct_mod_residue: {
        id: int,
        auth_asym_id: str,
        auth_comp_id: str,
        auth_seq_id: int,
        PDB_ins_code: str,
        label_asym_id: str,
        label_comp_id: str,
        label_seq_id: int,
        parent_comp_id: str,
        details: str,
    },
    pdbx_struct_oper_list: {
        id: str,
        type: Aliased<'identity operation' | 'point symmetry operation' | 'helical symmetry operation' | 'crystal symmetry operation' | '3D crystal symmetry operation' | '2D crystal symmetry operation' | 'transform to point frame' | 'transform to helical frame' | 'transform to crystal frame' | 'transform to 2D crystal frame' | 'transform to 3D crystal frame' | 'build point asymmetric unit' | 'build helical asymmetric unit' | 'build 2D crystal asymmetric unit' | 'build 3D crystal asymmetric unit'>(str),
        name: str,
        symmetry_operation: str,
        matrix: Matrix(3, 3),
        vector: Vector(3),
    },
    pdbx_struct_assembly_gen: {
        asym_id_list: List(',', x => x),
        assembly_id: str,
        oper_expression: str,
    },
}

export type mmCIF_Schema = typeof mmCIF_Schema;
export type mmCIF_Database = Database<mmCIF_Schema>