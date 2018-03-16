/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated 'CCD' schema file
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
// const Matrix = Schema.Matrix;
// const Vector = Schema.Vector;
const List = Schema.List;

export const CCD_Schema = {
    chem_comp: {
        formula: str,
        formula_weight: float,
        id: str,
        mon_nstd_parent_comp_id: List(',', x => x),
        name: str,
        one_letter_code: str,
        three_letter_code: str,
        type: Aliased<'D-peptide linking' | 'L-peptide linking' | 'D-peptide NH3 amino terminus' | 'L-peptide NH3 amino terminus' | 'D-peptide COOH carboxy terminus' | 'L-peptide COOH carboxy terminus' | 'DNA linking' | 'RNA linking' | 'L-RNA linking' | 'L-DNA linking' | 'DNA OH 5 prime terminus' | 'RNA OH 5 prime terminus' | 'DNA OH 3 prime terminus' | 'RNA OH 3 prime terminus' | 'D-saccharide 1,4 and 1,4 linking' | 'L-saccharide 1,4 and 1,4 linking' | 'D-saccharide 1,4 and 1,6 linking' | 'L-saccharide 1,4 and 1,6 linking' | 'L-saccharide' | 'D-saccharide' | 'saccharide' | 'non-polymer' | 'peptide linking' | 'peptide-like' | 'L-gamma-peptide, C-delta linking' | 'D-gamma-peptide, C-delta linking' | 'L-beta-peptide, C-gamma linking' | 'D-beta-peptide, C-gamma linking' | 'other'>(str),
        pdbx_synonyms: str,
        pdbx_type: str,
        pdbx_ambiguous_flag: str,
        pdbx_replaced_by: str,
        pdbx_replaces: str,
        pdbx_formal_charge: int,
        pdbx_model_coordinates_details: str,
        pdbx_model_coordinates_db_code: str,
        pdbx_ideal_coordinates_details: str,
        pdbx_ideal_coordinates_missing_flag: Aliased<'Y' | 'N'>(str),
        pdbx_model_coordinates_missing_flag: Aliased<'Y' | 'N'>(str),
        pdbx_initial_date: str,
        pdbx_modified_date: str,
        pdbx_processing_site: Aliased<'PDBE' | 'EBI' | 'PDBJ' | 'RCSB'>(str),
    },
    chem_comp_atom: {
        alt_atom_id: str,
        atom_id: str,
        charge: int,
        model_Cartn_x: coord,
        model_Cartn_y: coord,
        model_Cartn_z: coord,
        comp_id: str,
        type_symbol: str,
        pdbx_align: int,
        pdbx_ordinal: int,
        pdbx_model_Cartn_x_ideal: coord,
        pdbx_model_Cartn_y_ideal: coord,
        pdbx_model_Cartn_z_ideal: coord,
        pdbx_stereo_config: Aliased<'R' | 'S' | 'N'>(str),
        pdbx_aromatic_flag: Aliased<'Y' | 'N'>(str),
        pdbx_leaving_atom_flag: Aliased<'Y' | 'N'>(str),
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
    pdbx_chem_comp_descriptor: {
        comp_id: str,
        descriptor: str,
        type: Aliased<'SMILES_CANNONICAL' | 'SMILES_CANONICAL' | 'SMILES' | 'SMILES' | 'InChI' | 'InChI_MAIN' | 'InChI_MAIN_FORMULA' | 'InChI_MAIN_CONNECT' | 'InChI_MAIN_HATOM' | 'InChI_CHARGE' | 'InChI_STEREO' | 'InChI_ISOTOPE' | 'InChI_FIXEDH' | 'InChI_RECONNECT' | 'InChIKey'>(str),
        program: str,
        program_version: str,
    },
    pdbx_chem_comp_identifier: {
        comp_id: str,
        identifier: str,
        type: Aliased<'COMMON NAME' | 'SYSTEMATIC NAME' | 'CAS REGISTRY NUMBER' | 'PUBCHEM Identifier' | 'MDL Identifier' | 'SYNONYM'>(str),
        program: str,
        program_version: str,
    },
}

export type CCD_Schema = typeof CCD_Schema;
export interface CCD_Database extends Database<CCD_Schema> {}