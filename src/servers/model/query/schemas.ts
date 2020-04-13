/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from '../../../mol-io/writer/cif';

const InteractionCategories = new Set([
    'entry',
    'entity',
    'exptl',
    'cell',
    'symmetry',
    'struct_conf',
    'struct_sheet_range',
    'entity_poly',
    'struct_asym',
    'struct_conn',
    'struct_conn_type',
    'pdbx_struct_mod_residue',
    'chem_comp_bond',
    'atom_sites',
    'atom_site',
    'pdbx_entity_branch',
    'pdbx_entity_branch_link',
    'pdbx_branch_scheme'
]);

const AssemblyCategories = new Set([
    'entry',
    'entity',
    'exptl',
    'cell',
    'symmetry',
    'struct_conf',
    'struct_sheet_range',
    'entity_poly',
    'entity_poly_seq',
    'pdbx_nonpoly_scheme',
    'struct_asym',
    'struct_conn',
    'struct_conn_type',
    'pdbx_struct_mod_residue',
    'chem_comp_bond',
    'atom_sites',
    'atom_site',
    'pdbx_entity_branch',
    'pdbx_entity_branch_link',
    'pdbx_branch_scheme'
]);

export const QuerySchemas = {
    interaction: <CifWriter.Category.Filter>{
        includeCategory(name) { return InteractionCategories.has(name); },
        includeField(cat, field) { return true; }
    },
    assembly: <CifWriter.Category.Filter>{
        includeCategory(name) { return AssemblyCategories.has(name); },
        includeField(cat, field) { return true; }
    }
};