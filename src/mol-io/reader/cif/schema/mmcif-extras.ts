/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { mmCIF_Schema } from './mmcif';
import { Column } from '../../../../mol-data/db';

export const mmCIF_residueId_schema = {
    label_comp_id: mmCIF_Schema.atom_site.label_comp_id,
    label_seq_id: mmCIF_Schema.atom_site.label_seq_id,
    pdbx_PDB_ins_code: mmCIF_Schema.atom_site.pdbx_PDB_ins_code,
    label_asym_id: mmCIF_Schema.atom_site.label_asym_id,
    label_entity_id: mmCIF_Schema.atom_site.label_entity_id,
    auth_comp_id: mmCIF_Schema.atom_site.auth_atom_id,
    auth_seq_id: mmCIF_Schema.atom_site.auth_seq_id,
    auth_asym_id: mmCIF_Schema.atom_site.auth_asym_id
};

export const mmCIF_chemCompBond_schema = {
    ...mmCIF_Schema.chem_comp_bond,
    /** Indicates if the bond entry was taken from the protonation variant dictionary */
    molstar_protonation_variant: Column.Schema.Str()
};