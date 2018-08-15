/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

 import { mmCIF_Schema } from './mmcif';

export const mmCIF_residueId_schema = {
    label_comp_id: mmCIF_Schema.atom_site.label_comp_id,
    label_seq_id: mmCIF_Schema.atom_site.label_seq_id,
    pdbx_PDB_ins_code: mmCIF_Schema.atom_site.pdbx_PDB_ins_code,
    label_asym_id: mmCIF_Schema.atom_site.label_asym_id,
    label_entity_id: mmCIF_Schema.atom_site.label_entity_id,
    auth_comp_id: mmCIF_Schema.atom_site.auth_atom_id,
    auth_seq_id: mmCIF_Schema.atom_site.auth_seq_id,
    auth_asym_id: mmCIF_Schema.atom_site.auth_asym_id
}