/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Yakov Pechersky <ffxen158@gmail.com>
 */

import { CifCategory, CifField } from '../../../mol-io/reader/cif';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { Tokens } from '../../../mol-io/reader/common/text/tokenizer';

export function parseConect(lines: Tokens, lineStart: number, lineEnd: number, sites: { [K in keyof mmCIF_Schema['atom_site']]?: CifField }): CifCategory {
    const idMap: { [k: string]: number } = {};
    for (let i = 0, il = sites.id!.rowCount; i < il; ++i) {
        idMap[sites.id!.str(i)] = i;
    }

    const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1]);

    const id: string[] = [];
    const conn_type_id: string[] = [];

    const ptnr1_label_asym_id: string[] = [];
    const ptnr1_label_seq_id: number[] = [];
    const ptnr1_auth_seq_id: number[] = [];
    const ptnr1_label_atom_id: string[] = [];
    const ptnr1_label_alt_id: string[] = [];
    const ptnr1_PDB_ins_code: string[] = [];

    const ptnr2_label_asym_id: string[] = [];
    const ptnr2_label_seq_id: number[] = [];
    const ptnr2_auth_seq_id: number[] = [];
    const ptnr2_label_atom_id: string[] = [];
    const ptnr2_label_alt_id: string[] = [];
    const ptnr2_PDB_ins_code: string[] = [];

    const pos = [11, 16, 21, 26];

    let k = 1;

    for (let i = lineStart; i < lineEnd; i++) {
        const line = getLine(i);
        const idxA = idMap[parseInt(line.substr(6, 5))];

        const bondIndex: {[k: number]: number} = {};

        if (idxA === undefined) continue;

        for (let j = 0; j < 4; ++j) {
            const idB = parseInt(line.substr(pos[j], 5));
            if (Number.isNaN(idB)) continue;

            const idxB = idMap[idB];
            if (idxB === undefined) continue;
            if (idxA > idxB) continue;

            // TODO: interpret records where a 'idxB' atom is given multiple times
            // as double/triple bonds, e.g. CONECT 1529 1528 1528 is a double bond
            if (bondIndex[idxB] !== undefined) continue;

            id.push(`covale${k}`);
            conn_type_id.push('covale');

            ptnr1_label_asym_id.push(sites.label_asym_id!.str(idxA));
            ptnr1_label_seq_id.push(sites.label_seq_id!.int(idxA));
            ptnr1_auth_seq_id.push(sites.auth_seq_id!.int(idxA));
            ptnr1_label_atom_id.push(sites.label_atom_id!.str(idxA));
            ptnr1_label_alt_id.push(sites.label_alt_id!.str(idxA));
            ptnr1_PDB_ins_code.push(sites.pdbx_PDB_ins_code!.str(idxA));

            ptnr2_label_asym_id.push(sites.label_asym_id!.str(idxB));
            ptnr2_label_seq_id.push(sites.label_seq_id!.int(idxB));
            ptnr2_auth_seq_id.push(sites.auth_seq_id!.int(idxB));
            ptnr2_label_atom_id.push(sites.label_atom_id!.str(idxB));
            ptnr2_label_alt_id.push(sites.label_alt_id!.str(idxB));
            ptnr2_PDB_ins_code.push(sites.pdbx_PDB_ins_code!.str(idxB));

            k += 1;
        }
    }

    const struct_conn: Partial<CifCategory.Fields<mmCIF_Schema['struct_conn']>> = {
        id: CifField.ofStrings(id),
        conn_type_id: CifField.ofStrings(conn_type_id),

        ptnr1_label_asym_id: CifField.ofStrings(ptnr1_label_asym_id),
        ptnr1_label_seq_id: CifField.ofNumbers(ptnr1_label_seq_id),
        ptnr1_auth_seq_id: CifField.ofNumbers(ptnr1_auth_seq_id),
        ptnr1_label_atom_id: CifField.ofStrings(ptnr1_label_atom_id),
        pdbx_ptnr1_label_alt_id: CifField.ofStrings(ptnr1_label_alt_id),
        pdbx_ptnr1_PDB_ins_code: CifField.ofStrings(ptnr1_PDB_ins_code),

        ptnr2_label_asym_id: CifField.ofStrings(ptnr2_label_asym_id),
        ptnr2_label_seq_id: CifField.ofNumbers(ptnr2_label_seq_id),
        ptnr2_auth_seq_id: CifField.ofNumbers(ptnr2_auth_seq_id),
        ptnr2_label_atom_id: CifField.ofStrings(ptnr2_label_atom_id),
        pdbx_ptnr2_label_alt_id: CifField.ofStrings(ptnr2_label_alt_id),
        pdbx_ptnr2_PDB_ins_code: CifField.ofStrings(ptnr2_PDB_ins_code),
    };

    return CifCategory.ofFields('struct_conn', struct_conn);
}