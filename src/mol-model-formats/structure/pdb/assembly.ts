/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifCategory, CifField } from 'mol-io/reader/cif';
import { mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif';

export function parseCryst1(id: string, record: string): CifCategory[] {
    // COLUMNS       DATA TYPE      CONTENTS
    // --------------------------------------------------------------------------------
    //  1 -  6       Record name    "CRYST1"
    //  7 - 15       Real(9.3)      a (Angstroms)
    // 16 - 24       Real(9.3)      b (Angstroms)
    // 25 - 33       Real(9.3)      c (Angstroms)
    // 34 - 40       Real(7.2)      alpha (degrees)
    // 41 - 47       Real(7.2)      beta (degrees)
    // 48 - 54       Real(7.2)      gamma (degrees)
    // 56 - 66       LString        Space group
    // 67 - 70       Integer        Z value

    const get = (s: number, l: number) => (record.substr(s, l) || '').trim()

    const cell: CifCategory.Fields<mmCIF_Schema['cell']> = {
        entry_id: CifField.ofString(id),
        length_a: CifField.ofString(get(6, 9)),
        length_b: CifField.ofString(get(15, 9)),
        length_c: CifField.ofString(get(24, 9)),
        angle_alpha: CifField.ofString(get(33, 7)),
        angle_beta: CifField.ofString(get(40, 7)),
        angle_gamma: CifField.ofString(get(47, 7)),
        Z_PDB: CifField.ofString(get(66, 4)),
        pdbx_unique_axis: CifField.ofString('?')
    };
    const symmetry: CifCategory.Fields<mmCIF_Schema['symmetry']> = {
        entry_id: CifField.ofString(id),
        'space_group_name_H-M': CifField.ofString(get(55, 11)),
        Int_Tables_number: CifField.ofString('?'),
        cell_setting: CifField.ofString('?'),
        space_group_name_Hall: CifField.ofString('?')
    }
    return [CifCategory.ofFields('cell', cell), CifCategory.ofFields('symmetry', symmetry)];
}