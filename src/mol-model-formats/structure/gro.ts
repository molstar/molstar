/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../mol-model/structure/model';
import { Task } from '../../mol-task';
import { ModelFormat } from './format';
import { _parse_mmCif } from './mmcif/parser';
import { GroFile, GroAtoms } from '../../mol-io/reader/gro/schema';
import { CifCategory, CifField } from '../../mol-io/reader/cif';
import { Column } from '../../mol-data/db';
import { mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { guessElementSymbolString } from './util';

// TODO multi model files
// TODO seperate chains
// TODO better entity handling
// TODO improve performance

function _entity(): { [K in keyof mmCIF_Schema['entity']]?: CifField } {
    return {
        id: CifField.ofStrings(['1', '2', '3']),
        type: CifField.ofStrings(['polymer', 'non-polymer', 'water'])
    }
}

function _atom_site(atoms: GroAtoms): { [K in keyof mmCIF_Schema['atom_site']]?: CifField } {
    const auth_asym_id = CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.str))
    const auth_atom_id = CifField.ofColumn(atoms.atomName)
    const auth_comp_id = CifField.ofColumn(atoms.residueName)
    const auth_seq_id = CifField.ofColumn(atoms.residueNumber)

    return {
        auth_asym_id,
        auth_atom_id,
        auth_comp_id,
        auth_seq_id,
        B_iso_or_equiv: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.float)),
        Cartn_x: CifField.ofNumbers(Column.mapToArray(atoms.x, x => x * 10, Float32Array)),
        Cartn_y: CifField.ofNumbers(Column.mapToArray(atoms.y, y => y * 10, Float32Array)),
        Cartn_z: CifField.ofNumbers(Column.mapToArray(atoms.z, z => z * 10, Float32Array)),
        group_PDB: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.str)),
        id: CifField.ofColumn(atoms.atomNumber),

        label_alt_id: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.str)),

        label_asym_id: auth_asym_id,
        label_atom_id: auth_atom_id,
        label_comp_id: auth_comp_id,
        label_seq_id: auth_seq_id,
        label_entity_id: CifField.ofColumn(Column.ofConst('1', atoms.count, Column.Schema.str)),

        occupancy: CifField.ofColumn(Column.ofConst(1, atoms.count, Column.Schema.float)),
        type_symbol: CifField.ofStrings(Column.mapToArray(atoms.atomName, s => guessElementSymbolString(s))),
        // type_symbol: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.str)),

        pdbx_PDB_ins_code: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.str)),
        pdbx_PDB_model_num: CifField.ofColumn(Column.ofConst('1', atoms.count, Column.Schema.str)),
    }
}

async function groToMmCif(gro: GroFile) {
    const categories = {
        entity: CifCategory.ofFields('entity', _entity()),
        atom_site: CifCategory.ofFields('atom_site', _atom_site(gro.structures[0].atoms))
    } as any;

    return {
        header: 'GRO',
        categoryNames: Object.keys(categories),
        categories
    };
}

export function trajectoryFromGRO(gro: GroFile): Task<Model.Trajectory> {
    return Task.create('Parse GRO', async ctx => {
        await ctx.update('Converting to mmCIF');
        const cif = await groToMmCif(gro);
        const format = ModelFormat.mmCIF(cif);
        return _parse_mmCif(format, ctx);
    })
}
