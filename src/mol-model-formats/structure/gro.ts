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
import { MoleculeType, getMoleculeType } from '../../mol-model/structure/model/types';
import { ComponentBuilder } from './common/component';
import { getChainId } from './common/util';
import { EntityBuilder } from './common/entity';

// TODO multi model files

function getCategories(atoms: GroAtoms) {
    const auth_atom_id = CifField.ofColumn(atoms.atomName)
    const auth_comp_id = CifField.ofColumn(atoms.residueName)

    const entityIds = new Array<string>(atoms.count)
    const asymIds = new Array<string>(atoms.count)
    const seqIds = new Uint32Array(atoms.count)
    const ids = new Uint32Array(atoms.count)

    const entityBuilder = new EntityBuilder()
    const componentBuilder = new ComponentBuilder(atoms.residueNumber, atoms.atomName)

    let currentEntityId = ''
    let currentAsymIndex = 0
    let currentAsymId = ''
    let currentSeqId = 0
    let prevMoleculeType = MoleculeType.unknown
    let prevResidueNumber = -1

    for (let i = 0, il = atoms.count; i < il; ++i) {
        const residueNumber = atoms.residueNumber.value(i)
        if (residueNumber !== prevResidueNumber) {
            const compId = atoms.residueName.value(i)
            const moleculeType = getMoleculeType(componentBuilder.add(compId, i).type, compId)

            if (moleculeType !== prevMoleculeType || (
                residueNumber !== prevResidueNumber + 1 && !(
                    // gro format allows only for 5 character residueNumbers, handle overflow here
                    prevResidueNumber === 99999 && residueNumber === 0
                )
            )) {
                currentAsymId = getChainId(currentAsymIndex)
                currentAsymIndex += 1
                currentSeqId = 0
            }

            currentEntityId = entityBuilder.getEntityId(compId, moleculeType, currentAsymId)
            currentSeqId += 1

            prevResidueNumber = residueNumber
            prevMoleculeType = moleculeType
        }

        entityIds[i] = currentEntityId
        asymIds[i] = currentAsymId
        seqIds[i] = currentSeqId
        ids[i] = i
    }

    const auth_asym_id = CifField.ofColumn(Column.ofStringArray(asymIds))

    const atom_site: CifCategory.SomeFields<mmCIF_Schema['atom_site']> = {
        auth_asym_id,
        auth_atom_id,
        auth_comp_id,
        auth_seq_id: CifField.ofColumn(atoms.residueNumber),
        B_iso_or_equiv: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.float)),
        Cartn_x: CifField.ofNumbers(Column.mapToArray(atoms.x, x => x * 10, Float32Array)),
        Cartn_y: CifField.ofNumbers(Column.mapToArray(atoms.y, y => y * 10, Float32Array)),
        Cartn_z: CifField.ofNumbers(Column.mapToArray(atoms.z, z => z * 10, Float32Array)),
        group_PDB: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.str)),
        id: CifField.ofColumn(Column.ofIntArray(ids)),

        label_alt_id: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.str)),

        label_asym_id: auth_asym_id,
        label_atom_id: auth_atom_id,
        label_comp_id: auth_comp_id,
        label_seq_id: CifField.ofColumn(Column.ofIntArray(seqIds)),
        label_entity_id: CifField.ofColumn(Column.ofStringArray(entityIds)),

        occupancy: CifField.ofColumn(Column.ofConst(1, atoms.count, Column.Schema.float)),
        type_symbol: CifField.ofStrings(Column.mapToArray(atoms.atomName, s => guessElementSymbolString(s))),

        pdbx_PDB_ins_code: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.str)),
        pdbx_PDB_model_num: CifField.ofColumn(Column.ofConst('1', atoms.count, Column.Schema.str)),
    }

    return {
        entity: entityBuilder.getEntityCategory(),
        chem_comp: componentBuilder.getChemCompCategory(),
        atom_site: CifCategory.ofFields('atom_site', atom_site)
    }
}

async function groToMmCif(gro: GroFile) {
    const categories = getCategories(gro.structures[0].atoms)

    return {
        header: gro.structures[0].header.title,
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
