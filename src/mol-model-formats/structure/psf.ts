/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PsfFile } from '../../mol-io/reader/psf/parser';
import { mmCIF_Schema } from '../../mol-io/reader/cif/schema/mmcif';
import { Column } from '../../mol-data/db';
import { EntityBuilder } from './common/entity';
import { ComponentBuilder } from './common/component';
import { CifCategory, CifField } from '../../mol-io/reader/cif';
import { guessElementSymbolString } from './util';
import { MoleculeType, getMoleculeType } from '../../mol-model/structure/model/types';
import { getChainId } from './common/util';
import { Task } from '../../mol-task';
import { ModelFormat } from './format';
import { Topology } from '../../mol-model/structure/topology/topology';

// TODO: shares most of the code with ./gro.ts#getCategories
function getCategories(atoms: PsfFile['atoms']) {
    const auth_atom_id = CifField.ofColumn(atoms.atomName)
    const auth_comp_id = CifField.ofColumn(atoms.residueName)

    const entityIds = new Array<string>(atoms.count)
    const asymIds = new Array<string>(atoms.count)
    const seqIds = new Uint32Array(atoms.count)
    const ids = new Uint32Array(atoms.count)

    const entityBuilder = new EntityBuilder()
    const componentBuilder = new ComponentBuilder(atoms.residueId, atoms.atomName)

    let currentEntityId = ''
    let currentAsymIndex = 0
    let currentAsymId = ''
    let currentSeqId = 0
    let prevMoleculeType = MoleculeType.Unknown
    let prevResidueNumber = -1

    for (let i = 0, il = atoms.count; i < il; ++i) {
        const residueNumber = atoms.residueId.value(i)
        if (residueNumber !== prevResidueNumber) {
            const compId = atoms.residueName.value(i)
            const moleculeType = getMoleculeType(componentBuilder.add(compId, i).type, compId)

            if (moleculeType !== prevMoleculeType || residueNumber !== prevResidueNumber + 1) {
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
        auth_seq_id: CifField.ofColumn(atoms.residueId),
        B_iso_or_equiv: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.float)),
        Cartn_x: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.float)),
        Cartn_y: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.float)),
        Cartn_z: CifField.ofColumn(Column.Undefined(atoms.count, Column.Schema.float)),
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

function psfToMmCif(psf: PsfFile) {
    const categories = getCategories(psf.atoms)

    return {
        header: psf.id,
        categoryNames: Object.keys(categories),
        categories
    };
}

export function topologyFromPsf(psf: PsfFile): Task<Topology> {
    return Task.create('Parse PSF', async ctx => {
        const label = psf.id
        const cif = psfToMmCif(psf);
        const format = ModelFormat.mmCIF(cif);

        const { atomIdA, atomIdB } = psf.bonds

        const bonds = {
            indexA: Column.ofLambda({
                value: (row: number) => atomIdA.value(row) - 1,
                rowCount: atomIdA.rowCount,
                schema: atomIdA.schema,
            }),
            indexB: Column.ofLambda({
                value: (row: number) => atomIdB.value(row) - 1,
                rowCount: atomIdB.rowCount,
                schema: atomIdB.schema,
            }),
            order: Column.ofConst(1, psf.bonds.count, Column.Schema.int)
        }

        return Topology.create(label, format, bonds)
    })
}