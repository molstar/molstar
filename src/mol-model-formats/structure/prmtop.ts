/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, Table } from '../../mol-data/db';
import { PrmtopFile } from '../../mol-io/reader/prmtop/parser';
import { getMoleculeType, MoleculeType } from '../../mol-model/structure/model/types';
import { Topology } from '../../mol-model/structure/topology/topology';
import { Task } from '../../mol-task';
import { ModelFormat } from '../format';
import { BasicSchema, createBasic } from './basic/schema';
import { ComponentBuilder } from './common/component';
import { EntityBuilder } from './common/entity';
import { getChainId } from './common/util';
import { guessElementSymbolString } from './util';

function getBasic(prmtop: PrmtopFile) {
    const { pointers, residuePointer, residueLabel, atomName } = prmtop;
    const atomCount = pointers.NATOM;
    const residueCount = pointers.NRES;

    //

    const residueIds = new Uint32Array(atomCount);
    const residueNames: string[] = [];

    const addResidue = (i: number, from: number, to: number) => {
        const rn = residueLabel.value(i);
        for (let j = from, jl = to; j < jl; ++j) {
            residueIds[j] = i + 1;
            residueNames[j] = rn;
        }
    };

    for (let i = 0, il = residueCount - 1; i < il; ++i) {
        addResidue(i, residuePointer.value(i) - 1, residuePointer.value(i + 1) - 1);

    }
    addResidue(residueCount - 1, residuePointer.value(residueCount - 1) - 1, atomCount);

    const residueId = Column.ofIntArray(residueIds);
    const residueName = Column.ofStringArray(residueNames);

    //

    const entityIds = new Array<string>(atomCount);
    const asymIds = new Array<string>(atomCount);
    const seqIds = new Uint32Array(atomCount);
    const ids = new Uint32Array(atomCount);

    const entityBuilder = new EntityBuilder();
    const componentBuilder = new ComponentBuilder(residueId, atomName);

    let currentEntityId = '';
    let currentAsymIndex = 0;
    let currentAsymId = '';
    let currentSeqId = 0;
    let prevMoleculeType = MoleculeType.Unknown;
    let prevResidueNumber = -1;

    for (let i = 0, il = atomCount; i < il; ++i) {
        const residueNumber = residueId.value(i);
        if (residueNumber !== prevResidueNumber) {
            const compId = residueName.value(i);
            const moleculeType = getMoleculeType(componentBuilder.add(compId, i).type, compId);

            if (moleculeType !== prevMoleculeType) {
                currentAsymId = getChainId(currentAsymIndex);
                currentAsymIndex += 1;
                currentSeqId = 0;
            }

            currentEntityId = entityBuilder.getEntityId(compId, moleculeType, currentAsymId);
            currentSeqId += 1;

            prevResidueNumber = residueNumber;
            prevMoleculeType = moleculeType;
        }

        entityIds[i] = currentEntityId;
        asymIds[i] = currentAsymId;
        seqIds[i] = currentSeqId;
        ids[i] = i;
    }

    const id = Column.ofIntArray(ids);
    const asym_id = Column.ofStringArray(asymIds);

    //

    const type_symbol = new Array<string>(atomCount);
    for (let i = 0; i < atomCount; ++i) {
        type_symbol[i] = guessElementSymbolString(atomName.value(i), residueName.value(i));
    }

    const atom_site = Table.ofPartialColumns(BasicSchema.atom_site, {
        auth_asym_id: asym_id,
        auth_atom_id: Column.asArrayColumn(atomName),
        auth_comp_id: residueName,
        auth_seq_id: residueId,
        id: Column.asArrayColumn(id),

        label_asym_id: asym_id,
        label_atom_id: Column.asArrayColumn(atomName),
        label_comp_id: residueName,
        label_seq_id: Column.ofIntArray(seqIds),
        label_entity_id: Column.ofStringArray(entityIds),

        occupancy: Column.ofConst(1, atomCount, Column.Schema.float),
        type_symbol: Column.ofStringArray(type_symbol),

        pdbx_PDB_model_num: Column.ofConst(1, atomCount, Column.Schema.int),
    }, atomCount);

    const basic = createBasic({
        entity: entityBuilder.getEntityTable(),
        chem_comp: componentBuilder.getChemCompTable(),
        atom_site
    });

    return basic;
}

//

export { PrmtopFormat };

type PrmtopFormat = ModelFormat<PrmtopFile>

namespace PrmtopFormat {
    export function is(x?: ModelFormat): x is PrmtopFormat {
        return x?.kind === 'prmtop';
    }

    export function fromPrmtop(prmtop: PrmtopFile): PrmtopFormat {
        return { kind: 'prmtop', name: prmtop.title.join(' ') || 'PRMTOP', data: prmtop };
    }
}

export function topologyFromPrmtop(prmtop: PrmtopFile): Task<Topology> {
    return Task.create('Parse PRMTOP', async ctx => {
        const format = PrmtopFormat.fromPrmtop(prmtop);
        const basic = getBasic(prmtop);

        const { pointers: { NBONH, NBONA }, bondsIncHydrogen, bondsWithoutHydrogen } = prmtop;
        const bondCount = NBONH + NBONA;

        const bonds = {
            indexA: Column.ofLambda({
                value: (row: number) => {
                    return row < NBONH
                        ? bondsIncHydrogen.value(row * 3) / 3
                        : bondsWithoutHydrogen.value((row - NBONH) * 3) / 3;
                },
                rowCount: bondCount,
                schema: Column.Schema.int,
            }),
            indexB: Column.ofLambda({
                value: (row: number) => {
                    return row < NBONH
                        ? bondsIncHydrogen.value(row * 3 + 1) / 3
                        : bondsWithoutHydrogen.value((row - NBONH) * 3 + 1) / 3;
                },
                rowCount: bondCount,
                schema: Column.Schema.int,
            }),
            order: Column.ofConst(1, bondCount, Column.Schema.int)
        };

        return Topology.create(prmtop.title.join(' ') || 'PRMTOP', basic, bonds, format);
    });
}