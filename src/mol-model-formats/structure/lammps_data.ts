/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from '../../mol-data/db';
import { Model } from '../../mol-model/structure/model';
import { LammpDataFile } from '../../mol-io/reader/lammps/parser';
import { Trajectory, ArrayTrajectory } from '../../mol-model/structure';
import { BondType, MoleculeType } from '../../mol-model/structure/model/types';
import { RuntimeContext, Task } from '../../mol-task';
import { ModelFormat } from '../format';
import { createModels } from './basic/parser';
import { BasicSchema, createBasic } from './basic/schema';
import { ComponentBuilder } from './common/component';
import { EntityBuilder } from './common/entity';
import { IndexPairBonds } from './property/bonds/index-pair';
import { AtomPartialCharge } from './property/partial-charge';

async function getModels(mol: LammpDataFile, ctx: RuntimeContext) {
    const { atoms, bonds } = mol;
    const models: Model[] = [];
    const count = atoms.count;

    const type_symbols = new Array<string>(count);
    const id = new Int32Array(count);
    const cx = new Float32Array(count);
    const cy = new Float32Array(count);
    const cz = new Float32Array(count);
    const model_num = new Int32Array(count);

    let offset = 0;
    for (let j = 0; j < count; j++) {
        type_symbols[offset] = 'C'; // atoms.atomType.value(j);
        cx[offset] = atoms.x.value(j);
        cy[offset] = atoms.y.value(j);
        cz[offset] = atoms.z.value(j);
        id[offset] = atoms.atomId.value(j);
        model_num[offset] = 0;
        offset++;
    }


    const MOL = Column.ofConst('MOL', count, Column.Schema.str);
    const A = Column.ofConst('A', count, Column.Schema.str);
    const seq_id = Column.ofConst(1, count, Column.Schema.int);

    const type_symbol = Column.ofStringArray(type_symbols);

    const atom_site = Table.ofPartialColumns(BasicSchema.atom_site, {
        auth_asym_id: A,
        auth_atom_id: type_symbol,
        auth_comp_id: MOL,
        auth_seq_id: seq_id,
        Cartn_x: Column.ofFloatArray(cx),
        Cartn_y: Column.ofFloatArray(cy),
        Cartn_z: Column.ofFloatArray(cz),
        id: Column.ofIntArray(id),

        label_asym_id: A,
        label_atom_id: type_symbol,
        label_comp_id: MOL,
        label_seq_id: seq_id,
        label_entity_id: Column.ofConst('1', count, Column.Schema.str),

        occupancy: Column.ofConst(1, count, Column.Schema.float),
        type_symbol,

        pdbx_PDB_model_num: Column.ofIntArray(model_num),
    }, count);

    const entityBuilder = new EntityBuilder();
    entityBuilder.setNames([['MOL', 'Unknown Entity']]);
    entityBuilder.getEntityId('MOL', MoleculeType.Unknown, 'A');

    const componentBuilder = new ComponentBuilder(seq_id, type_symbol);
    componentBuilder.setNames([['MOL', 'Unknown Molecule']]);
    componentBuilder.add('MOL', 0);

    const basic = createBasic({
        entity: entityBuilder.getEntityTable(),
        chem_comp: componentBuilder.getChemCompTable(),
        atom_site
    });
    const _models = await createModels(basic, LammpsDataFormat.create(mol), ctx);
    if (_models.frameCount > 0) {
        const indexA = Column.ofIntArray(Column.mapToArray(bonds.atomIdA, x => x - 1, Int32Array));
        const indexB = Column.ofIntArray(Column.mapToArray(bonds.atomIdB, x => x - 1, Int32Array));
        const key = bonds.bondId;
        const order = Column.ofIntArray(Column.mapToArray(bonds.bondType, x => {
            switch (x) {
                default:
                    return 1;
            }
        }, Int8Array));
        const flag = Column.ofIntArray(Column.mapToArray(bonds.bondType, x => {
            switch (x) {
                default:
                    return BondType.Flag.Covalent;
            }
        }, Int8Array));
        const pairBonds = IndexPairBonds.fromData(
            { pairs: { key, indexA, indexB, order, flag }, count: atoms.count },
            { maxDistance: Infinity }
        );

        const first = _models.representative;
        IndexPairBonds.Provider.set(first, pairBonds);

        AtomPartialCharge.Provider.set(first, {
            data: atoms.charge,
            type: 'NO_CHARGES'
        });

        models.push(first);
    }
    return new ArrayTrajectory(models);
}

//

export { LammpsDataFormat };

type LammpsDataFormat = ModelFormat<LammpDataFile>

namespace LammpsDataFormat {
    export function is(x?: ModelFormat): x is LammpsDataFormat {
        return x?.kind === 'xyz';
    }

    export function create(mol: LammpDataFile): LammpsDataFormat {
        return { kind: 'data', name: 'data', data: mol };
    }
}

export function trajectoryFromLammpsData(mol: LammpDataFile): Task<Trajectory> {
    return Task.create('Parse Lammps Data', ctx => getModels(mol, ctx));
}
