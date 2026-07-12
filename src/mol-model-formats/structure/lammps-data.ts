/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Column, Table } from '../../mol-data/db';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Model } from '../../mol-model/structure/model';
import { LammpsDataFile, lammpsUnitStyles, UnitStyle } from '../../mol-io/reader/lammps/schema';
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
import { ModelSymmetry } from './property/symmetry';

async function getModels(mol: LammpsDataFile, ctx: RuntimeContext, unitsStyle: UnitStyle = 'real') {
    const { atoms, bonds } = mol;
    const models: Model[] = [];
    const count = atoms.count;
    const scale = lammpsUnitStyles[unitsStyle].scale;
    const type_symbols = new Array<string>(count);
    const id = new Int32Array(count);
    const cx = new Float32Array(count);
    const cy = new Float32Array(count);
    const cz = new Float32Array(count);

    // A LAMMPS `.data` file may list atoms in any order (rows are not necessarily sorted by atom
    // id), while the Bonds section references atoms by id. Record each id's row so bonds connect
    // the correct atoms. Atom ids are positive and, for a data file, on the order of 1..count, so
    // an array indexed by id (with a -1 sentinel for holes) is smaller and faster than a Map; this
    // deliberately assumes `maxId` stays close to `count` (true for a real data file).
    let maxId = 0;
    for (let j = 0; j < count; j++) {
        const atomId = atoms.atomId.value(j);
        type_symbols[j] = atoms.atomType.value(j).toString();
        cx[j] = atoms.x.value(j) * scale;
        cy[j] = atoms.y.value(j) * scale;
        cz[j] = atoms.z.value(j) * scale;
        id[j] = atomId;
        if (atomId > maxId) maxId = atomId;
    }
    const rowOfId = new Int32Array(maxId + 1).fill(-1);
    for (let j = 0; j < count; j++) rowOfId[id[j]] = j;

    const MOL = Column.ofConst('MOL', count, Column.Schema.str);
    const asym_id = Column.ofLambda({
        value: (row: number) => atoms.moleculeId.value(row).toString(),
        rowCount: count,
        schema: Column.Schema.str,
    });
    const seq_id = Column.ofConst(1, count, Column.Schema.int);

    const type_symbol = Column.ofStringArray(type_symbols);

    const atom_site = Table.ofPartialColumns(BasicSchema.atom_site, {
        auth_asym_id: asym_id,
        auth_atom_id: type_symbol,
        auth_comp_id: MOL,
        auth_seq_id: seq_id,
        Cartn_x: Column.ofFloatArray(cx),
        Cartn_y: Column.ofFloatArray(cy),
        Cartn_z: Column.ofFloatArray(cz),
        id: Column.ofIntArray(id),

        label_asym_id: asym_id,
        label_atom_id: type_symbol,
        label_comp_id: MOL,
        label_seq_id: seq_id,
        label_entity_id: Column.ofConst('1', count, Column.Schema.str),

        occupancy: Column.ofConst(1, count, Column.Schema.float),
        type_symbol,

        pdbx_PDB_model_num: Column.ofConst(0, count, Column.Schema.int),
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
        const first = _models.representative;
        if (bonds.count !== 0) {
            // resolve bond atom ids to atom rows via `rowOfId`, dropping a bond that references an
            // unknown atom id (past the table, or a -1 hole in the id range) rather than silently
            // pointing it at row 0
            const ia: number[] = [], ib: number[] = [], keys: number[] = [];
            for (let i = 0; i < bonds.count; i++) {
                const idA = bonds.atomIdA.value(i), idB = bonds.atomIdB.value(i);
                if (idA > maxId || idB > maxId) continue;
                const a = rowOfId[idA], b = rowOfId[idB];
                if (a < 0 || b < 0) continue;
                ia.push(a); ib.push(b); keys.push(bonds.bondId.value(i));
            }
            const indexA = Column.ofIntArray(new Int32Array(ia));
            const indexB = Column.ofIntArray(new Int32Array(ib));
            const key = Column.ofIntArray(new Int32Array(keys));
            const order = Column.ofConst(1, ia.length, Column.Schema.int);
            const flag = Column.ofConst(BondType.Flag.Covalent, ia.length, Column.Schema.int);
            const pairBonds = IndexPairBonds.fromData(
                { pairs: { key, indexA, indexB, order, flag }, count: atoms.count },
                { maxDistance: Infinity }
            );
            IndexPairBonds.Provider.set(first, pairBonds);
        }

        AtomPartialCharge.Provider.set(first, {
            data: atoms.charge,
            type: 'NO_CHARGES'
        });

        if (mol.box) {
            const [lx, ly, lz] = mol.box.length;
            const symmetry = ModelSymmetry.fromCell(
                Vec3.create(lx * scale, ly * scale, lz * scale),
                Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2)
            );
            ModelSymmetry.Provider.set(first, symmetry);
        }

        models.push(first);
    }
    return new ArrayTrajectory(models);
}

//

export { LammpsDataFormat };

type LammpsDataFormat = ModelFormat<LammpsDataFile>

namespace LammpsDataFormat {
    export function is(x?: ModelFormat): x is LammpsDataFormat {
        return x?.kind === 'data';
    }

    export function create(mol: LammpsDataFile): LammpsDataFormat {
        return { kind: 'data', name: 'data', data: mol };
    }
}

export function trajectoryFromLammpsData(mol: LammpsDataFile, unitsStyle?: UnitStyle): Task<Trajectory> {
    if (unitsStyle === void 0) unitsStyle = 'real';
    return Task.create('Parse Lammps Data', ctx => getModels(mol, ctx, unitsStyle));
}
