/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Column, Table } from '../../mol-data/db';
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

async function getModels(mol: LammpsDataFile, ctx: RuntimeContext, unitsStyle: UnitStyle = 'real') {
    const { atoms, bonds, bounds } = mol;
    const models: Model[] = [];
    const count = atoms.count;
    const scale = lammpsUnitStyles[unitsStyle].scale;
    const type_symbols = new Array<string>(count);
    const id = new Int32Array(count);
    const cx = new Float32Array(count);
    const cy = new Float32Array(count);
    const cz = new Float32Array(count);
    const model_num = new Int32Array(count);

    const ll = bounds.length;
    const lo = bounds.lower;
    let offset = 0;
    for (let j = 0; j < count; j++) {
        type_symbols[offset] = atoms.atomType.value(j).toString();
        cx[offset] = ((atoms.x.value(j) - lo[0]) % ll[0] + lo[0]) * scale;
        cy[offset] = ((atoms.y.value(j) - lo[1]) % ll[1] + lo[1]) * scale;
        cz[offset] = ((atoms.z.value(j) - lo[2]) % ll[2] + lo[2]) * scale;
        id[offset] = atoms.atomId.value(j) - 1;
        model_num[offset] = 0;
        offset++;
    }

    const MOL = Column.ofConst('MOL', count, Column.Schema.str);
    const asym_id = (count > 200000) ? Column.ofConst('A', count, Column.Schema.str)
        : Column.ofLambda({
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
        const first = _models.representative;
        if (bonds.count !== 0) {
            const indexA = Column.ofIntArray(Column.mapToArray(bonds.atomIdA, x => x - 1, Int32Array));
            const indexB = Column.ofIntArray(Column.mapToArray(bonds.atomIdB, x => x - 1, Int32Array));
            const key = bonds.bondId;
            const order = Column.ofConst(1, bonds.count, Column.Schema.int);
            const flag = Column.ofConst(BondType.Flag.Covalent, bonds.count, Column.Schema.int);
            /*
            const pairBonds = IndexPairBonds.fromData(
                { pairs: { key, indexA, indexB, order, flag }, count: atoms.count },
                { maxDistance: Infinity }
            );
            IndexPairBonds.Provider.set(first, pairBonds);*/
            // Assuming 'atoms' contains the coordinates after wrapping
            const newBonds = [];
            const pos1 = [0, 0, 0];
            const pos2 = [0, 0, 0];
            for (let i = 0; i < bonds.count; i++) {
                const [atom1Index, atom2Index] = [indexA.value(i), indexB.value(i)];

                // Here, you would ideally have access to the wrapped positions of atom1 and atom2
                // If not, you might need to apply wrapping logic here:
                pos1[0] = atoms.x.value(atom1Index);
                pos1[1] = atoms.y.value(atom1Index);
                pos1[2] = atoms.z.value(atom1Index);
                pos2[0] = atoms.x.value(atom2Index);
                pos2[1] = atoms.y.value(atom2Index);
                pos2[2] = atoms.z.value(atom2Index);
                let crossBoundary = false;

                // Check periodicity and apply wrapping if necessary (pseudo-code, adjust according to your actual structure)
                for (let dim = 0; dim < 3; dim++) {
                    // Wrapping logic as previously described
                    pos1[dim] = ((pos1[dim] - bounds.lower[dim]) % bounds.length[dim] + bounds.length[dim]) % bounds.length[dim] + bounds.lower[dim];
                    pos2[dim] = ((pos2[dim] - bounds.lower[dim]) % bounds.length[dim] + bounds.length[dim]) % bounds.length[dim] + bounds.lower[dim];
                    if (Math.abs(pos1[dim] - pos2[dim]) > bounds.length[dim] / 2) {
                        crossBoundary = true;
                    }
                }

                // Decide whether to keep this bond based on whether it crosses the boundary
                if (!crossBoundary) {
                    newBonds.push({
                        key: key.value(i),
                        indexA: atom1Index,
                        indexB: atom2Index,
                        order: order.value(i),
                        flag: flag.value(i)
                    });
                }
            }

            const pairBonds = IndexPairBonds.fromData(
                {
                    pairs: {
                        key: Column.ofIntArray(newBonds.map(b => b.key)),
                        indexA: Column.ofIntArray(newBonds.map(b => b.indexA)),
                        indexB: Column.ofIntArray(newBonds.map(b => b.indexB)),
                        order: Column.ofIntArray(newBonds.map(b => b.order)),
                        flag: Column.ofIntArray(newBonds.map(b => b.flag))
                    },
                    count: atoms.count
                },
                { maxDistance: Infinity }
            );
            IndexPairBonds.Provider.set(first, pairBonds);
        }

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
