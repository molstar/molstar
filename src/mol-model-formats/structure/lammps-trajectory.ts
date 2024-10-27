/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Coordinates, Frame, Time } from '../../mol-model/structure/coordinates';
import { LammpsTrajectoryFile, lammpsUnitStyles, UnitStyle } from '../../mol-io/reader/lammps/schema';
import { Model } from '../../mol-model/structure/model';
import { RuntimeContext, Task } from '../../mol-task';
import { Column, Table } from '../../mol-data/db';
import { Trajectory } from '../../mol-model/structure';
import { BasicSchema, createBasic } from './basic/schema';
import { ComponentBuilder } from './common/component';
import { EntityBuilder } from './common/entity';
import { createModels } from './basic/parser';
import { MoleculeType } from '../../mol-model/structure/model/types';
import { ModelFormat } from '../format';

export function coordinatesFromLammpsTrajectory(file: LammpsTrajectoryFile, unitsStyle: UnitStyle = 'real'): Task<Coordinates> {
    return Task.create('Parse Lammps Trajectory', async ctx => {
        await ctx.update('Converting to coordinates');
        const scale = lammpsUnitStyles[unitsStyle].scale;
        const deltaTime = Time(file.deltaTime, 'step');
        const offsetTime = Time(file.timeOffset, deltaTime.unit);
        const offset_pos = { x: 0.0, y: 0.0, z: 0.0 };
        const offset_scale = { x: 1.0, y: 1.0, z: 1.0 };
        const atomsMode = file.frames[0].atomMode;
        const isScaled = atomsMode.includes('s');
        const frames: Frame[] = [];
        for (let i = 0, il = file.frames.length; i < il; ++i) {
            const box = file.bounds[i];
            if (isScaled) {
                offset_scale.x = box.length[0];
                offset_scale.y = box.length[1];
                offset_scale.z = box.length[2];
                offset_pos.x = box.lower[0];
                offset_pos.y = box.lower[1];
                offset_pos.z = box.lower[2];
            }
            const count = file.frames[i].count;
            const cx = new Float32Array(count);
            const cy = new Float32Array(count);
            const cz = new Float32Array(count);
            let offset = 0;
            for (let j = 0; j < count; j++) {
                cx[offset] = (file.frames[i].x.value(j) * offset_scale.x + offset_pos.x) * scale;
                cy[offset] = (file.frames[i].y.value(j) * offset_scale.y + offset_pos.y) * scale;
                cz[offset] = (file.frames[i].z.value(j) * offset_scale.z + offset_pos.z) * scale;
                offset++;
            }
            frames.push({
                elementCount: file.frames[i].count,
                x: cx,
                y: cy,
                z: cz,
                xyzOrdering: { isIdentity: true },
                time: Time(offsetTime.value + deltaTime.value * i, deltaTime.unit)
            });
        }

        return Coordinates.create(frames, deltaTime, offsetTime);
    });
}

async function getModels(mol: LammpsTrajectoryFile, ctx: RuntimeContext, unitsStyle: UnitStyle = 'real') {
    const atoms = mol.frames[0];
    const count = atoms.count;
    const atomsMode = atoms.atomMode;
    const box = mol.bounds[0];
    const offset_pos = { x: 0.0, y: 0.0, z: 0.0 };
    const offset_scale = { x: 1.0, y: 1.0, z: 1.0 };
    const scale = lammpsUnitStyles[unitsStyle].scale;
    // if caracter s in atomsMode, we need to scale the coordinates
    if (atomsMode.includes('s')) {
        offset_scale.x = box.length[0];
        offset_scale.y = box.length[1];
        offset_scale.z = box.length[2];
        offset_pos.x = box.lower[0];
        offset_pos.y = box.lower[1];
        offset_pos.z = box.lower[2];
    }
    const type_symbols = new Array<string>(count);
    const id = new Int32Array(count);
    const cx = new Float32Array(count);
    const cy = new Float32Array(count);
    const cz = new Float32Array(count);
    const model_num = new Int32Array(count);

    let offset = 0;
    for (let j = 0; j < count; j++) {
        type_symbols[offset] = atoms.atomType.value(j).toString();
        cx[offset] = (atoms.x.value(j) * offset_scale.x + offset_pos.x) * scale;
        cy[offset] = (atoms.y.value(j) * offset_scale.y + offset_pos.y) * scale;
        cz[offset] = (atoms.z.value(j) * offset_scale.z + offset_pos.z) * scale;
        id[offset] = atoms.atomId.value(j);
        model_num[offset] = 0;
        offset++;
    }

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
    const _models = await createModels(basic, LammpsTrajectoryFormat.create(mol), ctx);
    const first = _models.representative;
    const coordinates = await coordinatesFromLammpsTrajectory(mol, unitsStyle).runInContext(ctx);
    return Model.trajectoryFromModelAndCoordinates(first, coordinates);
}

export { LammpsTrajectoryFormat };

type LammpsTrajectoryFormat = ModelFormat<LammpsTrajectoryFile>

namespace LammpsTrajectoryFormat {
    export function is(x?: ModelFormat): x is LammpsTrajectoryFormat {
        return x?.kind === 'lammpstrj';
    }

    export function create(mol: LammpsTrajectoryFile): LammpsTrajectoryFormat {
        return { kind: 'lammpstrj', name: 'lammpstrj', data: mol };
    }
}

export function trajectoryFromLammpsTrajectory(mol: LammpsTrajectoryFile, unitsStyle?: UnitStyle): Task<Trajectory> {
    if (unitsStyle === void 0) unitsStyle = 'real';
    return Task.create('Parse Lammps Traj Data', ctx => getModels(mol, ctx, unitsStyle));
}
