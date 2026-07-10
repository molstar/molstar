/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Coordinates, Frame, Time } from '../../mol-model/structure/coordinates';
import { LammpsFrame, LammpsTrajectoryFile, lammpsUnitStyles, UnitStyle } from '../../mol-io/reader/lammps/schema';
import { Model } from '../../mol-model/structure/model';
import { RuntimeContext, Task } from '../../mol-task';
import { Column, Table } from '../../mol-data/db';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Trajectory } from '../../mol-model/structure';
import { BasicSchema, createBasic } from './basic/schema';
import { ComponentBuilder } from './common/component';
import { EntityBuilder } from './common/entity';
import { createModels } from './basic/parser';
import { MoleculeType } from '../../mol-model/structure/model/types';
import { ModelFormat } from '../format';
import { ModelSymmetry } from './property/symmetry';

/**
 * A LAMMPS dump may list atoms in a different order in every frame (atoms are
 * reordered as they migrate between MPI domains), and that order generally does
 * not match a separately-loaded topology (e.g. a sorted `.data` file). Each atom
 * row carries its `id`, so we scatter every frame into canonical id order
 * (destination index = id - 1) to keep all frames consistent with each other and
 * with the topology. Without this, coordinates land on the wrong atoms.
 *
 * Returns the source-row -> canonical-index permutation, or `undefined` when no
 * reordering is needed (rows already in id order) or the ids are not a contiguous
 * 1..count set - in both cases callers keep file order. Allocates the permutation
 * only once a row is found out of order, in a single pass (`atomId.value` parses a
 * token, so it is read exactly once per row).
 */
function getCanonicalOrder(frame: LammpsFrame): Int32Array | undefined {
    const { count, atomId } = frame;
    let order: Int32Array | undefined;
    let seen: Uint8Array | undefined;
    for (let j = 0; j < count; j++) {
        const dst = atomId.value(j) - 1;
        if (dst < 0 || dst >= count) return undefined;
        if (!order) {
            // in order so far (id === row + 1): nothing to permute yet
            if (dst === j) continue;
            // out of order, but dst falls in the already-filled prefix [0, j): duplicate id, reject
            if (dst < j) return undefined;
            // first out-of-order row: allocate and backfill the identity prefix [0, j)
            order = new Int32Array(count);
            seen = new Uint8Array(count);
            for (let k = 0; k < j; k++) { order[k] = k; seen[k] = 1; }
        }
        if (seen![dst]) return undefined; // seen is set with order
        seen![dst] = 1;
        order[j] = dst;
    }
    return order;
}

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
            const frame = file.frames[i];
            const count = frame.count;
            const order = getCanonicalOrder(frame);
            const cx = new Float32Array(count);
            const cy = new Float32Array(count);
            const cz = new Float32Array(count);
            // fold the loop-invariant `* scale` into the per-frame scale/offset factors
            const sx = offset_scale.x * scale, sy = offset_scale.y * scale, sz = offset_scale.z * scale;
            const ox = offset_pos.x * scale, oy = offset_pos.y * scale, oz = offset_pos.z * scale;
            for (let j = 0; j < count; j++) {
                const dst = order ? order[j] : j;
                cx[dst] = frame.x.value(j) * sx + ox;
                cy[dst] = frame.y.value(j) * sy + oy;
                cz[dst] = frame.z.value(j) * sz + oz;
            }
            frames.push({
                elementCount: count,
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
    const molecule_ids = new Int32Array(count);

    // fold the loop-invariant `* scale` into the scale/offset factors
    const sx = offset_scale.x * scale, sy = offset_scale.y * scale, sz = offset_scale.z * scale;
    const ox = offset_pos.x * scale, oy = offset_pos.y * scale, oz = offset_pos.z * scale;

    // Build the topology in canonical id order so it stays aligned with the
    // per-frame coordinates produced by coordinatesFromLammpsTrajectory.
    const order = getCanonicalOrder(atoms);
    for (let j = 0; j < count; j++) {
        const dst = order ? order[j] : j;
        type_symbols[dst] = atoms.atomType.value(j).toString();
        cx[dst] = atoms.x.value(j) * sx + ox;
        cy[dst] = atoms.y.value(j) * sy + oy;
        cz[dst] = atoms.z.value(j) * sz + oz;
        id[dst] = atoms.atomId.value(j);
        molecule_ids[dst] = atoms.moleculeId.value(j);
        model_num[dst] = 0;
    }

    const MOL = Column.ofConst('MOL', count, Column.Schema.str);
    const asym_id = Column.ofLambda({
        value: (row: number) => molecule_ids[row].toString(),
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

    if (mol.bounds.length > 0) {
        const bounds = mol.bounds[0];
        const [lx, ly, lz] = bounds.length;
        const symmetry = ModelSymmetry.fromCell(
            Vec3.create(lx * scale, ly * scale, lz * scale),
            Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2)
        );
        ModelSymmetry.Provider.set(first, symmetry);
    }

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
