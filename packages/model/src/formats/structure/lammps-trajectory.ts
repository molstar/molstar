/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Coordinates, type Frame, Time } from '@molstar/model/model/structure/coordinates';
import {
  type LammpsBox,
  type LammpsFrame,
  type LammpsTrajectoryFile,
  lammpsUnitStyles,
  type UnitStyle,
} from '@molstar/io/reader/lammps/schema';
import { Model } from '@molstar/model/model/structure/model';
import { type RuntimeContext, Task } from '@molstar/core/task';
import { Column, Table } from '@molstar/core/data/db';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { Trajectory } from '@molstar/model/model/structure';
import { BasicSchema, createBasic } from './basic/schema.js';
import { ComponentBuilder } from './common/component.js';
import { EntityBuilder } from './common/entity.js';
import { createModels } from './basic/parser.js';
import { MoleculeType } from '@molstar/model/model/structure/model/types';
import type { ModelFormat } from '../format.js';
import { ModelSymmetry } from './property/symmetry.js';
import { lammpsPeriodicBox } from './common/lammps-box.js';

/**
 * A LAMMPS dump may list atoms in a different order in every frame: LAMMPS splits
 * the simulation box across processors (spatial domain decomposition) and writes
 * each processor's atoms as they arrive, so the per-frame order follows how atoms
 * are distributed and migrate between domains rather than their id. That order also
 * generally does not match a separately-loaded topology (e.g. a sorted `.data`
 * file). Each atom row carries its `id`, so we scatter every frame into canonical
 * id order (row for `id` lands at 0-based array slot `id - 1`) to keep all frames
 * consistent with each other and with the topology. Without this, coordinates land
 * on the wrong atoms.
 * (Dumps written with `dump_modify sort id` are already ordered and hit the
 * identity fast-path below.)
 *
 * Returns the source-row -> canonical-index permutation, or `undefined` when no
 * reordering is needed (rows already in id order) or the ids are not a contiguous
 * 1..count set - in both cases callers keep file order. Allocates the permutation
 * only once a row is found out of order, in a single pass (`atomId.value` parses a
 * token, so it is read exactly once per row).
 */
export function getCanonicalOrder(frame: Pick<LammpsFrame, 'count' | 'atomId'>): Int32Array | undefined {
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
      for (let k = 0; k < j; k++) {
        order[k] = k;
        seen[k] = 1;
      }
    }
    if (seen![dst]) return undefined; // seen is set with order
    seen![dst] = 1;
    order[j] = dst;
  }
  return order;
}

/**
 * A frame's coordinates in Å in atom id order (`getCanonicalOrder`): scaled ones (`xs`) through the box and its tilt;
 * with image flags (`ix iy iz`) on wrapped coordinates also the unwrapped ones, `x + ix a + iy b + iz c`: the atoms as
 * they moved, exact (model.md §4 "Periodic frames").
 */
export function lammpsFrameCoordinates(frame: LammpsFrame, box: LammpsBox | undefined, scale: number) {
  const count = frame.count;
  const order = getCanonicalOrder(frame);
  const x = new Float32Array(count);
  const y = new Float32Array(count);
  const z = new Float32Array(count);
  const isScaled = frame.atomMode.includes('s');
  const [lx, ly, lz] = box?.length ?? [1, 1, 1];
  const [xy, xz, yz] = box?.tilt ?? [0, 0, 0];
  const [ox, oy, oz] = box?.lower ?? [0, 0, 0];
  for (let j = 0; j < count; j++) {
    const dst = order ? order[j] : j;
    const fx = frame.x.value(j),
      fy = frame.y.value(j),
      fz = frame.z.value(j);
    if (isScaled && box) {
      // fractional coordinates of the box: corner + xs a + ys b + zs c (a triclinic box's tilt in b and c)
      x[dst] = (ox + fx * lx + fy * xy + fz * xz) * scale;
      y[dst] = (oy + fy * ly + fz * yz) * scale;
      z[dst] = (oz + fz * lz) * scale;
    } else {
      x[dst] = fx * scale;
      y[dst] = fy * scale;
      z[dst] = fz * scale;
    }
  }
  // image flags mean something only on wrapped coordinates (`x`, `xs`; `xu` is unwrapped already) and with a box
  const image = frame.image;
  if (!image || !box || frame.atomMode.includes('u')) return { x, y, z };
  const ux = new Float32Array(count);
  const uy = new Float32Array(count);
  const uz = new Float32Array(count);
  for (let j = 0; j < count; j++) {
    const dst = order ? order[j] : j;
    const na = image.ix.value(j),
      nb = image.iy.value(j),
      nc = image.iz.value(j);
    ux[dst] = x[dst] + (na * lx + nb * xy + nc * xz) * scale;
    uy[dst] = y[dst] + (nb * ly + nc * yz) * scale;
    uz[dst] = z[dst] + nc * lz * scale;
  }
  return { x, y, z, unwrapped: { x: ux, y: uy, z: uz } };
}

export function coordinatesFromLammpsTrajectory(
  file: LammpsTrajectoryFile,
  unitsStyle: UnitStyle = 'real',
): Task<Coordinates> {
  return Task.create('Parse Lammps Trajectory', async (ctx) => {
    await ctx.update('Converting to coordinates');
    const scale = lammpsUnitStyles[unitsStyle].scale;
    const deltaTime = Time(file.deltaTime, 'step');
    const offsetTime = Time(file.timeOffset, deltaTime.unit);
    const frames: Frame[] = [];
    for (let i = 0, il = file.frames.length; i < il; ++i) {
      const box = file.bounds[i];
      const { x, y, z, unwrapped } = lammpsFrameCoordinates(file.frames[i], box, scale);
      // the frame's simulation box (model.md §4 "Periodic frames"): with its corner, which LAMMPS need not put at 0
      const periodic = box ? lammpsPeriodicBox(box, scale) : undefined;
      frames.push({
        elementCount: file.frames[i].count,
        x,
        y,
        z,
        ...(periodic ? { cell: periodic.cell, box: periodic.box } : {}),
        ...(periodic && unwrapped ? { unwrapped: { ...unwrapped, xyzOrdering: { isIdentity: true } } } : {}),
        xyzOrdering: { isIdentity: true },
        time: Time(offsetTime.value + deltaTime.value * i, deltaTime.unit),
      });
    }

    return Coordinates.create(frames, deltaTime, offsetTime);
  });
}

async function getModels(mol: LammpsTrajectoryFile, ctx: RuntimeContext, unitsStyle: UnitStyle = 'real') {
  const atoms = mol.frames[0];
  const count = atoms.count;
  const box = mol.bounds[0];
  const scale = lammpsUnitStyles[unitsStyle].scale;
  // the topology's coordinates as the frames' (scaled ones through the tilt): bonds and ranges are found on them
  const topology = lammpsFrameCoordinates(atoms, box, scale);
  const type_symbols = new Array<string>(count);
  const asym_ids = new Array<string>(count);
  const id = new Int32Array(count);
  // Build the topology in canonical id order so it stays aligned with the
  // per-frame coordinates produced by coordinatesFromLammpsTrajectory.
  const order = getCanonicalOrder(atoms);
  for (let j = 0; j < count; j++) {
    const dst = order ? order[j] : j;
    type_symbols[dst] = atoms.atomType.value(j).toString();
    asym_ids[dst] = atoms.moleculeId.value(j).toString();
    id[dst] = atoms.atomId.value(j);
  }

  const MOL = Column.ofConst('MOL', count, Column.Schema.str);
  const asym_id = Column.ofStringArray(asym_ids);
  const seq_id = Column.ofConst(1, count, Column.Schema.int);

  const type_symbol = Column.ofStringArray(type_symbols);

  const atom_site = Table.ofPartialColumns(
    BasicSchema.atom_site,
    {
      auth_asym_id: asym_id,
      auth_atom_id: type_symbol,
      auth_comp_id: MOL,
      auth_seq_id: seq_id,
      Cartn_x: Column.ofFloatArray(topology.x),
      Cartn_y: Column.ofFloatArray(topology.y),
      Cartn_z: Column.ofFloatArray(topology.z),
      id: Column.ofIntArray(id),

      label_asym_id: asym_id,
      label_atom_id: type_symbol,
      label_comp_id: MOL,
      label_seq_id: seq_id,
      label_entity_id: Column.ofConst('1', count, Column.Schema.str),

      occupancy: Column.ofConst(1, count, Column.Schema.float),
      type_symbol,

      pdbx_PDB_model_num: Column.ofConst(0, count, Column.Schema.int),
    },
    count,
  );

  const entityBuilder = new EntityBuilder();
  entityBuilder.setNames([['MOL', 'Unknown Entity']]);
  entityBuilder.getEntityId('MOL', MoleculeType.Unknown, 'A');

  const componentBuilder = new ComponentBuilder(seq_id, type_symbol);
  componentBuilder.setNames([['MOL', 'Unknown Molecule']]);
  componentBuilder.add('MOL', 0);

  const basic = createBasic({
    entity: entityBuilder.getEntityTable(),
    chem_comp: componentBuilder.getChemCompTable(),
    atom_site,
  });
  const _models = await createModels(basic, LammpsTrajectoryFormat.create(mol), ctx);
  const first = _models.representative;

  if (mol.bounds.length > 0) {
    const periodic = lammpsPeriodicBox(mol.bounds[0], scale);
    const [lx, ly, lz] = mol.bounds[0].length;
    const cell = periodic?.cell;
    ModelSymmetry.Provider.set(
      first,
      cell
        ? ModelSymmetry.fromCell(cell.size, cell.anglesInRadians)
        : ModelSymmetry.fromCell(
            Vec3.create(lx * scale, ly * scale, lz * scale),
            Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2),
          ),
    );
  }

  const coordinates = await coordinatesFromLammpsTrajectory(mol, unitsStyle).runInContext(ctx);
  return Model.trajectoryFromModelAndCoordinates(first, coordinates);
}

export { LammpsTrajectoryFormat };

type LammpsTrajectoryFormat = ModelFormat<LammpsTrajectoryFile>;

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
  return Task.create('Parse Lammps Traj Data', (ctx) => getModels(mol, ctx, unitsStyle));
}
