/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Model } from '@molstar/model/model/structure/model';
import { Task } from '@molstar/core/task';
import type { ModelFormat } from '../format.js';
import type { GroFile, GroAtoms } from '@molstar/io/reader/gro/schema';
import { Column, Table } from '@molstar/core/data/db';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { guessElementSymbolString } from './util.js';
import { MoleculeType, getMoleculeType } from '@molstar/model/model/structure/model/types';
import { ComponentBuilder } from './common/component.js';
import { getChainId } from './common/util.js';
import { EntityBuilder } from './common/entity.js';
import { type BasicData, BasicSchema, createBasic } from './basic/schema.js';
import { createModels } from './basic/parser.js';
import type { Trajectory } from '@molstar/model/model/structure';
import { ArrayTrajectory } from '@molstar/model/model/structure/trajectory';
import { ModelSymmetry } from './property/symmetry.js';
import { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import { PeriodicBox } from '@molstar/core/math/geometry/periodic-box';

function getBasic(atoms: GroAtoms, modelNum: number): BasicData {
  const auth_atom_id = atoms.atomName;
  const auth_comp_id = atoms.residueName;

  const entityIds = new Array<string>(atoms.count);
  const asymIds = new Array<string>(atoms.count);
  const seqIds = new Uint32Array(atoms.count);
  const ids = new Uint32Array(atoms.count);
  const typeSymbol = new Array<string>(atoms.count);

  const entityBuilder = new EntityBuilder();
  const componentBuilder = new ComponentBuilder(atoms.residueNumber, atoms.atomName);

  let currentEntityId = '';
  let currentAsymIndex = 0;
  let currentAsymId = '';
  let currentSeqId = 0;
  let prevMoleculeType = MoleculeType.Unknown;
  let prevResidueNumber = -1;

  for (let i = 0, il = atoms.count; i < il; ++i) {
    const residueNumber = atoms.residueNumber.value(i);
    if (residueNumber !== prevResidueNumber) {
      const compId = atoms.residueName.value(i);
      const moleculeType = getMoleculeType(componentBuilder.add(compId, i).type, compId);

      if (
        moleculeType !== prevMoleculeType ||
        (residueNumber !== prevResidueNumber + 1 &&
          !(
            // gro format allows only for 5 character residueNumbers, handle overflow here
            prevResidueNumber === 99999 && residueNumber === 0
          ))
      ) {
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

    typeSymbol[i] = guessElementSymbolString(atoms.atomName.value(i), atoms.residueName.value(i));
  }

  const auth_asym_id = Column.ofStringArray(asymIds);

  const atom_site = Table.ofPartialColumns(
    BasicSchema.atom_site,
    {
      auth_asym_id,
      auth_atom_id,
      auth_comp_id,
      auth_seq_id: atoms.residueNumber,
      Cartn_x: Column.ofFloatArray(Column.mapToArray(atoms.x, (x) => x * 10, Float32Array)),
      Cartn_y: Column.ofFloatArray(Column.mapToArray(atoms.y, (y) => y * 10, Float32Array)),
      Cartn_z: Column.ofFloatArray(Column.mapToArray(atoms.z, (z) => z * 10, Float32Array)),
      id: Column.ofIntArray(ids),

      label_asym_id: auth_asym_id,
      label_atom_id: auth_atom_id,
      label_comp_id: auth_comp_id,
      label_seq_id: Column.ofIntArray(seqIds),
      label_entity_id: Column.ofStringArray(entityIds),

      occupancy: Column.ofConst(1, atoms.count, Column.Schema.float),
      type_symbol: Column.ofStringArray(typeSymbol),

      pdbx_PDB_model_num: Column.ofConst(modelNum, atoms.count, Column.Schema.int),
    },
    atoms.count,
  );

  return createBasic({
    entity: entityBuilder.getEntityTable(),
    chem_comp: componentBuilder.getChemCompTable(),
    atom_site,
  });
}

//

export { GroFormat };

type GroFormat = ModelFormat<GroFile>;

namespace GroFormat {
  export function is(x?: ModelFormat): x is GroFormat {
    return x?.kind === 'gro';
  }

  export function fromGro(gro: GroFile): GroFormat {
    return { kind: 'gro', name: gro.structures[0].header.title, data: gro };
  }
}

// TODO: reuse static model parts when hierarchy is identical
//       need to pass all gro.structures as one table into createModels

export function trajectoryFromGRO(gro: GroFile): Task<Trajectory> {
  return Task.create('Parse GRO', async (ctx) => {
    const format = GroFormat.fromGro(gro);
    const models: Model[] = [];
    for (let i = 0, il = gro.structures.length; i < il; ++i) {
      const basic = getBasic(gro.structures[i].atoms, i + 1);
      const m = await createModels(basic, format, ctx);
      if (m.frameCount === 1) {
        const model = m.representative;
        // nm to Å; a triclinic box has its off-diagonal values (the rhombic dodecahedron, the truncated octahedron)
        const v = gro.structures[i].header.boxVectors;
        const a = Vec3.create(v[0] * 10, v[1] * 10, v[2] * 10);
        const b = Vec3.create(v[3] * 10, v[4] * 10, v[5] * 10);
        const c = Vec3.create(v[6] * 10, v[7] * 10, v[8] * 10);
        const box = PeriodicBox.create(a, b, c);
        if (PeriodicBox.isValid(box)) {
          const cell = Cell.fromBasis(a, b, c);
          ModelSymmetry.Provider.set(model, ModelSymmetry.fromCell(cell.size, cell.anglesInRadians));
          Model.PeriodicBox.set(model, box);
        }
        models.push(model);
      }
    }
    return new ArrayTrajectory(models);
  });
}
