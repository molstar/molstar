/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Model, ResidueIndex, ElementIndex } from './model.js';
import {
  MoleculeType,
  type AtomRole,
  PolymerTypeAtomRoleId,
  getMoleculeType,
  type PolymerType,
} from './model/types.js';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Unit } from './structure.js';
import type { NumberArray } from '@molstar/core/util/type-helpers';

export function getCoarseBegCompId(unit: Unit.Spheres | Unit.Gaussians, element: ElementIndex) {
  const entityKey = unit.coarseElements.entityKey[element];
  const seq = unit.model.sequence.byEntityKey[entityKey].sequence;
  const seq_id_begin = unit.coarseElements.seq_id_begin.value(element);
  return seq.compId.value(seq_id_begin - 1); // 1-indexed
}

export function getElementMoleculeType(unit: Unit, element: ElementIndex): MoleculeType {
  switch (unit.kind) {
    case Unit.Kind.Atomic:
      return unit.model.atomicHierarchy.derived.residue.moleculeType[unit.residueIndex[element]];
    case Unit.Kind.Spheres:
    case Unit.Kind.Gaussians:
      // TODO add unit.model.coarseHierarchy.derived.residue.moleculeType
      const compId = getCoarseBegCompId(unit, element);
      const cc = unit.model.properties.chemicalComponentMap.get(compId);
      if (cc) return getMoleculeType(cc.type, compId);
  }
  return MoleculeType.Unknown;
}

export function getAtomicMoleculeType(model: Model, rI: ResidueIndex): MoleculeType {
  return model.atomicHierarchy.derived.residue.moleculeType[rI];
}

const EmptyAtomIds = new Set<string>();
export function getAtomIdForAtomRole(polymerType: PolymerType, atomRole: AtomRole) {
  const p = PolymerTypeAtomRoleId[polymerType];
  if (p !== undefined) {
    const a = p[atomRole];
    if (a !== undefined) return a;
  }
  return EmptyAtomIds;
}

const tmpPositionsVec = Vec3.zero();
export function getPositions(unit: Unit, indices: ArrayLike<number>): NumberArray {
  const c = unit.conformation;
  const positions = new Float32Array(indices.length * 3);
  const { elements } = unit;
  for (let i = 0, il = indices.length; i < il; ++i) {
    c.position(elements[indices[i]], tmpPositionsVec);
    Vec3.toArray(tmpPositionsVec, positions, i * 3);
  }
  return positions;
}
