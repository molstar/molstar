/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Structure } from '../structure.js';
import { Unit } from '../unit.js';
import type { ElementIndex } from '../../model/indexing.js';
import { StructureElement } from '../element.js';
import { type Coordination, type CoordinationIndex, EmptyCoordination } from './data.js';
import { cantorPairing } from '@molstar/core/data/util';
import { BondType } from '../../model/types.js';
import { AtomicNumbers } from '../../model/properties/atomic/measures.js';

/** Minimum number of bonds for an atom to be considered a coordination site */
const MinCoordination = 4;

function interBondCount(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex): number {
  let count = 0;
  const indices = structure.interUnitBonds.getEdgeIndices(index, unit.id);
  for (let i = 0, il = indices.length; i < il; ++i) {
    const b = structure.interUnitBonds.edges[indices[i]];
    if (BondType.isCovalent(b.props.flag) || BondType.is(b.props.flag, BondType.Flag.MetallicCoordination)) count += 1;
  }
  return count;
}

function intraBondCount(unit: Unit.Atomic, index: StructureElement.UnitIndex): number {
  let count = 0;
  const {
    offset,
    edgeProps: { flags },
  } = unit.bonds;
  for (let i = offset[index], il = offset[index + 1]; i < il; ++i) {
    if (BondType.isCovalent(flags[i]) || BondType.is(flags[i], BondType.Flag.MetallicCoordination)) count += 1;
  }
  return count;
}

function bondCount(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex): number {
  return interBondCount(structure, unit, index) + intraBondCount(unit, index);
}

//

const _H = AtomicNumbers['H'];

function eachInterBondedAtom(
  structure: Structure,
  unit: Unit.Atomic,
  index: StructureElement.UnitIndex,
  cb: (unit: Unit.Atomic, index: StructureElement.UnitIndex) => void,
): void {
  const indices = structure.interUnitBonds.getEdgeIndices(index, unit.id);
  for (let i = 0, il = indices.length; i < il; ++i) {
    const b = structure.interUnitBonds.edges[indices[i]];
    const uB = structure.unitMap.get(b.unitB) as Unit.Atomic;
    if (
      (BondType.isCovalent(b.props.flag) || BondType.is(b.props.flag, BondType.Flag.MetallicCoordination)) &&
      uB.model.atomicHierarchy.derived.atom.atomicNumber[uB.elements[b.indexB]] !== _H
    )
      cb(uB, b.indexB);
  }
}

function eachIntraBondedAtom(
  unit: Unit.Atomic,
  index: StructureElement.UnitIndex,
  cb: (unit: Unit.Atomic, index: StructureElement.UnitIndex) => void,
): void {
  const { elements } = unit;
  const { atomicNumber } = unit.model.atomicHierarchy.derived.atom;
  const {
    offset,
    b,
    edgeProps: { flags },
  } = unit.bonds;
  for (let i = offset[index], il = offset[index + 1]; i < il; ++i) {
    if (
      (BondType.isCovalent(flags[i]) || BondType.is(flags[i], BondType.Flag.MetallicCoordination)) &&
      atomicNumber[elements[b[i]]] !== _H
    )
      cb(unit, b[i] as StructureElement.UnitIndex);
  }
}

function eachBondedAtom(
  structure: Structure,
  unit: Unit.Atomic,
  index: StructureElement.UnitIndex,
  cb: (unit: Unit.Atomic, index: StructureElement.UnitIndex) => void,
): void {
  eachInterBondedAtom(structure, unit, index, cb);
  eachIntraBondedAtom(unit, index, cb);
}

//

export function computeCoordination(structure: Structure): Coordination {
  const unitIds: number[] = [];
  const indices: StructureElement.UnitIndex[] = [];
  const numbers: number[] = [];
  const siteIndex = new Map<number, CoordinationIndex>();

  for (let ui = 0, uil = structure.units.length; ui < uil; ++ui) {
    const unit = structure.units[ui];
    if (!Unit.isAtomic(unit)) continue;

    const { elements } = unit;
    for (let ei = 0, eil = elements.length; ei < eil; ++ei) {
      const element = elements[ei];
      const unitIndex = ei as StructureElement.UnitIndex;
      const _bondCount = bondCount(structure, unit, unitIndex);
      if (_bondCount >= MinCoordination) {
        siteIndex.set(coordinationKey(unit.id, element), siteIndex.size as CoordinationIndex);
        unitIds.push(unit.id);
        indices.push(unitIndex);
        numbers.push(_bondCount);
      }
    }
  }

  if (siteIndex.size === 0) return EmptyCoordination;

  const l = StructureElement.Location.create(structure);

  return {
    sites: {
      unitIds,
      indices,
      numbers,
      count: siteIndex.size,
    },
    getSiteIndex: (unit: Unit.Atomic, element: ElementIndex) => {
      return siteIndex.get(coordinationKey(unit.id, element)) ?? -1;
    },
    eachLigand: (siteIndex: CoordinationIndex, cb: (l: StructureElement.Location) => void) => {
      const unitId = unitIds[siteIndex];
      const element = indices[siteIndex];
      const unit = structure.unitMap.get(unitId) as Unit.Atomic;
      eachBondedAtom(structure, unit, element, (u, i) => {
        l.unit = u;
        l.element = u.elements[i];
        cb(l);
      });
    },
  };
}

function coordinationKey(unitId: number, element: ElementIndex) {
  return cantorPairing(unitId, element);
}
