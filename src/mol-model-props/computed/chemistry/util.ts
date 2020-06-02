/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Unit } from '../../../mol-model/structure';
import { StructureElement } from '../../../mol-model/structure/structure';
import { Elements } from '../../../mol-model/structure/model/properties/atomic/types';
import { BondType } from '../../../mol-model/structure/model/types';
import { SortedArray } from '../../../mol-data/int';

export function typeSymbol(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return unit.model.atomicHierarchy.atoms.type_symbol.value(unit.elements[index]);
}

export function formalCharge(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return unit.model.atomicHierarchy.atoms.pdbx_formal_charge.value(unit.elements[index]);
}

export function atomId(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return unit.model.atomicHierarchy.atoms.label_atom_id.value(unit.elements[index]);
}

export function altLoc(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return unit.model.atomicHierarchy.atoms.label_alt_id.value(unit.elements[index]);
}

export function compId(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return unit.model.atomicHierarchy.atoms.label_comp_id.value(unit.elements[index]);
}

//

export function interBondCount(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex): number {
    let count = 0;
    const indices = structure.interUnitBonds.getEdgeIndices(index, unit);
    for (let i = 0, il = indices.length; i < il; ++i) {
        const b = structure.interUnitBonds.edges[indices[i]];
        if (BondType.isCovalent(b.props.flag)) count += 1;
    }
    return count;
}

export function intraBondCount(unit: Unit.Atomic, index: StructureElement.UnitIndex): number {
    let count = 0;
    const { offset, edgeProps: { flags } } = unit.bonds;
    for (let i = offset[index], il = offset[index + 1]; i < il; ++i) {
        if (BondType.isCovalent(flags[i])) count += 1;
    }
    return count;
}

export function bondCount(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex): number {
    return interBondCount(structure, unit, index) + intraBondCount(unit, index);
}

export function bondToElementCount(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex, element: Elements): number {
    let count = 0;
    eachBondedAtom(structure, unit, index, (unit: Unit.Atomic, index: StructureElement.UnitIndex) => {
        if (typeSymbol(unit, index) === element) count += 1;
    });
    return count;
}

//

export function intraConnectedTo(unit: Unit.Atomic, indexA: StructureElement.UnitIndex, indexB: StructureElement.UnitIndex) {
    const { offset, b, edgeProps: { flags } } = unit.bonds;
    BondType.is;
    for (let i = offset[indexA], il = offset[indexA + 1]; i < il; ++i) {
        if (b[i] === indexB && BondType.isCovalent(flags[i])) return true;
    }
    return false;
}

export function interConnectedTo(structure: Structure, unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex) {
    const b = structure.interUnitBonds.getEdge(indexA, unitA, indexB, unitB);
    return b && BondType.isCovalent(b.props.flag);
}

export function connectedTo(structure: Structure, unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex) {
    return unitA === unitB ? intraConnectedTo(unitA, indexA, indexB) : interConnectedTo(structure, unitA, indexA, unitB, indexB);
}

//

export function eachInterBondedAtom(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex, cb: (unit: Unit.Atomic, index: StructureElement.UnitIndex) => void): void {
    const indices = structure.interUnitBonds.getEdgeIndices(index, unit);
    for (let i = 0, il = indices.length; i < il; ++i) {
        const b = structure.interUnitBonds.edges[indices[i]];
        if (BondType.isCovalent(b.props.flag)) cb(b.unitB, b.indexB);
    }
}

export function eachIntraBondedAtom(unit: Unit.Atomic, index: StructureElement.UnitIndex, cb: (unit: Unit.Atomic, index: StructureElement.UnitIndex) => void): void {
    const { offset, b, edgeProps: { flags } } = unit.bonds;
    for (let i = offset[index], il = offset[index + 1]; i < il; ++i) {
        if (BondType.isCovalent(flags[i])) cb(unit, b[i] as StructureElement.UnitIndex);
    }
}

export function eachBondedAtom(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex, cb: (unit: Unit.Atomic, index: StructureElement.UnitIndex) => void): void {
    eachInterBondedAtom(structure, unit, index, cb);
    eachIntraBondedAtom(unit, index, cb);
}

//

export function eachResidueAtom(unit: Unit.Atomic, index: StructureElement.UnitIndex, cb: (index: StructureElement.UnitIndex) => void): void {
    const { offsets } = unit.model.atomicHierarchy.residueAtomSegments;
    const rI = unit.getResidueIndex(index);
    for (let i = offsets[rI], il = offsets[rI + 1]; i < il; ++i) {
        // TODO optimize, avoid search with .indexOf
        const idx = SortedArray.indexOf(unit.elements, i);
        if (idx !== -1) cb(idx as StructureElement.UnitIndex);
    }
}