/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, ElementIndex, StructureElement } from 'mol-model/structure';
import { Mat4 } from 'mol-math/linear-algebra';
import { TransformData, createTransform } from 'mol-geo/geometry/transform-data';
import { OrderedSet, SortedArray } from 'mol-data/int';
import { EmptyLoci, Loci } from 'mol-model/loci';

/** Return a Loci for the elements of a whole residue the elementIndex belongs to. */
export function getResidueLoci(structure: Structure, unit: Unit.Atomic, elementIndex: ElementIndex): Loci {
    const { elements, model } = unit
    if (OrderedSet.indexOf(elements, elementIndex) !== -1) {
        const { index, offsets } = model.atomicHierarchy.residueAtomSegments
        const rI = index[elementIndex]
        const _indices: number[] = []
        for (let i = offsets[rI], il = offsets[rI + 1]; i < il; ++i) {
            const unitIndex = OrderedSet.indexOf(elements, i)
            if (unitIndex !== -1) _indices.push(unitIndex)
        }
        const indices = OrderedSet.ofSortedArray<StructureElement.UnitIndex>(SortedArray.ofSortedArray(_indices))
        return StructureElement.Loci(structure, [{ unit, indices }])
    }
    return EmptyLoci
}

//

export function createUnitsTransform({ units }: Unit.SymmetryGroup, transformData?: TransformData) {
    const unitCount = units.length
    const n = unitCount * 16
    const array = transformData && transformData.aTransform.ref.value.length >= n ? transformData.aTransform.ref.value : new Float32Array(n)
    for (let i = 0; i < unitCount; i++) {
        Mat4.toArray(units[i].conformation.operator.matrix, array, i * 16)
    }
    return createTransform(array, unitCount, transformData)
}

export const UnitKindInfo = {
    'atomic': {},
    'spheres': {},
    'gaussians': {},
}
export type UnitKind = keyof typeof UnitKindInfo
export const UnitKindNames = Object.keys(UnitKindInfo)
export const UnitKindOptions = UnitKindNames.map(n => [n, n] as [UnitKind, string])

export function includesUnitKind(unitKinds: UnitKind[], unit: Unit) {
    for (let i = 0, il = unitKinds.length; i < il; ++i) {
        if (Unit.isAtomic(unit) && unitKinds[i] === 'atomic') return true
        if (Unit.isSpheres(unit) && unitKinds[i] === 'spheres') return true
        if (Unit.isGaussians(unit) && unitKinds[i] === 'gaussians') return true
    }
    return false
}