/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, ElementIndex, StructureElement, ResidueIndex } from '../../../../mol-model/structure';
import { Mat4 } from '../../../../mol-math/linear-algebra';
import { TransformData, createTransform } from '../../../../mol-geo/geometry/transform-data';
import { OrderedSet, SortedArray } from '../../../../mol-data/int';
import { EmptyLoci, Loci } from '../../../../mol-model/loci';
import { PhysicalSizeTheme } from '../../../../mol-theme/size/physical';
import { AtomicNumbers, AtomNumber } from '../../../../mol-model/structure/model/properties/atomic';
import { fillSerial } from '../../../../mol-util/array';

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

/**
 * Return a Loci for the elements of a whole residue the elementIndex belongs to but
 * restrict to elements that have the same label_alt_id or none
 */
export function getAltResidueLoci(structure: Structure, unit: Unit.Atomic, elementIndex: ElementIndex) {
    const { elements, model } = unit
    const { label_alt_id } = model.atomicHierarchy.atoms
    const elementAltId = label_alt_id.value(elementIndex)
    if (OrderedSet.indexOf(elements, elementIndex) !== -1) {
        const { index } = model.atomicHierarchy.residueAtomSegments
        const rI = index[elementIndex]
        return getAltResidueLociFromId(structure, unit, rI, elementAltId)
    }
    return StructureElement.Loci(structure, [])
}

export function getAltResidueLociFromId(structure: Structure, unit: Unit.Atomic, residueIndex: ResidueIndex, elementAltId: string) {
    const { elements, model } = unit
    const { label_alt_id } = model.atomicHierarchy.atoms
    const { offsets } = model.atomicHierarchy.residueAtomSegments

    const _indices: number[] = []
    for (let i = offsets[residueIndex], il = offsets[residueIndex + 1]; i < il; ++i) {
        const unitIndex = OrderedSet.indexOf(elements, i)
        if (unitIndex !== -1) {
            const altId = label_alt_id.value(i)
            if (elementAltId === altId || altId === '') {
                _indices.push(unitIndex)
            }
        }
    }
    const indices = OrderedSet.ofSortedArray<StructureElement.UnitIndex>(SortedArray.ofSortedArray(_indices))
    return StructureElement.Loci(structure, [{ unit, indices }])
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

//

export function getConformation(unit: Unit) {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return unit.model.atomicConformation
        case Unit.Kind.Spheres: return unit.model.coarseConformation.spheres
        case Unit.Kind.Gaussians: return unit.model.coarseConformation.gaussians
    }
}

export function getUnitConformationAndRadius(structure: Structure, unit: Unit, ignoreHydrogens = false) {
    const conformation = getConformation(unit)
    const { elements } = unit

    let indices: SortedArray<ElementIndex>
    let id: ArrayLike<number>

    if (ignoreHydrogens) {
        const _indices = []
        const _id = []
        for (let i = 0, il = elements.length; i < il; ++i) {
            if (isHydrogen(unit, elements[i])) continue
            _indices.push(elements[i])
            _id.push(i)
        }
        indices = SortedArray.ofSortedArray(_indices)
        id = _id
    } else {
        indices = elements
        id = fillSerial(new Uint32Array(indices.length))
    }

    const position = {
        indices,
        x: conformation.x,
        y: conformation.y,
        z: conformation.z,
        id
    }

    const l = StructureElement.Location.create(structure, unit)
    const sizeTheme = PhysicalSizeTheme({}, {})
    const radius = (index: number) => {
        l.element = index as ElementIndex
        return sizeTheme.size(l)
    }

    return { position, radius }
}

export function getStructureConformationAndRadius(structure: Structure, ignoreHydrogens = false) {
    const l = StructureElement.Location.create(structure)
    const sizeTheme = PhysicalSizeTheme({}, {})

    let xs: ArrayLike<number>
    let ys: ArrayLike<number>
    let zs: ArrayLike<number>
    let rs: ArrayLike<number>
    let id: ArrayLike<number>

    if (ignoreHydrogens) {
        const _xs: number[] = []
        const _ys: number[] = []
        const _zs: number[] = []
        const _rs: number[] = []
        const _id: number[] = []
        for (let i = 0, m = 0, il = structure.units.length; i < il; ++i) {
            const unit = structure.units[i]
            const { elements } = unit
            const { x, y, z } = unit.conformation
            l.unit = unit
            for (let j = 0, jl = elements.length; j < jl; ++j) {
                const eI = elements[j]
                if (ignoreHydrogens && isHydrogen(unit, eI)) continue

                _xs.push(x(eI))
                _ys.push(y(eI))
                _zs.push(z(eI))
                l.element = eI
                _rs.push(sizeTheme.size(l))
                _id.push(m + j)
            }
            m += elements.length
        }
        xs = _xs, ys = _ys, zs = _zs, rs = _rs
        id = _id
    } else {
        const { elementCount } = structure
        const _xs = new Float32Array(elementCount)
        const _ys = new Float32Array(elementCount)
        const _zs = new Float32Array(elementCount)
        const _rs = new Float32Array(elementCount)
        for (let i = 0, m = 0, il = structure.units.length; i < il; ++i) {
            const unit = structure.units[i]
            const { elements } = unit
            const { x, y, z } = unit.conformation
            l.unit = unit
            for (let j = 0, jl = elements.length; j < jl; ++j) {
                const eI = elements[j]

                const mj = m + j
                _xs[mj] = x(eI)
                _ys[mj] = y(eI)
                _zs[mj] = z(eI)
                l.element = eI
                _rs[mj] = sizeTheme.size(l)
            }
            m += elements.length
        }
        xs = _xs, ys = _ys, zs = _zs, rs = _rs
        id = fillSerial(new Uint32Array(elementCount))
    }

    const position = { indices: OrderedSet.ofRange(0, id.length), x: xs, y: ys, z: zs, id }
    const radius = (index: number) => rs[index]

    return { position, radius }
}

const _H = AtomicNumbers['H']
export function isHydrogen(unit: Unit, element: ElementIndex) {
    if (Unit.isCoarse(unit)) return false
    if (AtomNumber(unit.model.atomicHierarchy.atoms.type_symbol.value(element)) === _H) return true
    return false
}