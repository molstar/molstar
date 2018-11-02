/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, ElementIndex, StructureElement, Link } from 'mol-model/structure';
import SortedRanges from 'mol-data/int/sorted-ranges';
import { OrderedSet, Interval } from 'mol-data/int';
import { EmptyLoci, Loci } from 'mol-model/loci';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';

export * from './polymer/backbone-iterator'
export * from './polymer/gap-iterator'
export * from './polymer/trace-iterator'
export * from './polymer/curve-segment'

export function getPolymerRanges(unit: Unit): SortedRanges<ElementIndex> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return unit.model.atomicHierarchy.polymerRanges
        case Unit.Kind.Spheres: return unit.model.coarseHierarchy.spheres.polymerRanges
        case Unit.Kind.Gaussians: return unit.model.coarseHierarchy.gaussians.polymerRanges
    }
}

export function getGapRanges(unit: Unit): SortedRanges<ElementIndex> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return unit.model.atomicHierarchy.gapRanges
        case Unit.Kind.Spheres: return unit.model.coarseHierarchy.spheres.gapRanges
        case Unit.Kind.Gaussians: return unit.model.coarseHierarchy.gaussians.gapRanges
    }
}

export namespace PolymerLocationIterator {
    export function fromGroup(group: Unit.SymmetryGroup): LocationIterator {
        const polymerElements = group.units[0].polymerElements
        const groupCount = polymerElements.length
        const instanceCount = group.units.length
        const location = StructureElement.create()
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex]
            location.unit = unit
            location.element = polymerElements[groupIndex]
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }
}

export namespace PolymerGapLocationIterator {
    export function fromGroup(group: Unit.SymmetryGroup): LocationIterator {
        const gapElements = group.units[0].gapElements
        const groupCount = gapElements.length
        const instanceCount = group.units.length
        const location = StructureElement.create()
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex]
            location.unit = unit
            location.element = gapElements[groupIndex]
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }
}

export function getPolymerElementLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const unit = group.units[instanceId]
        if (unit === undefined) {
            console.log(id, { objectId, instanceId, groupId }, group.units)
        }
        const unitIndex = OrderedSet.indexOf(unit.elements, unit.polymerElements[groupId]) as StructureElement.UnitIndex
        if (unitIndex !== -1) {
            const indices = OrderedSet.ofSingleton(unitIndex)
            return StructureElement.Loci([{ unit, indices }])
        }
    }
    return EmptyLoci
}

export function markPolymerElement(loci: Loci, group: Unit.SymmetryGroup, apply: (interval: Interval) => boolean) {
    const groupCount = group.units[0].polymerElements.length

    let changed = false
    if (StructureElement.isLoci(loci)) {
        for (const e of loci.elements) {
            const unitIdx = group.unitIndexMap.get(e.unit.id)
            if (unitIdx !== undefined) {
                if (Interval.is(e.indices)) {
                    const min =  + OrderedSet.indexOf(e.unit.polymerElements, e.unit.elements[Interval.min(e.indices)])
                    const max = OrderedSet.indexOf(e.unit.polymerElements, e.unit.elements[Interval.max(e.indices)])
                    if (min !== -1 && max !== -1) {
                        if (apply(Interval.ofRange(unitIdx * groupCount + min, unitIdx * groupCount + max))) changed = true
                    }
                } else {
                    for (let i = 0, _i = e.indices.length; i < _i; i++) {
                        const idx = OrderedSet.indexOf(e.unit.polymerElements, e.unit.elements[e.indices[i]])
                        if (idx !== -1) {
                            if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true
                        }
                    }
                }
            }
        }
    }
    return changed
}

export function getPolymerGapElementLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const unit = group.units[instanceId]
        const unitIndexA = OrderedSet.indexOf(unit.elements, unit.gapElements[groupId]) as StructureElement.UnitIndex
        const unitIndexB = OrderedSet.indexOf(unit.elements, unit.gapElements[groupId % 2 ? groupId - 1 : groupId + 1]) as StructureElement.UnitIndex
        if (unitIndexA !== -1 && unitIndexB !== -1) {
            return Link.Loci([ Link.Location(unit, unitIndexA, unit, unitIndexB) ])
        }
    }
    return EmptyLoci
}

export function markPolymerGapElement(loci: Loci, group: Unit.SymmetryGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    if (Link.isLoci(loci)) {
        const groupCount = group.units[0].gapElements.length
        for (const b of loci.links) {
            const unitIdx = group.unitIndexMap.get(b.aUnit.id)
            if (unitIdx !== undefined) {
                const idxA = OrderedSet.indexOf(b.aUnit.gapElements, b.aUnit.elements[b.aIndex])
                const idxB = OrderedSet.indexOf(b.bUnit.gapElements, b.bUnit.elements[b.bIndex])
                if (idxA !== -1 && idxB !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idxA))) changed = true
                }
            }
        }
    }
    return changed
}