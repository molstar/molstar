/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, ElementIndex, StructureElement } from 'mol-model/structure';
import SortedRanges from 'mol-data/int/sorted-ranges';
import { LocationIterator } from '../../../../util/location-iterator';
import { PickingId } from '../../../../util/picking';
import { OrderedSet, Interval } from 'mol-data/int';
import { EmptyLoci, Loci } from 'mol-model/loci';

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
        const unitIndex = OrderedSet.findPredecessorIndex(unit.elements, unit.polymerElements[groupId]) as StructureElement.UnitIndex
        const indices = OrderedSet.ofSingleton(unitIndex);
        return StructureElement.Loci([{ unit, indices }])
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
                    const start = unitIdx * groupCount + OrderedSet.findPredecessorIndex(e.unit.polymerElements, e.unit.elements[Interval.start(e.indices)])
                    const end = unitIdx * groupCount + OrderedSet.findPredecessorIndex(e.unit.polymerElements, e.unit.elements[Interval.end(e.indices)])
                    if (apply(Interval.ofBounds(start, end))) changed = true
                } else {
                    for (let i = 0, _i = e.indices.length; i < _i; i++) {
                        const idx = unitIdx * groupCount + OrderedSet.findPredecessorIndex(e.unit.polymerElements, e.unit.elements[e.indices[i]])
                        if (apply(Interval.ofSingleton(idx))) changed = true
                    }
                }
            }
        }
    }
    return changed
}