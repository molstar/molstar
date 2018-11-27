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
import { StructureGroup } from 'mol-repr/structure/units-visual';
import { getResidueLoci } from './common';

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

/** Return a Loci for the elements of the whole residue of a polymer element. */
export function getPolymerElementLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const { structure, group } = structureGroup
        const unit = group.units[instanceId]
        return getResidueLoci(structure, unit, unit.polymerElements[groupId])
    }
    return EmptyLoci
}

/** Mark a polymer element (e.g. part of a cartoon trace) when all its residue's elements are in a loci. */
export function markPolymerElement(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    if (!StructureElement.isLoci(loci)) return false
    const { structure, group } = structureGroup
    if (loci.structure !== structure) return false
    const { polymerElements, model, elements } = group.units[0]
    const { index, offsets } = model.atomicHierarchy.residueAtomSegments
    const { traceElementIndex } = model.atomicHierarchy.derived.residue
    const groupCount = polymerElements.length
    for (const e of loci.elements) {
        const unitIdx = group.unitIndexMap.get(e.unit.id)
        if (unitIdx !== undefined) {
            // TODO optimized implementation for intervals
            OrderedSet.forEach(e.indices, v => {
                const rI = index[elements[v]]
                const unitIndexMin = OrderedSet.findPredecessorIndex(elements, offsets[rI])
                const unitIndexMax = OrderedSet.findPredecessorIndex(elements, offsets[rI + 1] - 1)
                const unitIndexInterval = Interval.ofRange(unitIndexMin, unitIndexMax)
                if (!OrderedSet.isSubset(e.indices, unitIndexInterval)) return
                const eI = traceElementIndex[rI]
                const idx = OrderedSet.indexOf(e.unit.polymerElements, eI)
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true
                }
            })
        }
    }
    return changed
}

/** Return a Loci for both directions of the polymer gap element. */
export function getPolymerGapElementLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const { structure, group } = structureGroup
        const unit = group.units[instanceId]
        const unitIndexA = OrderedSet.indexOf(unit.elements, unit.gapElements[groupId]) as StructureElement.UnitIndex
        const unitIndexB = OrderedSet.indexOf(unit.elements, unit.gapElements[groupId % 2 ? groupId - 1 : groupId + 1]) as StructureElement.UnitIndex
        if (unitIndexA !== -1 && unitIndexB !== -1) {
            return Link.Loci(structure, [
                Link.Location(unit, unitIndexA, unit, unitIndexB),
                Link.Location(unit, unitIndexB, unit, unitIndexA)
            ])
        }
    }
    return EmptyLoci
}

export function markPolymerGapElement(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    if (!Link.isLoci(loci)) return false
    const { structure, group } = structureGroup
    if (loci.structure !== structure) return false
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
    return changed
}