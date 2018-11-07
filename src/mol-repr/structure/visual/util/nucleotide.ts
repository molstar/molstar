/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement } from 'mol-model/structure';
import { getNucleotideElements } from 'mol-model/structure/structure/util/nucleotide';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { OrderedSet, Interval } from 'mol-data/int';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { PickingId } from 'mol-geo/geometry/picking';
import { StructureGroup } from 'mol-repr/structure/units-visual';

export namespace NucleotideLocationIterator {
    export function fromGroup(group: Unit.SymmetryGroup): LocationIterator {
        const u = group.units[0]
        const nucleotideElementIndices = Unit.isAtomic(u) ? getNucleotideElements(u) : []
        const groupCount = nucleotideElementIndices.length
        const instanceCount = group.units.length
        const location = StructureElement.create()
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex]
            location.unit = unit
            location.element = nucleotideElementIndices[groupIndex]
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }
}

export function getNucleotideElementLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const { structure, group } = structureGroup
        const unit = group.units[instanceId]
        if (Unit.isAtomic(unit)) {
            const unitIndex = OrderedSet.indexOf(unit.elements, unit.nucleotideElements[groupId]) as StructureElement.UnitIndex
            if (unitIndex !== -1) {
                const indices = OrderedSet.ofSingleton(unitIndex)
                return StructureElement.Loci(structure, [{ unit, indices }])
            }
        }
    }
    return EmptyLoci
}

export function markNucleotideElement(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    if (!StructureElement.isLoci(loci)) return false
    const { structure, group } = structureGroup
    if (loci.structure !== structure) return false
    const unit = group.units[0]
    if (!Unit.isAtomic(unit)) return false
    const groupCount = unit.nucleotideElements.length
    for (const e of loci.elements) {
        const unitIdx = group.unitIndexMap.get(e.unit.id)
        if (unitIdx !== undefined && Unit.isAtomic(e.unit)) {
            if (Interval.is(e.indices)) {
                const min = OrderedSet.indexOf(e.unit.nucleotideElements, e.unit.elements[Interval.min(e.indices)])
                const max = OrderedSet.indexOf(e.unit.nucleotideElements, e.unit.elements[Interval.max(e.indices)])
                if (min !== -1 && max !== -1) {
                    if (apply(Interval.ofRange(unitIdx * groupCount + min, unitIdx * groupCount + max))) changed = true
                }
            } else {
                for (let i = 0, _i = e.indices.length; i < _i; i++) {
                    const idx = OrderedSet.indexOf(e.unit.nucleotideElements, e.unit.elements[e.indices[i]])
                    if (idx !== -1) {
                        if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true
                    }
                }
            }
        }
    }
    return changed
}