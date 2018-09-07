/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement } from 'mol-model/structure';
import { LocationIterator } from '../../../../util/location-iterator';
import { getNucleotideElements } from 'mol-model/structure/structure/util/nucleotide';
import { PickingId } from '../../../../util/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { OrderedSet, Interval } from 'mol-data/int';

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

export function getNucleotideElementLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const unit = group.units[instanceId]
        if (Unit.isAtomic(unit)) {
            const unitIndex = OrderedSet.indexOf(unit.elements, unit.nucleotideElements[groupId]) as StructureElement.UnitIndex
            if (unitIndex !== -1) {
                const indices = OrderedSet.ofSingleton(unitIndex)
                return StructureElement.Loci([{ unit, indices }])
            }
        }
    }
    return EmptyLoci
}

export function markNucleotideElement(loci: Loci, group: Unit.SymmetryGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    const u = group.units[0]
    if (StructureElement.isLoci(loci) && Unit.isAtomic(u)) {
        const groupCount = u.nucleotideElements.length
        for (const e of loci.elements) {
            const unitIdx = group.unitIndexMap.get(e.unit.id)
            if (unitIdx !== undefined && Unit.isAtomic(e.unit)) {
                if (Interval.is(e.indices)) {
                    const min = unitIdx * groupCount + OrderedSet.indexOf(e.unit.nucleotideElements, e.unit.elements[Interval.min(e.indices)])
                    const max = unitIdx * groupCount + OrderedSet.indexOf(e.unit.nucleotideElements, e.unit.elements[Interval.max(e.indices)])
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
    }
    return changed
}