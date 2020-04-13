/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, Structure } from '../../../../mol-model/structure';
import { Loci, EmptyLoci } from '../../../../mol-model/loci';
import { Interval } from '../../../../mol-data/int';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { PickingId } from '../../../../mol-geo/geometry/picking';
import { StructureGroup } from '../../../../mol-repr/structure/units-visual';
import { getResidueLoci } from './common';
import { eachAtomicUnitTracedElement } from './polymer';

export namespace NucleotideLocationIterator {
    export function fromGroup(structureGroup: StructureGroup): LocationIterator {
        const { group, structure } = structureGroup;
        const u = group.units[0];
        const nucleotideElementIndices = Unit.isAtomic(u) ? u.nucleotideElements : [];
        const groupCount = nucleotideElementIndices.length;
        const instanceCount = group.units.length;
        const location = StructureElement.Location.create(structure);
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex];
            location.unit = unit;
            location.element = nucleotideElementIndices[groupIndex];
            return location;
        };
        return LocationIterator(groupCount, instanceCount, getLocation);
    }
}

export function getNucleotideElementLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId;
    if (id === objectId) {
        const { structure, group } = structureGroup;
        const unit = group.units[instanceId];
        if (Unit.isAtomic(unit)) {
            return getResidueLoci(structure, unit, unit.nucleotideElements[groupId]);
        }
    }
    return EmptyLoci;
}

function selectNuclotideElements(u: Unit.Atomic) { return u.nucleotideElements; }

/**
 * Mark a nucleotide element (e.g. part of a cartoon block)
 * - mark only when all its residue's elements are in a loci
 */
export function eachNucleotideElement(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    const { structure, group } = structureGroup;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;
    const unit = group.units[0];
    if (!Unit.isAtomic(unit)) return false;
    const { nucleotideElements } = unit;
    const groupCount = nucleotideElements.length;
    for (const e of loci.elements) {
        if (!Unit.isAtomic(e.unit)) continue;
        if (!group.unitIndexMap.has(e.unit.id)) continue;

        const intervalOffset = group.unitIndexMap.get(e.unit.id) * groupCount;

        if (Unit.isAtomic(e.unit)) {
            changed = eachAtomicUnitTracedElement(intervalOffset, groupCount, selectNuclotideElements, apply, e) || changed;
        }
    }
    return changed;
}