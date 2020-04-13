/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, ElementIndex, StructureElement, Bond, Structure, ResidueIndex } from '../../../../mol-model/structure';
import SortedRanges from '../../../../mol-data/int/sorted-ranges';
import { OrderedSet, Interval, SortedArray } from '../../../../mol-data/int';
import { EmptyLoci, Loci } from '../../../../mol-model/loci';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { PickingId } from '../../../../mol-geo/geometry/picking';
import { StructureGroup } from '../../../structure/units-visual';
import { getResidueLoci } from './common';

export * from './polymer/backbone-iterator';
export * from './polymer/gap-iterator';
export * from './polymer/trace-iterator';
export * from './polymer/curve-segment';

export const StandardTension = 0.5;
export const HelixTension = 0.9;
export const StandardShift = 0.5;
export const NucleicShift = 0.3;
export const OverhangFactor = 2;

export function getPolymerRanges(unit: Unit): SortedRanges<ElementIndex> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return unit.model.atomicRanges.polymerRanges;
        case Unit.Kind.Spheres: return unit.model.coarseHierarchy.spheres.polymerRanges;
        case Unit.Kind.Gaussians: return unit.model.coarseHierarchy.gaussians.polymerRanges;
    }
}

export function getGapRanges(unit: Unit): SortedRanges<ElementIndex> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return unit.model.atomicRanges.gapRanges;
        case Unit.Kind.Spheres: return unit.model.coarseHierarchy.spheres.gapRanges;
        case Unit.Kind.Gaussians: return unit.model.coarseHierarchy.gaussians.gapRanges;
    }
}

export namespace PolymerLocationIterator {
    export function fromGroup(structureGroup: StructureGroup): LocationIterator {
        const { group, structure } = structureGroup;
        const polymerElements = group.units[0].polymerElements;
        const groupCount = polymerElements.length;
        const instanceCount = group.units.length;
        const location = StructureElement.Location.create(structure);
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex];
            location.unit = unit;
            location.element = polymerElements[groupIndex];
            return location;
        };
        return LocationIterator(groupCount, instanceCount, getLocation);
    }
}

export namespace PolymerGapLocationIterator {
    export function fromGroup(structureGroup: StructureGroup): LocationIterator {
        const { group, structure } = structureGroup;
        const gapElements = group.units[0].gapElements;
        const groupCount = gapElements.length;
        const instanceCount = group.units.length;
        const location = StructureElement.Location.create(structure);
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex];
            location.unit = unit;
            location.element = gapElements[groupIndex];
            return location;
        };
        return LocationIterator(groupCount, instanceCount, getLocation);
    }
}

/** Return a Loci for the elements of the whole residue of a polymer element. */
export function getPolymerElementLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId;
    if (id === objectId) {
        const { structure, group } = structureGroup;
        const unit = group.units[instanceId];
        if (Unit.isAtomic(unit)) {
            return getResidueLoci(structure, unit, unit.polymerElements[groupId]);
        } else {
            const { elements } = unit;
            const elementIndex = unit.polymerElements[groupId];
            const unitIndex = OrderedSet.indexOf(elements, elementIndex) as StructureElement.UnitIndex | -1;
            if (unitIndex !== -1) {
                const indices = OrderedSet.ofSingleton(unitIndex);
                return StructureElement.Loci(structure, [{ unit, indices }]);
            }
        }
    }
    return EmptyLoci;
}


function tryApplyResidueInterval(offset: number, elements: SortedArray<ElementIndex>, traceElementIndex: ArrayLike<ElementIndex | -1>, apply: (interval: Interval) => boolean, r1: ResidueIndex, r2: ResidueIndex) {
    let start = -1, startIdx = -1;

    for (let rI = r1; rI <= r2; rI++) {
        const eI = traceElementIndex[rI];
        if (eI < 0) continue;
        start = OrderedSet.indexOf(elements, eI);
        if (start >= 0) {
            startIdx = rI;
            break;
        }
    }

    if (start < 0) {
        return false;
    }

    let end = start;

    for (let rI = r2; rI > startIdx; rI--) {
        const eI = traceElementIndex[rI];
        if (eI < 0) continue;
        const e = OrderedSet.indexOf(elements, eI);
        if (e >= 0) {
            end = e;
            break;
        }
    }

    return apply(Interval.ofRange(offset + start, offset + end));
}

export function eachAtomicUnitTracedElement(offset: number, groupSize: number, elementsSelector: (u: Unit.Atomic) => SortedArray<ElementIndex>, apply: (interval: Interval) => boolean, e: StructureElement.Loci['elements'][0]) {
    let changed = false;

    const { elements } = e.unit;
    const { traceElementIndex } = e.unit.model.atomicHierarchy.derived.residue;
    const { index: resIndex } = e.unit.model.atomicHierarchy.residueAtomSegments;
    const tracedElements = elementsSelector(e.unit as Unit.Atomic);

    if (Interval.is(e.indices)) {
        if (Interval.start(e.indices) === 0 && Interval.end(e.indices) === e.unit.elements.length) {
            // full unit here
            changed = apply(Interval.ofBounds(offset, offset + groupSize)) || changed;
        } else {
            let r1 = resIndex[elements[Interval.min(e.indices)]];
            let r2 = resIndex[elements[Interval.max(e.indices)]];
            changed = tryApplyResidueInterval(offset, tracedElements, traceElementIndex, apply, r1, r2) || changed;
        }
    } else {
        const { indices } = e;

        for (let i = 0, _i = indices.length; i < _i; i++) {
            const r1 = resIndex[elements[indices[i]]];
            let r2 = r1;

            let endI = i + 1;
            while (endI < _i) {
                const _r = resIndex[elements[indices[endI]]];
                if (_r - r2 > 1) break;
                r2 = _r;
                endI++;
            }
            i = endI - 1;
            changed = tryApplyResidueInterval(offset, tracedElements, traceElementIndex, apply, r1, r2) || changed;
        }
    }

    return changed;
}

function selectPolymerElements(u: Unit) { return u.polymerElements; }

/** Mark a polymer element (e.g. part of a cartoon trace) */
export function eachPolymerElement(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    const { structure, group } = structureGroup;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;
    const groupCount = group.units[0].polymerElements.length;
    for (const e of loci.elements) {
        if (!group.unitIndexMap.has(e.unit.id)) continue;

        const offset = group.unitIndexMap.get(e.unit.id) * groupCount; // to target unit instance

        if (Unit.isAtomic(e.unit)) {
            changed = eachAtomicUnitTracedElement(offset, groupCount, selectPolymerElements, apply, e) || changed;
        } else {
            if (Interval.is(e.indices)) {
                const start = offset + Interval.start(e.indices);
                const end = offset + Interval.end(e.indices);
                changed = apply(Interval.ofBounds(start, end)) || changed;
            } else {
                for (let i = 0, _i = e.indices.length; i < _i; i++) {
                    const start = e.indices[i];
                    let endI = i + 1;
                    while (endI < _i && e.indices[endI] === start) endI++;
                    i = endI - 1;
                    const end = e.indices[i];
                    changed = apply(Interval.ofRange(offset + start, offset + end)) || changed;
                }
            }
        }
    }
    return changed;
}

/** Return a Loci for both directions of the polymer gap element. */
export function getPolymerGapElementLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId;
    if (id === objectId) {
        const { structure, group } = structureGroup;
        const unit = group.units[instanceId];
        const unitIndexA = OrderedSet.indexOf(unit.elements, unit.gapElements[groupId]) as StructureElement.UnitIndex;
        const unitIndexB = OrderedSet.indexOf(unit.elements, unit.gapElements[groupId % 2 ? groupId - 1 : groupId + 1]) as StructureElement.UnitIndex;
        if (unitIndexA !== -1 && unitIndexB !== -1) {
            return Bond.Loci(structure, [
                Bond.Location(structure, unit, unitIndexA, structure, unit, unitIndexB),
                Bond.Location(structure, unit, unitIndexB, structure, unit, unitIndexA)
            ]);
        }
    }
    return EmptyLoci;
}

export function eachPolymerGapElement(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (Bond.isLoci(loci)) {
        const { structure, group } = structureGroup;
        if (!Structure.areRootsEquivalent(loci.structure, structure)) return false;
        loci = Bond.remapLoci(loci, structure);
        const groupCount = group.units[0].gapElements.length;
        for (const b of loci.bonds) {
            const unitIdx = group.unitIndexMap.get(b.aUnit.id);
            if (unitIdx !== undefined) {
                const idxA = OrderedSet.indexOf(b.aUnit.gapElements, b.aUnit.elements[b.aIndex]);
                const idxB = OrderedSet.indexOf(b.bUnit.gapElements, b.bUnit.elements[b.bIndex]);
                if (idxA !== -1 && idxB !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idxA))) changed = true;
                }
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        const { structure, group } = structureGroup;
        if (!Structure.areRootsEquivalent(loci.structure, structure)) return false;
        loci = StructureElement.Loci.remap(loci, structure);
        const groupCount = group.units[0].gapElements.length;
        for (const e of loci.elements) {
            const unitIdx = group.unitIndexMap.get(e.unit.id);
            if (unitIdx !== undefined) {
                OrderedSet.forEach(e.indices, v => {
                    const idx = OrderedSet.indexOf(e.unit.gapElements, e.unit.elements[v]);
                    if (idx !== -1) {
                        if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true;
                    }
                });
            }
        }
    }
    return changed;
}