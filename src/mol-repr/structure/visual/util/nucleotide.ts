/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, Structure, ResidueIndex, ElementIndex } from '../../../../mol-model/structure';
import { Loci, EmptyLoci } from '../../../../mol-model/loci';
import { Interval } from '../../../../mol-data/int';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { PickingId } from '../../../../mol-geo/geometry/picking';
import { getResidueLoci, StructureGroup } from './common';
import { eachAtomicUnitTracedElement } from './polymer';
import { isPurineBase, isPyrimidineBase } from '../../../../mol-model/structure/model/types';
import { Vec3 } from '../../../../mol-math/linear-algebra/3d/vec3';

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
        return LocationIterator(groupCount, instanceCount, 1, getLocation);
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

//

const pC4 = Vec3();
const pN9 = Vec3();

export function getNucleotideBaseType(unit: Unit.Atomic, residueIndex: ResidueIndex) {
    const { model, conformation: c } = unit;
    const { residueAtomSegments, atoms, index: atomicIndex } = model.atomicHierarchy;
    const { label_comp_id } = atoms;

    const compId = label_comp_id.value(residueAtomSegments.offsets[residueIndex]);

    let isPurine = isPurineBase(compId);
    let isPyrimidine = isPyrimidineBase(compId);

    if (!isPurine && !isPyrimidine) {
        // detect Purine or Pyrimidin based on geometry
        const idxC4 = atomicIndex.findAtomOnResidue(residueIndex, 'C4');
        const idxN9 = atomicIndex.findAtomOnResidue(residueIndex, 'N9');
        if (idxC4 !== -1 && idxN9 !== -1 && Vec3.distance(c.invariantPosition(idxC4, pC4), c.invariantPosition(idxN9, pN9)) < 1.6) {
            isPurine = true;
        } else {
            isPyrimidine = true;
        }
    }

    return { isPurine, isPyrimidine };
}

export function createNucleicIndices() {
    return {
        trace: -1 as ElementIndex | -1,
        N1: -1 as ElementIndex | -1,
        C2: -1 as ElementIndex | -1,
        N3: -1 as ElementIndex | -1,
        C4: -1 as ElementIndex | -1,
        C5: -1 as ElementIndex | -1,
        C6: -1 as ElementIndex | -1,
        N7: -1 as ElementIndex | -1,
        C8: -1 as ElementIndex | -1,
        N9: -1 as ElementIndex | -1,
        C1_1: -1 as ElementIndex | -1,
        C2_1: -1 as ElementIndex | -1,
        C3_1: -1 as ElementIndex | -1,
        C4_1: -1 as ElementIndex | -1,
        O4_1: -1 as ElementIndex | -1,
    };
}
export type NucleicIndices = ReturnType<typeof createNucleicIndices>

export function setPurinIndices(idx: NucleicIndices, unit: Unit.Atomic, residueIndex: ResidueIndex) {
    const atomicIndex = unit.model.atomicHierarchy.index;
    const { traceElementIndex } = unit.model.atomicHierarchy.derived.residue;

    idx.trace = traceElementIndex[residueIndex];
    idx.N1 = atomicIndex.findAtomOnResidue(residueIndex, 'N1');
    idx.C2 = atomicIndex.findAtomOnResidue(residueIndex, 'C2');
    idx.N3 = atomicIndex.findAtomOnResidue(residueIndex, 'N3');
    idx.C4 = atomicIndex.findAtomOnResidue(residueIndex, 'C4');
    idx.C5 = atomicIndex.findAtomOnResidue(residueIndex, 'C5');
    if (idx.C5 === -1) {
        // modified ring, e.g. DP
        idx.C5 = atomicIndex.findAtomOnResidue(residueIndex, 'N5');
    }
    idx.C6 = atomicIndex.findAtomOnResidue(residueIndex, 'C6');
    idx.N7 = atomicIndex.findAtomOnResidue(residueIndex, 'N7');
    if (idx.N7 === -1) {
        // modified ring, e.g. DP
        idx.N7 = atomicIndex.findAtomOnResidue(residueIndex, 'C7');
    }
    idx.C8 = atomicIndex.findAtomOnResidue(residueIndex, 'C8');
    idx.N9 = atomicIndex.findAtomOnResidue(residueIndex, 'N9');

    return idx;
}

export function hasPurinIndices(idx: NucleicIndices): idx is NucleicIndices & { trace: ElementIndex, N1: ElementIndex, C2: ElementIndex, N3: ElementIndex, C4: ElementIndex, C5: ElementIndex, C6: ElementIndex, N7: ElementIndex, C8: ElementIndex, N9: ElementIndex } {
    return idx.trace !== -1 && idx.N1 !== -1 && idx.C2 !== -1 && idx.N3 !== -1 && idx.C4 !== -1 && idx.C5 !== -1 && idx.C6 !== -1 && idx.N7 !== -1 && idx.C8 !== -1 && idx.N9 !== -1;
}

export function setPyrimidineIndices(idx: NucleicIndices, unit: Unit.Atomic, residueIndex: ResidueIndex) {
    const atomicIndex = unit.model.atomicHierarchy.index;
    const { traceElementIndex } = unit.model.atomicHierarchy.derived.residue;

    idx.trace = traceElementIndex[residueIndex];
    idx.N1 = atomicIndex.findAtomOnResidue(residueIndex, 'N1');
    if (idx.N1 === -1) {
        // modified ring, e.g. DZ
        idx.N1 = atomicIndex.findAtomOnResidue(residueIndex, 'C1');
    }
    idx.C2 = atomicIndex.findAtomOnResidue(residueIndex, 'C2');
    idx.N3 = atomicIndex.findAtomOnResidue(residueIndex, 'N3');
    idx.C4 = atomicIndex.findAtomOnResidue(residueIndex, 'C4');
    idx.C5 = atomicIndex.findAtomOnResidue(residueIndex, 'C5');
    idx.C6 = atomicIndex.findAtomOnResidue(residueIndex, 'C6');

    return idx;
}

export function hasPyrimidineIndices(idx: NucleicIndices): idx is NucleicIndices & { trace: ElementIndex, N1: ElementIndex, C2: ElementIndex, N3: ElementIndex, C4: ElementIndex, C5: ElementIndex, C6: ElementIndex } {
    return idx.trace !== -1 && idx.N1 !== -1 && idx.C2 !== -1 && idx.N3 !== -1 && idx.C4 !== -1 && idx.C5 !== -1 && idx.C6 !== -1;
}

export function setSugarIndices(idx: NucleicIndices, unit: Unit.Atomic, residueIndex: ResidueIndex) {
    const atomicIndex = unit.model.atomicHierarchy.index;
    const { traceElementIndex } = unit.model.atomicHierarchy.derived.residue;

    idx.trace = traceElementIndex[residueIndex];
    idx.C1_1 = atomicIndex.findAtomOnResidue(residueIndex, "C1'");
    idx.C2_1 = atomicIndex.findAtomOnResidue(residueIndex, "C2'");
    idx.C3_1 = atomicIndex.findAtomOnResidue(residueIndex, "C3'");
    idx.C4_1 = atomicIndex.findAtomOnResidue(residueIndex, "C4'");
    idx.O4_1 = atomicIndex.findAtomOnResidue(residueIndex, "O4'");

    return idx;
}

export function hasSugarIndices(idx: NucleicIndices): idx is NucleicIndices & { trace: ElementIndex, C1_1: ElementIndex, C2_1: ElementIndex, C3_1: ElementIndex, C4_1: ElementIndex, O4_1: ElementIndex } {
    return idx.trace !== -1 && idx.C1_1 !== -1 && idx.C2_1 !== -1 && idx.C3_1 !== -1 && idx.C4_1 !== -1 && idx.O4_1 !== -1;
}
