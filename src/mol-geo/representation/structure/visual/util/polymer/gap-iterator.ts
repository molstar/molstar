/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, ElementIndex, ResidueIndex } from 'mol-model/structure';
import { AtomRole } from 'mol-model/structure/model/types';
import Iterator from 'mol-data/iterator';
import SortedRanges from 'mol-data/int/sorted-ranges';
import { getElementIndexForAtomRole } from 'mol-model/structure/util';
import { getGapRanges } from '../polymer';

/** Iterates over gaps, i.e. the stem residues/coarse elements adjacent to gaps */
export function PolymerGapIterator(unit: Unit): Iterator<PolymerGapPair> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return new AtomicPolymerGapIterator(unit)
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return new CoarsePolymerGapIterator(unit)
    }
}

interface PolymerGapPair {
    centerA: StructureElement
    centerB: StructureElement
}

function createPolymerGapPair (unit: Unit) {
    return {
        centerA: StructureElement.create(unit),
        centerB: StructureElement.create(unit),
    }
}

export class AtomicPolymerGapIterator implements Iterator<PolymerGapPair> {
    private value: PolymerGapPair
    private gapIt: SortedRanges.Iterator<ElementIndex, ResidueIndex>
    hasNext: boolean = false;

    private getElementIndex(residueIndex: ResidueIndex, atomRole: AtomRole) {
        const elementIndex = getElementIndexForAtomRole(this.unit.model, residueIndex, atomRole)
        return elementIndex === -1 ? this.unit.model.atomicHierarchy.residueAtomSegments.offsets[residueIndex] : elementIndex
    }

    move() {
        const { elements, residueIndex } = this.unit
        const gapSegment = this.gapIt.move();
        this.value.centerA.element = this.getElementIndex(residueIndex[elements[gapSegment.start]], 'trace')
        this.value.centerB.element = this.getElementIndex(residueIndex[elements[gapSegment.end - 1]], 'trace')
        this.hasNext = this.gapIt.hasNext
        return this.value;
    }

    constructor(private unit: Unit.Atomic) {
        this.gapIt = SortedRanges.transientSegments(getGapRanges(unit), unit.elements);
        this.value = createPolymerGapPair(unit)
        this.hasNext = this.gapIt.hasNext
    }
}

export class CoarsePolymerGapIterator implements Iterator<PolymerGapPair> {
    private value: PolymerGapPair
    private gapIt: SortedRanges.Iterator<ElementIndex, ElementIndex>
    hasNext: boolean = false;

    move() {
        const gapSegment = this.gapIt.move();
        this.value.centerA.element = this.unit.elements[gapSegment.start]
        this.value.centerB.element = this.unit.elements[gapSegment.end - 1]
        this.hasNext = this.gapIt.hasNext
        return this.value;
    }

    constructor(private unit: Unit.Spheres | Unit.Gaussians) {
        this.gapIt = SortedRanges.transientSegments(getGapRanges(unit), unit.elements);
        this.value = createPolymerGapPair(unit)
        this.hasNext = this.gapIt.hasNext
    }
}
