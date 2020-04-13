/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure, StructureElement, ElementIndex, ResidueIndex } from '../../../../../mol-model/structure';
import { Segmentation } from '../../../../../mol-data/int';
import Iterator from '../../../../../mol-data/iterator';
import SortedRanges from '../../../../../mol-data/int/sorted-ranges';
import { getPolymerRanges } from '../polymer';

/** Iterates over consecutive pairs of residues/coarse elements in polymers */
export function PolymerBackboneIterator(structure: Structure, unit: Unit): Iterator<PolymerBackbonePair> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return new AtomicPolymerBackboneIterator(structure, unit);
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return new CoarsePolymerBackboneIterator(structure, unit);
    }
}

interface PolymerBackbonePair {
    centerA: StructureElement.Location
    centerB: StructureElement.Location
}

function createPolymerBackbonePair (structure: Structure, unit: Unit) {
    return {
        centerA: StructureElement.Location.create(structure, unit),
        centerB: StructureElement.Location.create(structure, unit),
    };
}

const enum AtomicPolymerBackboneIteratorState { nextPolymer, firstResidue, nextResidue, cycle }

export class AtomicPolymerBackboneIterator implements Iterator<PolymerBackbonePair> {
    private traceElementIndex: ArrayLike<ElementIndex>
    private value: PolymerBackbonePair
    private polymerIt: SortedRanges.Iterator<ElementIndex, ResidueIndex>
    private residueIt: Segmentation.SegmentIterator<ResidueIndex>
    private state: AtomicPolymerBackboneIteratorState = AtomicPolymerBackboneIteratorState.nextPolymer
    private residueSegment: Segmentation.Segment<ResidueIndex>
    hasNext: boolean = false;

    move() {
        if (this.state === AtomicPolymerBackboneIteratorState.nextPolymer) {
            while (this.polymerIt.hasNext) {
                this.residueIt.setSegment(this.polymerIt.move());
                if (this.residueIt.hasNext) {
                    this.residueSegment = this.residueIt.move();
                    this.value.centerB.element = this.traceElementIndex[this.residueSegment.index];
                    this.state = AtomicPolymerBackboneIteratorState.nextResidue;
                    break;
                }
            }
        }

        if (this.state === AtomicPolymerBackboneIteratorState.nextResidue) {
            this.residueSegment = this.residueIt.move();
            this.value.centerA.element = this.value.centerB.element;
            this.value.centerB.element = this.traceElementIndex[this.residueSegment.index];
            if (!this.residueIt.hasNext) {
                if (this.unit.model.atomicRanges.cyclicPolymerMap.has(this.residueSegment.index)) {
                    this.state = AtomicPolymerBackboneIteratorState.cycle;
                } else {
                    // TODO need to advance to a polymer that has two or more residues (can't assume it has)
                    this.state = AtomicPolymerBackboneIteratorState.nextPolymer;
                }
            }
        } else if (this.state === AtomicPolymerBackboneIteratorState.cycle) {
            const { cyclicPolymerMap } = this.unit.model.atomicRanges;
            this.value.centerA.element = this.value.centerB.element;
            this.value.centerB.element = this.traceElementIndex[cyclicPolymerMap.get(this.residueSegment.index)!];
            // TODO need to advance to a polymer that has two or more residues (can't assume it has)
            this.state = AtomicPolymerBackboneIteratorState.nextPolymer;
        }

        this.hasNext = this.residueIt.hasNext || this.polymerIt.hasNext || this.state === AtomicPolymerBackboneIteratorState.cycle;
        return this.value;
    }

    constructor(structure: Structure, private unit: Unit.Atomic) {
        this.traceElementIndex = unit.model.atomicHierarchy.derived.residue.traceElementIndex as ArrayLike<ElementIndex>; // can assume it won't be -1 for polymer residues
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);
        this.residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
        this.value = createPolymerBackbonePair(structure, unit);
        this.hasNext = this.residueIt.hasNext && this.polymerIt.hasNext;
    }
}

const enum CoarsePolymerBackboneIteratorState { nextPolymer, nextElement }

export class CoarsePolymerBackboneIterator implements Iterator<PolymerBackbonePair> {
    private value: PolymerBackbonePair
    private polymerIt: SortedRanges.Iterator<ElementIndex, ResidueIndex>
    private polymerSegment: Segmentation.Segment<ResidueIndex>
    private state: CoarsePolymerBackboneIteratorState = CoarsePolymerBackboneIteratorState.nextPolymer
    private elementIndex: number
    hasNext: boolean = false;

    move() {
        if (this.state === CoarsePolymerBackboneIteratorState.nextPolymer) {
            if (this.polymerIt.hasNext) {
                this.polymerSegment = this.polymerIt.move();
                this.elementIndex = this.polymerSegment.start;
                if (this.elementIndex + 1 < this.polymerSegment.end) {
                    this.value.centerB.element = this.unit.elements[this.elementIndex];
                    this.state = CoarsePolymerBackboneIteratorState.nextElement;
                } else {
                    this.state = CoarsePolymerBackboneIteratorState.nextPolymer;
                }
            }
        }

        if (this.state === CoarsePolymerBackboneIteratorState.nextElement) {
            this.elementIndex += 1;
            this.value.centerA.element = this.value.centerB.element;
            this.value.centerB.element = this.unit.elements[this.elementIndex];
            if (this.elementIndex + 1 >= this.polymerSegment.end) {
                this.state = CoarsePolymerBackboneIteratorState.nextPolymer;
            }
        }

        this.hasNext = this.elementIndex + 1 < this.polymerSegment.end || this.polymerIt.hasNext;
        return this.value;
    }

    constructor(structure: Structure, private unit: Unit.Spheres | Unit.Gaussians) {
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);
        this.value = createPolymerBackbonePair(structure, unit);
        this.hasNext = this.polymerIt.hasNext;
    }
}