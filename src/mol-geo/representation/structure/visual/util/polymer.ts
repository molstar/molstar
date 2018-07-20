/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, ElementIndex, ResidueIndex } from 'mol-model/structure';
import { Segmentation, OrderedSet, Interval, SortedArray } from 'mol-data/int';
import { MoleculeType, SecondaryStructureType, AtomRole } from 'mol-model/structure/model/types';
import Iterator from 'mol-data/iterator';
import { Vec3 } from 'mol-math/linear-algebra';
import SortedRanges from 'mol-data/int/sorted-ranges';
import { CoarseSphereConformation, CoarseGaussianConformation } from 'mol-model/structure/model/properties/coarse';
import { getMoleculeType, getElementIndexForAtomRole } from 'mol-model/structure/util';

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

export function getPolymerElementCount(unit: Unit) {
    let count = 0
    const { elements } = unit
    const polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), elements)
    switch (unit.kind) {
        case Unit.Kind.Atomic:
            const residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements)
            while (polymerIt.hasNext) {
                const polymerSegment = polymerIt.move()
                residueIt.setSegment(polymerSegment)
                while (residueIt.hasNext) {
                    const residueSegment = residueIt.move()
                    const { start, end } = residueSegment
                    if (OrderedSet.areIntersecting(Interval.ofBounds(elements[start], elements[end - 1]), elements)) ++count
                }
            }
            break
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            while (polymerIt.hasNext) {
                const { start, end } = polymerIt.move()
                count += OrderedSet.intersectionSize(Interval.ofBounds(elements[start], elements[end - 1]), elements)
            }
            break
    }
    return count
}

export function getPolymerGapCount(unit: Unit) {
    let count = 0
    const { elements } = unit
    const gapIt = SortedRanges.transientSegments(getGapRanges(unit), elements)
    while (gapIt.hasNext) {
        const { start, end } = gapIt.move()
        if (OrderedSet.areIntersecting(Interval.ofBounds(elements[start], elements[end - 1]), elements)) ++count
    }
    return count
}



/** Iterates over consecutive pairs of residues/coarse elements in polymers */
export function PolymerBackboneIterator(unit: Unit): Iterator<PolymerBackbonePair> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return new AtomicPolymerBackboneIterator(unit)
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return new CoarsePolymerBackboneIterator(unit)
    }
}

interface PolymerBackbonePair {
    centerA: StructureElement
    centerB: StructureElement
}

function createPolymerBackbonePair (unit: Unit) {
    return {
        centerA: StructureElement.create(unit),
        centerB: StructureElement.create(unit),
    }
}

const enum AtomicPolymerBackboneIteratorState { nextPolymer, firstResidue, nextResidue, cycle }

export class AtomicPolymerBackboneIterator implements Iterator<PolymerBackbonePair> {
    private value: PolymerBackbonePair
    private polymerIt: SortedRanges.Iterator<ElementIndex, ResidueIndex>
    private residueIt: Segmentation.SegmentIterator<ResidueIndex>
    private state: AtomicPolymerBackboneIteratorState = AtomicPolymerBackboneIteratorState.nextPolymer
    private residueSegment: Segmentation.Segment<ResidueIndex>
    hasNext: boolean = false;

    private getElementIndex(residueIndex: ResidueIndex, atomRole: AtomRole) {
        const index = getElementIndexForAtomRole(this.unit.model, residueIndex, atomRole)
        // TODO handle case when it returns -1
        const elementIndex = SortedArray.indexOf(this.unit.elements, index) as ElementIndex
        if (elementIndex === -1) {
            console.log('-1', residueIndex, atomRole, index)
        }
        return elementIndex === -1 ? 0 as ElementIndex : elementIndex
    }

    move() {
        if (this.state === AtomicPolymerBackboneIteratorState.nextPolymer) {
            while (this.polymerIt.hasNext) {
                this.residueIt.setSegment(this.polymerIt.move());
                if (this.residueIt.hasNext) {
                    this.residueSegment = this.residueIt.move()
                    this.value.centerB.element = this.getElementIndex(this.residueSegment.index, 'trace')
                    this.state = AtomicPolymerBackboneIteratorState.nextResidue
                    break
                }
            }
        }

        if (this.state === AtomicPolymerBackboneIteratorState.nextResidue) {
            this.residueSegment = this.residueIt.move()
            this.value.centerA.element = this.value.centerB.element
            this.value.centerB.element = this.getElementIndex(this.residueSegment.index, 'trace')
            if (!this.residueIt.hasNext) {
                if (this.unit.model.atomicHierarchy.cyclicPolymerMap.has(this.residueSegment.index)) {
                    this.state = AtomicPolymerBackboneIteratorState.cycle
                } else {
                    // TODO need to advance to a polymer that has two or more residues (can't assume it has)
                    this.state = AtomicPolymerBackboneIteratorState.nextPolymer
                }
            }
        } else if (this.state === AtomicPolymerBackboneIteratorState.cycle) {
            const { cyclicPolymerMap } = this.unit.model.atomicHierarchy
            this.value.centerA.element = this.value.centerB.element
            this.value.centerB.element = this.getElementIndex(cyclicPolymerMap.get(this.residueSegment.index)!, 'trace')
            // TODO need to advance to a polymer that has two or more residues (can't assume it has)
            this.state = AtomicPolymerBackboneIteratorState.nextPolymer
        }

        this.hasNext = this.residueIt.hasNext || this.polymerIt.hasNext || this.state === AtomicPolymerBackboneIteratorState.cycle
        return this.value;
    }

    constructor(private unit: Unit.Atomic) {
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements)
        this.residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements)
        this.value = createPolymerBackbonePair(unit)
        this.hasNext = this.residueIt.hasNext && this.polymerIt.hasNext
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
                this.elementIndex = this.polymerSegment.start
                if (this.elementIndex + 1 < this.polymerSegment.end) {
                    this.value.centerB.element = this.unit.elements[this.elementIndex]
                    this.state = CoarsePolymerBackboneIteratorState.nextElement
                } else {
                    this.state = CoarsePolymerBackboneIteratorState.nextPolymer
                }
            }
        }

        if (this.state === CoarsePolymerBackboneIteratorState.nextElement) {
            this.elementIndex += 1
            this.value.centerA.element = this.value.centerB.element
            this.value.centerB.element = this.unit.elements[this.elementIndex]
            if (this.elementIndex + 1 >= this.polymerSegment.end) {
                this.state = CoarsePolymerBackboneIteratorState.nextPolymer
            }
        }

        this.hasNext = this.elementIndex + 1 < this.polymerSegment.end || this.polymerIt.hasNext
        return this.value;
    }

    constructor(private unit: Unit.Spheres | Unit.Gaussians) {
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);
        this.value = createPolymerBackbonePair(unit)
        this.hasNext = this.polymerIt.hasNext
    }
}

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
        const index = getElementIndexForAtomRole(this.unit.model, residueIndex, atomRole)
        // TODO handle case when it returns -1
        const elementIndex = SortedArray.indexOf(this.unit.elements, index) as ElementIndex
        if (elementIndex === -1) {
            console.log('-1', residueIndex, atomRole, index)
        }
        return elementIndex === -1 ? 0 as ElementIndex : elementIndex
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

/**
 * Iterates over individual residues/coarse elements in polymers of a unit while
 * providing information about the neighbourhood in the underlying model for drawing splines
 */
export function PolymerTraceIterator(unit: Unit): Iterator<PolymerTraceElement> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return new AtomicPolymerTraceIterator(unit)
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return new CoarsePolymerTraceIterator(unit)
    }
}

interface PolymerTraceElement {
    center: StructureElement
    first: boolean, last: boolean
    secStrucType: SecondaryStructureType
    secStrucChange: boolean
    moleculeType: MoleculeType
    t0: Vec3, t1: Vec3, t2: Vec3, t3: Vec3, t4: Vec3
    d12: Vec3, d23: Vec3
}

function createPolymerTraceElement (unit: Unit): PolymerTraceElement {
    return {
        center: StructureElement.create(unit),
        first: false, last: false,
        secStrucType: SecondaryStructureType.create(SecondaryStructureType.Flag.NA),
        secStrucChange: false,
        moleculeType: MoleculeType.unknown,
        t0: Vec3.zero(), t1: Vec3.zero(), t2: Vec3.zero(), t3: Vec3.zero(), t4: Vec3.zero(),
        d12: Vec3.create(1, 0, 0), d23: Vec3.create(1, 0, 0),
    }
}

const enum AtomicPolymerTraceIteratorState { nextPolymer, nextResidue }

export class AtomicPolymerTraceIterator implements Iterator<PolymerTraceElement> {
    private value: PolymerTraceElement
    private polymerIt: SortedRanges.Iterator<ElementIndex, ResidueIndex>
    private residueIt: Segmentation.SegmentIterator<ResidueIndex>
    private polymerSegment: Segmentation.Segment<ResidueIndex>
    private residueSegmentMin: ResidueIndex
    private residueSegmentMax: ResidueIndex
    private state: AtomicPolymerTraceIteratorState = AtomicPolymerTraceIteratorState.nextPolymer
    private residueAtomSegments: Segmentation<ElementIndex, ResidueIndex>

    private p0 = Vec3.zero();
    private p1 = Vec3.zero();
    private p2 = Vec3.zero();
    private p3 = Vec3.zero();
    private p4 = Vec3.zero();
    private p5 = Vec3.zero();
    private p6 = Vec3.zero();

    // private v01 = Vec3.zero();
    private v12 = Vec3.zero();
    private v23 = Vec3.zero();
    // private v34 = Vec3.zero();

    hasNext: boolean = false;

    private pos(target: Vec3, index: number) {
        target[0] = this.unit.model.atomicConformation.x[index]
        target[1] = this.unit.model.atomicConformation.y[index]
        target[2] = this.unit.model.atomicConformation.z[index]
    }

    private updateResidueSegmentRange(polymerSegment: Segmentation.Segment<ResidueIndex>) {
        const { index } = this.unit.model.atomicHierarchy.residueAtomSegments
        this.residueSegmentMin = index[this.unit.elements[polymerSegment.start]]
        this.residueSegmentMax = index[this.unit.elements[polymerSegment.end - 1]]
    }

    private getAtomIndex(residueIndex: ResidueIndex, atomRole: AtomRole) {
        const { cyclicPolymerMap } = this.unit.model.atomicHierarchy
        if (residueIndex < this.residueSegmentMin) {
            const cyclicIndex = cyclicPolymerMap.get(this.residueSegmentMin)
            if (cyclicIndex !== undefined) {

                residueIndex = cyclicIndex - (this.residueSegmentMin - residueIndex - 1) as ResidueIndex
            } else {
                residueIndex = this.residueSegmentMin
            }
        } else if (residueIndex > this.residueSegmentMax) {
            const cyclicIndex = cyclicPolymerMap.get(this.residueSegmentMax)
            if (cyclicIndex !== undefined) {
                residueIndex = cyclicIndex + (residueIndex - this.residueSegmentMax - 1) as ResidueIndex
            } else {
                residueIndex = this.residueSegmentMax
            }
        }
        return getElementIndexForAtomRole(this.unit.model, residueIndex as ResidueIndex, atomRole)
    }

    private getElementIndex(residueIndex: ResidueIndex, atomRole: AtomRole) {
        const index = this.getAtomIndex(residueIndex, atomRole)
        // TODO handle case when it returns -1
        const elementIndex = SortedArray.indexOf(this.unit.elements, index) as ElementIndex
        if (elementIndex === -1) {
            console.log('-1', residueIndex, atomRole, index)
        }
        return elementIndex === -1 ? 0 as ElementIndex : elementIndex
    }

    private setControlPoint(out: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, residueIndex: ResidueIndex) {
        const ss = this.unit.model.properties.secondaryStructure.type[residueIndex]
        if (SecondaryStructureType.is(ss, SecondaryStructureType.Flag.Beta)) {
            Vec3.scale(out, Vec3.add(out, p1, Vec3.add(out, p3, Vec3.add(out, p2, p2))), 1/4)
        } else {
            Vec3.copy(out, p2)
        }
    }

    move() {
        const { residueIt, polymerIt, value } = this

        if (this.state === AtomicPolymerTraceIteratorState.nextPolymer) {
            while (polymerIt.hasNext) {
                this.polymerSegment = polymerIt.move();
                residueIt.setSegment(this.polymerSegment);
                this.updateResidueSegmentRange(this.polymerSegment)
                if (residueIt.hasNext) {
                    this.state = AtomicPolymerTraceIteratorState.nextResidue
                    break
                }
            }
        }

        if (this.state === AtomicPolymerTraceIteratorState.nextResidue) {
            const { index: residueIndex } = residueIt.move();
            value.center.element = this.getElementIndex(residueIndex, 'trace')

            this.pos(this.p0, this.getAtomIndex(residueIndex - 3 as ResidueIndex, 'trace'))
            this.pos(this.p1, this.getAtomIndex(residueIndex - 2 as ResidueIndex, 'trace'))
            this.pos(this.p2, this.getAtomIndex(residueIndex - 1 as ResidueIndex, 'trace'))
            this.pos(this.p3, this.getAtomIndex(residueIndex, 'trace'))
            this.pos(this.p4, this.getAtomIndex(residueIndex + 1 as ResidueIndex, 'trace'))
            this.pos(this.p5, this.getAtomIndex(residueIndex + 2 as ResidueIndex, 'trace'))
            this.pos(this.p6, this.getAtomIndex(residueIndex + 3 as ResidueIndex, 'trace'))

            // this.pos(this.v01, this.getAtomIndex(residueIndex - 2 as ResidueIndex, 'direction'))
            this.pos(this.v12, this.getAtomIndex(residueIndex - 1 as ResidueIndex, 'direction'))
            this.pos(this.v23, this.getAtomIndex(residueIndex, 'direction'))
            // this.pos(this.v34, this.getAtomIndex(residueIndex + 1 as ResidueIndex, 'direction'))

            this.value.secStrucType = this.unit.model.properties.secondaryStructure.type[residueIndex]

            this.setControlPoint(value.t0, this.p0, this.p1, this.p2, residueIndex - 2 as ResidueIndex)
            this.setControlPoint(value.t1, this.p1, this.p2, this.p3, residueIndex - 1 as ResidueIndex)
            this.setControlPoint(value.t2, this.p2, this.p3, this.p4, residueIndex)
            this.setControlPoint(value.t3, this.p3, this.p4, this.p5, residueIndex + 1 as ResidueIndex)
            this.setControlPoint(value.t4, this.p4, this.p5, this.p6, residueIndex + 2 as ResidueIndex)

            Vec3.copy(value.d12, this.v12)
            Vec3.copy(value.d23, this.v23)

            value.first = residueIndex === this.residueSegmentMin
            value.last = residueIndex === this.residueSegmentMax
            value.secStrucChange = this.unit.model.properties.secondaryStructure.key[residueIndex] !== this.unit.model.properties.secondaryStructure.key[residueIndex + 1]
            value.moleculeType = getMoleculeType(this.unit.model, residueIndex)

            if (!residueIt.hasNext) {
                this.state = AtomicPolymerTraceIteratorState.nextPolymer
            }
        }

        this.hasNext = residueIt.hasNext || polymerIt.hasNext

        return this.value;
    }

    constructor(private unit: Unit.Atomic) {
        this.residueAtomSegments = unit.model.atomicHierarchy.residueAtomSegments
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements)
        this.residueIt = Segmentation.transientSegments(this.residueAtomSegments, unit.elements);
        this.value = createPolymerTraceElement(unit)
        this.hasNext = this.residueIt.hasNext && this.polymerIt.hasNext
    }
}

const enum CoarsePolymerTraceIteratorState { nextPolymer, nextElement }

export class CoarsePolymerTraceIterator implements Iterator<PolymerTraceElement> {
    private value: PolymerTraceElement
    private polymerIt: SortedRanges.Iterator<ElementIndex, ResidueIndex>
    private polymerSegment: Segmentation.Segment<ResidueIndex>
    private state: CoarsePolymerTraceIteratorState = CoarsePolymerTraceIteratorState.nextPolymer
    private conformation: CoarseSphereConformation | CoarseGaussianConformation
    private elementIndex: number
    hasNext: boolean = false;

    private pos(target: Vec3, elementIndex: number) {
        elementIndex = Math.min(Math.max(this.polymerSegment.start, elementIndex), this.polymerSegment.end - 1)
        const index = this.unit.elements[elementIndex]
        target[0] = this.conformation.x[index]
        target[1] = this.conformation.y[index]
        target[2] = this.conformation.z[index]
    }

    move() {
        if (this.state === CoarsePolymerTraceIteratorState.nextPolymer) {
            while (this.polymerIt.hasNext) {
                this.polymerSegment = this.polymerIt.move();
                this.elementIndex = this.polymerSegment.start

                if (this.elementIndex + 1 < this.polymerSegment.end) {
                    this.state = CoarsePolymerTraceIteratorState.nextElement
                    break
                }
            }
        }

        if (this.state === CoarsePolymerTraceIteratorState.nextElement) {
            this.elementIndex += 1
            this.value.center.element = this.value.center.unit.elements[this.elementIndex]

            this.pos(this.value.t0, this.elementIndex - 2)
            this.pos(this.value.t1, this.elementIndex - 1)
            this.pos(this.value.t2, this.elementIndex)
            this.pos(this.value.t3, this.elementIndex + 1)
            this.pos(this.value.t4, this.elementIndex + 2)

            this.value.first = this.elementIndex === this.polymerSegment.start
            this.value.last = this.elementIndex === this.polymerSegment.end - 1

            if (this.elementIndex + 1 >= this.polymerSegment.end) {
                this.state = CoarsePolymerTraceIteratorState.nextPolymer
            }
        }

        this.hasNext = this.elementIndex + 1 < this.polymerSegment.end || this.polymerIt.hasNext
        return this.value;
    }

    constructor(private unit: Unit.Spheres | Unit.Gaussians) {
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);
        this.value = createPolymerTraceElement(unit)
        switch (unit.kind) {
            case Unit.Kind.Spheres: this.conformation = unit.model.coarseConformation.spheres; break
            case Unit.Kind.Gaussians: this.conformation = unit.model.coarseConformation.gaussians; break
        }
        this.hasNext = this.polymerIt.hasNext
    }
}