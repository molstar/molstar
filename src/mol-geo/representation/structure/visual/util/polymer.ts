/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, Model, ElementIndex, ResidueIndex } from 'mol-model/structure';
import { Segmentation, OrderedSet, Interval, SortedArray } from 'mol-data/int';
import { MoleculeType, SecondaryStructureType } from 'mol-model/structure/model/types';
import Iterator from 'mol-data/iterator';
import { Vec3 } from 'mol-math/linear-algebra';
import SortedRanges from 'mol-data/int/sorted-ranges';
import { CoarseSphereConformation, CoarseGaussianConformation } from 'mol-model/structure/model/properties/coarse';

export function getPolymerRanges(unit: Unit): SortedRanges<ElementIndex> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return unit.model.atomicHierarchy.polymerRanges
        case Unit.Kind.Spheres: return unit.model.coarseHierarchy.spheres.polymerRanges
        case Unit.Kind.Gaussians: return unit.model.coarseHierarchy.gaussians.polymerRanges
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

function getResidueTypeAtomId(moleculeType: MoleculeType, atomType: 'trace' | 'direction') {
    switch (moleculeType) {
        case MoleculeType.protein:
            switch (atomType) {
                case 'trace': return 'CA'
                case 'direction': return 'O'
            }
            break
        case MoleculeType.RNA:
            switch (atomType) {
                case 'trace': return 'C4\''
                case 'direction': return 'C3\''
            }
            break
        case MoleculeType.DNA:
            switch (atomType) {
                case 'trace': return 'C3\''
                case 'direction': return 'O4\''
            }
            break
    }
    return ''
}

function getMoleculeType(model: Model, residueIndex: ResidueIndex) {
    const compId = model.atomicHierarchy.residues.label_comp_id.value(residueIndex)
    const chemCompMap = model.properties.chemicalComponentMap
    const cc = chemCompMap.get(compId)
    return cc ? cc.moleculeType : MoleculeType.unknown
}

function getElementIndexForAtomId(model: Model, rI: ResidueIndex, atomId: string): ElementIndex {
    const { offsets } = model.atomicHierarchy.residueAtomSegments
    const { label_atom_id } = model.atomicHierarchy.atoms
    for (let j = offsets[rI], _j = offsets[rI + 1]; j < _j; j++) {
        if (label_atom_id.value(j) === atomId) return j as ElementIndex
    }
    return offsets[rI] as ElementIndex
}

function getElementIndexForResidueTypeAtomId(model: Model, rI: ResidueIndex, atomType: 'trace' | 'direction') {
    const atomId = getResidueTypeAtomId(getMoleculeType(model, rI), atomType)
    return getElementIndexForAtomId(model, rI, atomId)
}

// function residueLabel(model: Model, rI: number) {
//     const { residues, chains, residueSegments, chainSegments } = model.atomicHierarchy
//     const { label_comp_id, label_seq_id } = residues
//     const { label_asym_id } = chains
//     const cI = chainSegments.segmentMap[residueSegments.segments[rI]]
//     return `${label_asym_id.value(cI)} ${label_comp_id.value(rI)} ${label_seq_id.value(rI)}`
// }

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

const enum AtomicPolymerBackboneIteratorState { nextPolymer, firstResidue, nextResidue }

export class AtomicPolymerBackboneIterator implements Iterator<PolymerBackbonePair> {
    private value: PolymerBackbonePair
    private polymerIt: SortedRanges.Iterator<ElementIndex, ResidueIndex>
    private residueIt: Segmentation.SegmentIterator<ResidueIndex>
    private state: AtomicPolymerBackboneIteratorState = AtomicPolymerBackboneIteratorState.nextPolymer
    hasNext: boolean = false;

    getElementIndex(residueIndex: ResidueIndex, atomType: 'trace' | 'direction') {
        const index = getElementIndexForResidueTypeAtomId(this.unit.model, residueIndex, atomType)
        // // TODO handle case when it returns -1
        // return SortedArray.indexOf(this.unit.elements, index) as ElementIndex

        const elementIndex = SortedArray.indexOf(this.unit.elements, index) as ElementIndex
        if (elementIndex === -1) {
            console.log('-1', residueIndex, atomType, index)
        }
        return elementIndex === -1 ? 0 as ElementIndex : elementIndex
    }

    move() {
        if (this.state === AtomicPolymerBackboneIteratorState.nextPolymer) {
            while (this.polymerIt.hasNext) {
                const residueSegment = this.polymerIt.move()
                this.residueIt.setSegment(residueSegment);
                if (this.residueIt.hasNext) {
                    this.value.centerB.element = this.getElementIndex(this.residueIt.move().index, 'trace')
                    this.state = AtomicPolymerBackboneIteratorState.nextResidue
                    break
                }
            }
        }

        if (this.state === AtomicPolymerBackboneIteratorState.nextResidue) {
            this.value.centerA.element = this.value.centerB.element
            this.value.centerB.element = this.getElementIndex(this.residueIt.move().index, 'trace')
            if (!this.residueIt.hasNext) {
                // TODO need to advance to a polymer that has two or more residues (can't assume it has)
                this.state = AtomicPolymerBackboneIteratorState.nextPolymer
            }
        }

        this.hasNext = this.residueIt.hasNext || this.polymerIt.hasNext
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
    t0: Vec3, t1: Vec3, t2: Vec3, t3: Vec3, t4: Vec3
    d12: Vec3, d23: Vec3
}

function createPolymerTraceElement (unit: Unit): PolymerTraceElement {
    return {
        center: StructureElement.create(unit),
        first: false, last: false,
        secStrucType: SecondaryStructureType.create(SecondaryStructureType.Flag.NA),
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

    hasNext: boolean = false;

    private pos(target: Vec3, index: number) {
        target[0] = this.unit.model.atomicConformation.x[index]
        target[1] = this.unit.model.atomicConformation.y[index]
        target[2] = this.unit.model.atomicConformation.z[index]
    }

    updateResidueSegmentRange(polymerSegment: Segmentation.Segment<ResidueIndex>) {
        const { index } = this.unit.model.atomicHierarchy.residueAtomSegments
        this.residueSegmentMin = index[this.unit.elements[polymerSegment.start]]
        this.residueSegmentMax = index[this.unit.elements[polymerSegment.end - 1]]
    }

    getAtomIndex(residueIndex: number, atomType: 'trace' | 'direction') {
        const index = Math.min(Math.max(this.residueSegmentMin, residueIndex), this.residueSegmentMax)
        return getElementIndexForResidueTypeAtomId(this.unit.model, index as ResidueIndex, atomType)
    }

    getElementIndex(residueIndex: number, atomType: 'trace' | 'direction') {
        const index = this.getAtomIndex(residueIndex, atomType)
        // TODO handle case when it returns -1
        const elementIndex = SortedArray.indexOf(this.unit.elements, index) as ElementIndex
        if (elementIndex === -1) {
            console.log('-1', residueIndex, atomType, index)
        }
        return elementIndex === -1 ? 0 as ElementIndex : elementIndex
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

            this.pos(value.t0, this.getAtomIndex(residueIndex - 2, 'trace'))
            this.pos(value.t1, this.getAtomIndex(residueIndex - 1, 'trace'))
            this.pos(value.t2, this.getAtomIndex(residueIndex, 'trace'))
            this.pos(value.t3, this.getAtomIndex(residueIndex + 1, 'trace'))
            this.pos(value.t4, this.getAtomIndex(residueIndex + 2, 'trace'))

            this.pos(value.d12, this.getAtomIndex(residueIndex - 1, 'direction'))
            this.pos(value.d23, this.getAtomIndex(residueIndex, 'direction'))

            this.value.secStrucType = this.unit.model.properties.secondaryStructure.type[residueIndex]

            value.first = residueIndex === this.polymerSegment.start
            value.last = residueIndex === this.polymerSegment.end - 1

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