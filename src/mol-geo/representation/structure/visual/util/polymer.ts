/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, Model, ElementIndex, ResidueIndex } from 'mol-model/structure';
import { Segmentation, OrderedSet, Interval } from 'mol-data/int';
import { MoleculeType, SecondaryStructureType } from 'mol-model/structure/model/types';
import Iterator from 'mol-data/iterator';
import { Vec3 } from 'mol-math/linear-algebra';
import SortedRanges from 'mol-data/int/sorted-ranges';

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

function getMoleculeType(model: Model, residueIndex: number) {
    const compId = model.atomicHierarchy.residues.label_comp_id.value(residueIndex)
    const chemCompMap = model.properties.chemicalComponentMap
    const cc = chemCompMap.get(compId)
    return cc ? cc.moleculeType : MoleculeType.unknown
}

function getElementIndexForAtomId(unit: Unit.Atomic, residueSegment: Segmentation.Segment, atomId: string) {
    const elements = unit.elements
    const { label_atom_id } = unit.model.atomicHierarchy.atoms
    for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
        if (label_atom_id.value(elements[j]) === atomId) return j as ElementIndex
    }
    return residueSegment.end - 1 as ElementIndex
}

function getResidueTypeAtomIdElementIndex(unit: Unit.Atomic, residueSegment: Segmentation.Segment, type: 'trace' | 'direction') {
    const atomId = getResidueTypeAtomId(getMoleculeType(unit.model, residueSegment.index), type)
    return getElementIndexForAtomId(unit, residueSegment, atomId)
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

export class AtomicPolymerBackboneIterator<T extends number = number> implements Iterator<PolymerBackbonePair> {
    private value: PolymerBackbonePair
    private polymerIt: SortedRanges.Iterator<ElementIndex>
    private residueIt: Segmentation.SegmentIterator<ResidueIndex>
    private state: AtomicPolymerBackboneIteratorState = AtomicPolymerBackboneIteratorState.nextPolymer
    hasNext: boolean = false;

    move() {
        if (this.state === AtomicPolymerBackboneIteratorState.nextPolymer) {
            while (this.polymerIt.hasNext) {
                const residueSegment = this.polymerIt.move()
                this.residueIt.setSegment(residueSegment);
                if (this.residueIt.hasNext) {
                    this.value.centerB.element = getResidueTypeAtomIdElementIndex(this.unit, residueSegment, 'trace')
                    // setTraceElement(this.value.centerB, this.residueIt.move())
                    this.state = AtomicPolymerBackboneIteratorState.nextResidue
                    break
                }
            }
        }

        if (this.state === AtomicPolymerBackboneIteratorState.nextResidue) {
            this.value.centerA.element = this.value.centerB.element
            this.value.centerB.element = getResidueTypeAtomIdElementIndex(this.unit, this.residueIt.move(), 'trace')
            // setTraceElement(this.value.centerB, this.residueIt.move())
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

export class CoarsePolymerBackboneIterator<T extends number = number> implements Iterator<PolymerBackbonePair> {
    private value: PolymerBackbonePair
    private polymerIt: SortedRanges.Iterator<ElementIndex>
    private polymerSegment: Segmentation.Segment<ElementIndex>
    private state: CoarsePolymerBackboneIteratorState = CoarsePolymerBackboneIteratorState.nextPolymer
    private elementIndex: number
    hasNext: boolean = false;

    move() {
        if (this.state === CoarsePolymerBackboneIteratorState.nextPolymer) {
            if (this.polymerIt.hasNext) {
                this.polymerSegment = this.polymerIt.move();
                this.elementIndex = this.polymerSegment.start
                // this.elementIndex += 1
                if (this.elementIndex + 1 < this.polymerSegment.end) {
                    this.value.centerB.element = this.value.centerB.unit.elements[this.elementIndex]
                    this.state = CoarsePolymerBackboneIteratorState.nextElement
                } else {
                    this.state = CoarsePolymerBackboneIteratorState.nextPolymer
                }
            }
        }

        if (this.state === CoarsePolymerBackboneIteratorState.nextElement) {
            this.elementIndex += 1
            this.value.centerA.element = this.value.centerB.element
            this.value.centerB.element = this.value.centerB.unit.elements[this.elementIndex]
            if (this.elementIndex + 1 >= this.polymerSegment.end) {
                this.state = CoarsePolymerBackboneIteratorState.nextPolymer
            }
        }

        this.hasNext = this.elementIndex + 1 < this.polymerSegment.end || this.polymerIt.hasNext
        return this.value;
    }

    constructor(unit: Unit.Spheres | Unit.Gaussians) {
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
        d12: Vec3.zero(), d23: Vec3.zero(),
    }
}

const enum AtomicPolymerTraceIteratorState { nextPolymer, nextResidue }

function setSegment (outSegment: Segmentation.Segment<number>, index: number, segments: Segmentation<number>, min: number, max: number): Segmentation.Segment<number> {
    // index = Math.min(Math.max(0, index), segments.segments.length - 2)
    const _index = Math.min(Math.max(min, index), max)
    if (isNaN(_index)) console.log(_index, index, min, max)
    outSegment.index = _index
    outSegment.start = segments.offsets[_index]
    outSegment.end = segments.offsets[_index + 1]
    // console.log(index, {...outSegment}, {...boundingSegment}, segments.segments[boundingSegment.index])
    return outSegment
}

export class AtomicPolymerTraceIterator<T extends number = number> implements Iterator<PolymerTraceElement> {
    private value: PolymerTraceElement
    private polymerIt: SortedRanges.Iterator<ElementIndex>
    private residueIt: Segmentation.SegmentIterator<ResidueIndex>
    private residueAtomSegmentMin: number
    private residueAtomSegmentMax: number
    private state: AtomicPolymerTraceIteratorState = AtomicPolymerTraceIteratorState.nextPolymer
    private residueAtomSegments: Segmentation<ElementIndex, ResidueIndex>
    private tmpSegment: Segmentation.Segment<ResidueIndex>

    hasNext: boolean = false;

    private pos(target: Vec3, index: number) {
        target[0] = this.unit.model.atomicConformation.x[index]
        target[1] = this.unit.model.atomicConformation.y[index]
        target[2] = this.unit.model.atomicConformation.z[index]
    }

    updateResidueSegmentRange(polymerSegment: Segmentation.Segment<ElementIndex>) {
        const { polymerRanges, residueAtomSegments } = this.unit.model.atomicHierarchy
        const sMin = polymerRanges[polymerSegment.index * 2]
        const sMax = polymerRanges[polymerSegment.index * 2 + 1]
        this.residueAtomSegmentMin = residueAtomSegments.index[sMin]
        this.residueAtomSegmentMax = residueAtomSegments.index[sMax]
    }

    move() {
        const { residueIt, polymerIt, value } = this

        if (this.state === AtomicPolymerTraceIteratorState.nextPolymer) {
            while (polymerIt.hasNext) {
                const polymerSegment = polymerIt.move();
                // console.log('polymerSegment', {...polymerSegment})
                residueIt.setSegment(polymerSegment);
                this.updateResidueSegmentRange(polymerSegment)
                if (residueIt.hasNext) {
                    this.state = AtomicPolymerTraceIteratorState.nextResidue
                    break
                }
            }
        }

        if (this.state === AtomicPolymerTraceIteratorState.nextResidue) {
            const { tmpSegment, residueAtomSegments, residueAtomSegmentMin, residueAtomSegmentMax } = this
            const residueSegment = residueIt.move();
            const resSegIdx = residueSegment.index
            // console.log(residueLabel(this.unit.model, resSegIdx), resSegIdx, this.unit.model.properties.secondaryStructure.type[resSegIdx])
            value.center.element = getResidueTypeAtomIdElementIndex(this.unit, residueSegment, 'trace')

            setSegment(tmpSegment, resSegIdx - 2, residueAtomSegments, residueAtomSegmentMin, residueAtomSegmentMax)
            this.pos(value.t0, getResidueTypeAtomIdElementIndex(this.unit, tmpSegment, 'trace'))

            setSegment(tmpSegment, resSegIdx - 1, residueAtomSegments, residueAtomSegmentMin, residueAtomSegmentMax)
            this.pos(value.t1, getResidueTypeAtomIdElementIndex(this.unit, tmpSegment, 'trace'))
            this.pos(value.d12, getResidueTypeAtomIdElementIndex(this.unit, tmpSegment, 'direction'))

            setSegment(tmpSegment, resSegIdx, residueAtomSegments, residueAtomSegmentMin, residueAtomSegmentMax)
            value.secStrucType = this.unit.model.properties.secondaryStructure.type[resSegIdx]
            this.pos(value.t2, getResidueTypeAtomIdElementIndex(this.unit, tmpSegment, 'trace'))
            this.pos(value.d23, getResidueTypeAtomIdElementIndex(this.unit, tmpSegment, 'direction'))

            setSegment(tmpSegment, resSegIdx + 1, residueAtomSegments, residueAtomSegmentMin, residueAtomSegmentMax)
            this.pos(value.t3, getResidueTypeAtomIdElementIndex(this.unit, tmpSegment, 'trace'))

            setSegment(tmpSegment, resSegIdx + 2, residueAtomSegments, residueAtomSegmentMin, residueAtomSegmentMax)
            this.pos(value.t4, getResidueTypeAtomIdElementIndex(this.unit, tmpSegment, 'trace'))

            value.first = resSegIdx === residueAtomSegmentMin
            value.last = resSegIdx === residueAtomSegmentMax

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
        this.tmpSegment = { index: 0 as ResidueIndex, start: 0 as ElementIndex, end: 0 as ElementIndex }
        this.hasNext = this.residueIt.hasNext && this.polymerIt.hasNext
    }
}

export class CoarsePolymerTraceIterator<T extends number = number> implements Iterator<PolymerTraceElement> {
    private value: PolymerTraceElement

    hasNext: boolean = false;

    move() {
        return this.value;
    }

    constructor(unit: Unit.Spheres | Unit.Gaussians) {
        this.value = createPolymerTraceElement(unit)
        this.hasNext = false
    }
}