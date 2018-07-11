/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, StructureProperties, ElementIndex } from 'mol-model/structure';
import { Segmentation } from 'mol-data/int';
import { MoleculeType } from 'mol-model/structure/model/types';
import Iterator from 'mol-data/iterator';
import { SegmentIterator } from 'mol-data/int/impl/segmentation';
import { Vec3 } from 'mol-math/linear-algebra';
import { SymmetryOperator } from 'mol-math/geometry';

export function getPolymerElementCount(unit: Unit) {
    let count = 0
    const { elements } = unit
    const l = StructureElement.create(unit)
    if (Unit.isAtomic(unit)) {
        const { polymerAtomSegments, residueAtomSegments } = unit.model.atomicHierarchy
        const polymerIt = Segmentation.transientSegments(polymerAtomSegments, elements);
        const residuesIt = Segmentation.transientSegments(residueAtomSegments, elements);
        while (polymerIt.hasNext) {
            residuesIt.setSegment(polymerIt.move());
            while (residuesIt.hasNext) {
                residuesIt.move();
                count++
            }
        }
    } else if (Unit.isSpheres(unit)) {
        for (let i = 0, il = elements.length; i < il; ++i) {
            l.element = elements[i]
            if (StructureProperties.entity.type(l) === 'polymer') count++
        }
    }
    return count
}

function getTraceName(l: StructureElement) {
    const compId = StructureProperties.residue.label_comp_id(l)
    const chemCompMap = l.unit.model.properties.chemicalComponentMap
    const cc = chemCompMap.get(compId)
    const moleculeType = cc ? cc.moleculeType : MoleculeType.unknown
    let traceName = ''
    if (moleculeType === MoleculeType.protein) {
        traceName = 'CA'
    } else if (moleculeType === MoleculeType.DNA || moleculeType === MoleculeType.RNA) {
        traceName = 'P'
    }
    return traceName
}

function setTraceElement(l: StructureElement, residueSegment: Segmentation.Segment<ElementIndex>) {
    const elements = l.unit.elements
    l.element = elements[residueSegment.start]
    const traceName = getTraceName(l)

    for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
        l.element = elements[j];
        if (StructureProperties.atom.label_atom_id(l) === traceName) return j
    }
    return residueSegment.end - 1
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
    indexA: number
    indexB: number
    posA: Vec3
    posB: Vec3
}

function createPolymerBackbonePair (unit: Unit) {
    return {
        centerA: StructureElement.create(unit),
        centerB: StructureElement.create(unit),
        indexA: 0,
        indexB: 0,
        posA: Vec3.zero(),
        posB: Vec3.zero()
    }
}

const enum AtomicPolymerBackboneIteratorState { nextPolymer, firstResidue, nextResidue }

export class AtomicPolymerBackboneIterator<T extends number = number> implements Iterator<PolymerBackbonePair> {
    private value: PolymerBackbonePair

    private polymerIt: SegmentIterator<ElementIndex>
    private residueIt: SegmentIterator<ElementIndex>
    private polymerSegment: Segmentation.Segment<ElementIndex>
    private state: AtomicPolymerBackboneIteratorState = AtomicPolymerBackboneIteratorState.nextPolymer
    private pos: SymmetryOperator.CoordinateMapper

    hasNext: boolean = false;

    move() {
        const { residueIt, polymerIt, value, pos } = this

        if (this.state === AtomicPolymerBackboneIteratorState.nextPolymer) {
            if (polymerIt.hasNext) {
                this.polymerSegment = polymerIt.move();
                residueIt.setSegment(this.polymerSegment);
                this.state = AtomicPolymerBackboneIteratorState.firstResidue
            }
        }

        if (this.state === AtomicPolymerBackboneIteratorState.firstResidue) {
            const residueSegment = residueIt.move();
            if (residueIt.hasNext) {
                value.indexB = setTraceElement(value.centerB, residueSegment)
                pos(value.centerB.element, value.posB)
                this.state = AtomicPolymerBackboneIteratorState.nextResidue
            } else {
                this.state = AtomicPolymerBackboneIteratorState.nextPolymer
            }

        }

        if (this.state === AtomicPolymerBackboneIteratorState.nextResidue) {
            const residueSegment = residueIt.move();
            value.centerA.element = value.centerB.element
            value.indexA = value.indexB
            Vec3.copy(value.posA, value.posB)
            value.indexB = setTraceElement(value.centerB, residueSegment)
            pos(value.centerB.element, value.posB)

            if (!residueIt.hasNext) {
                this.state = AtomicPolymerBackboneIteratorState.nextPolymer
            }
        }

        this.hasNext = residueIt.hasNext || polymerIt.hasNext

        return this.value;
    }

    constructor(unit: Unit.Atomic) {
        const { polymerAtomSegments, residueAtomSegments } = unit.model.atomicHierarchy
        this.polymerIt = Segmentation.transientSegments(polymerAtomSegments, unit.elements);
        this.residueIt = Segmentation.transientSegments(residueAtomSegments, unit.elements);
        this.pos = unit.conformation.invariantPosition
        this.value = createPolymerBackbonePair(unit)
        this.hasNext = this.residueIt.hasNext || this.polymerIt.hasNext
    }
}

const enum CoarsePolymerBackboneIteratorState { nextPolymer, firstElement, nextElement }

export class CoarsePolymerBackboneIterator<T extends number = number> implements Iterator<PolymerBackbonePair> {
    private value: PolymerBackbonePair

    private polymerIt: SegmentIterator<ElementIndex>
    private polymerSegment: Segmentation.Segment<ElementIndex>
    private state: CoarsePolymerBackboneIteratorState = CoarsePolymerBackboneIteratorState.nextPolymer
    private pos: SymmetryOperator.CoordinateMapper
    private elementIndex: number

    hasNext: boolean = false;

    move() {
        const { polymerIt, value, pos } = this

        if (this.state === CoarsePolymerBackboneIteratorState.nextPolymer) {
            if (polymerIt.hasNext) {
                this.polymerSegment = polymerIt.move();
                this.elementIndex = this.polymerSegment.start
                this.state = CoarsePolymerBackboneIteratorState.firstElement
            }
        }

        if (this.state === CoarsePolymerBackboneIteratorState.firstElement) {
            this.elementIndex += 1
            if (this.elementIndex + 1 < this.polymerSegment.end) {
                value.centerB.element = value.centerB.unit.elements[this.elementIndex]
                value.indexB = this.elementIndex
                pos(value.centerB.element, value.posB)

                this.state = CoarsePolymerBackboneIteratorState.nextElement
            } else {
                this.state = CoarsePolymerBackboneIteratorState.nextPolymer
            }

        }

        if (this.state === CoarsePolymerBackboneIteratorState.nextElement) {
            this.elementIndex += 1
            value.centerA.element = value.centerB.element
            value.indexA = value.indexB
            Vec3.copy(value.posA, value.posB)

            value.centerB.element = value.centerB.unit.elements[this.elementIndex]
            value.indexB = this.elementIndex
            pos(value.centerB.element, value.posB)

            if (this.elementIndex + 1 >= this.polymerSegment.end) {
                this.state = CoarsePolymerBackboneIteratorState.nextPolymer
            }
        }

        this.hasNext = this.elementIndex + 1 < this.polymerSegment.end || polymerIt.hasNext

        return this.value;
    }

    constructor(unit: Unit.Spheres | Unit.Gaussians) {
        const { polymerElementSegments } = Unit.isSpheres(unit)
            ? unit.model.coarseHierarchy.spheres
            : unit.model.coarseHierarchy.gaussians
        this.polymerIt = Segmentation.transientSegments(polymerElementSegments, unit.elements);

        this.pos = unit.conformation.invariantPosition
        this.value = createPolymerBackbonePair(unit)

        this.hasNext = this.polymerIt.hasNext
    }
}




/**
 * Iterates over individual residues/coarse elements in polymers while providing information
 * about the neighbourhood in the underlying model for drawing splines
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
    index: number
    pos: Vec3
    posPrev: Vec3
    posNext: Vec3
    posNextNext: Vec3
}

function createPolymerTraceElement (unit: Unit) {
    return {
        center: StructureElement.create(unit),
        index: 0,
        pos: Vec3.zero(),
        posPrev: Vec3.zero(),
        posNext: Vec3.zero(),
        posNextNext: Vec3.zero()
    }
}

// const enum AtomicPolymerTraceIteratorState { nextPolymer, firstResidue, nextResidue }

export class AtomicPolymerTraceIterator<T extends number = number> implements Iterator<PolymerTraceElement> {
    private value: PolymerTraceElement

    private polymerIt: SegmentIterator<ElementIndex>
    private residueIt: SegmentIterator<ElementIndex>
    // private polymerSegment: Segmentation.Segment<Element>
    // private state: AtomicPolymerTraceIteratorState = AtomicPolymerTraceIteratorState.nextPolymer
    // private pos: SymmetryOperator.CoordinateMapper

    hasNext: boolean = false;

    move() {
        // const { residueIt, polymerIt, value, pos } = this

        // if (this.state === AtomicPolymerTraceIteratorState.nextPolymer) {
        //     if (polymerIt.hasNext) {
        //         this.polymerSegment = polymerIt.move();
        //         residueIt.setSegment(this.polymerSegment);
        //         this.state = AtomicPolymerTraceIteratorState.firstResidue
        //     }
        // }

        // if (this.state === AtomicPolymerTraceIteratorState.firstResidue) {
        //     const residueSegment = residueIt.move();
        //     if (residueIt.hasNext) {
        //         value.indexB = setTraceElement(value.centerB, residueSegment)
        //         pos(value.centerB.element, value.posB)
        //         this.state = AtomicPolymerTraceIteratorState.nextResidue
        //     } else {
        //         this.state = AtomicPolymerTraceIteratorState.nextPolymer
        //     }

        // }

        // if (this.state === AtomicPolymerTraceIteratorState.nextResidue) {
        //     const residueSegment = residueIt.move();
        //     value.centerA.element = value.centerB.element
        //     value.indexA = value.indexB
        //     Vec3.copy(value.posA, value.posB)
        //     value.indexB = setTraceElement(value.centerB, residueSegment)
        //     pos(value.centerB.element, value.posB)

        //     if (!residueIt.hasNext) {
        //         this.state = AtomicPolymerTraceIteratorState.nextPolymer
        //     }
        // }

        // this.hasNext = residueIt.hasNext || polymerIt.hasNext

        return this.value;
    }

    constructor(unit: Unit.Atomic) {
        const { polymerAtomSegments, residueAtomSegments } = unit.model.atomicHierarchy
        this.polymerIt = Segmentation.transientSegments(polymerAtomSegments, unit.elements);
        this.residueIt = Segmentation.transientSegments(residueAtomSegments, unit.elements);
        // this.pos = unit.conformation.invariantPosition
        this.value = createPolymerTraceElement(unit)
        this.hasNext = this.residueIt.hasNext || this.polymerIt.hasNext
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