/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Element, StructureProperties, Model } from 'mol-model/structure';
import { Segmentation } from 'mol-data/int';
import { MoleculeType } from 'mol-model/structure/model/types';
import Iterator from 'mol-data/iterator';
import { SegmentIterator } from 'mol-data/int/impl/segmentation';
import { Vec3 } from 'mol-math/linear-algebra';
import { SymmetryOperator } from 'mol-math/geometry';

// type TraceMap = Map<number, number>

// interface TraceMaps {
//     atomic: TraceMap
//     spheres: TraceMap
//     gaussians: TraceMap
// }

// function calculateTraceMaps (model: Model): TraceMaps {

// }

export function getPolymerElementCount(unit: Unit) {
    let count = 0
    const { elements } = unit
    const l = Element.Location(unit)
    if (Unit.isAtomic(unit)) {
        const { polymerSegments, residueSegments } = unit.model.atomicHierarchy
        const polymerIt = Segmentation.transientSegments(polymerSegments, elements);
        const residuesIt = Segmentation.transientSegments(residueSegments, elements);
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

function getTraceName(l: Element.Location) {
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

function setTraceElement(l: Element.Location, residueSegment: Segmentation.Segment<Element>) {
    const elements = l.unit.elements
    l.element = elements[residueSegment.start]
    const traceName = getTraceName(l)

    for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
        l.element = elements[j];
        if (StructureProperties.atom.label_atom_id(l) === traceName) return j
    }
    return residueSegment.end - 1
}


function getTraceName2(model: Model, residueModelIndex: number) {
    const compId = model.atomicHierarchy.residues.label_comp_id.value(residueModelIndex)
    const chemCompMap = model.properties.chemicalComponentMap
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

function getTraceElement2(model: Model, residueModelSegment: Segmentation.Segment<Element>) {
    const traceName = getTraceName2(model, residueModelSegment.index)

    for (let j = residueModelSegment.start, _j = residueModelSegment.end; j < _j; j++) {
        if (model.atomicHierarchy.atoms.label_atom_id.value(j) === traceName) return j
    }
    console.log('trace name element not found', { ...residueModelSegment })
    return residueModelSegment.start
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
    centerA: Element.Location
    centerB: Element.Location
    indexA: number
    indexB: number
    posA: Vec3
    posB: Vec3
}

function createPolymerBackbonePair (unit: Unit) {
    return {
        centerA: Element.Location(unit),
        centerB: Element.Location(unit),
        indexA: 0,
        indexB: 0,
        posA: Vec3.zero(),
        posB: Vec3.zero()
    }
}

const enum AtomicPolymerBackboneIteratorState { nextPolymer, firstResidue, nextResidue }

export class AtomicPolymerBackboneIterator<T extends number = number> implements Iterator<PolymerBackbonePair> {
    private value: PolymerBackbonePair

    private polymerIt: SegmentIterator<Element>
    private residueIt: SegmentIterator<Element>
    private polymerSegment: Segmentation.Segment<Element>
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
        const { polymerSegments, residueSegments } = unit.model.atomicHierarchy
        this.polymerIt = Segmentation.transientSegments(polymerSegments, unit.elements);
        this.residueIt = Segmentation.transientSegments(residueSegments, unit.elements);
        this.pos = unit.conformation.invariantPosition
        this.value = createPolymerBackbonePair(unit)
        this.hasNext = this.residueIt.hasNext || this.polymerIt.hasNext
    }
}

const enum CoarsePolymerBackboneIteratorState { nextPolymer, firstElement, nextElement }

export class CoarsePolymerBackboneIterator<T extends number = number> implements Iterator<PolymerBackbonePair> {
    private value: PolymerBackbonePair

    private polymerIt: SegmentIterator<Element>
    private polymerSegment: Segmentation.Segment<Element>
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
        const { polymerSegments } = Unit.isSpheres(unit)
            ? unit.model.coarseHierarchy.spheres
            : unit.model.coarseHierarchy.gaussians
        this.polymerIt = Segmentation.transientSegments(polymerSegments, unit.elements);

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
    center: Element.Location
    index: number
    first: boolean
    last: boolean
    c0: Vec3
    c1: Vec3
    c2: Vec3
    c3: Vec3
    c4: Vec3
}

function createPolymerTraceElement (unit: Unit) {
    return {
        center: Element.Location(unit),
        index: 0,
        first: false,
        last: false,
        c0: Vec3.zero(),
        c1: Vec3.zero(),
        c2: Vec3.zero(),
        c3: Vec3.zero(),
        c4: Vec3.zero()
    }
}

const enum AtomicPolymerTraceIteratorState { nextPolymer, nextResidue }

function setSegment (outSegment: Segmentation.Segment<Element>, index: number, segments: Segmentation<Element>, boundingSegment: Segmentation.Segment<Element>): Segmentation.Segment<Element> {
    index = Math.min(Math.max(0, index), segments.segments.length - 2)
    outSegment.index = index
    outSegment.start = segments.segments[index]
    outSegment.end = segments.segments[index + 1]
    return outSegment
}

// const p0 = Vec3.zero()
// const p1 = Vec3.zero()
// const p2 = Vec3.zero()
// const p3 = Vec3.zero()
// const p4 = Vec3.zero()
// const p5 = Vec3.zero()
// const p6 = Vec3.zero()

export class AtomicPolymerTraceIterator<T extends number = number> implements Iterator<PolymerTraceElement> {
    private value: PolymerTraceElement

    private polymerIt: SegmentIterator<Element>
    private residueIt: SegmentIterator<Element>
    private polymerSegment: Segmentation.Segment<Element>
    private state: AtomicPolymerTraceIteratorState = AtomicPolymerTraceIteratorState.nextPolymer
    private residueSegments: Segmentation<Element>

    private tmpSegment: Segmentation.Segment<Element>

    private unit: Unit.Atomic

    hasNext: boolean = false;

    private pos(target: Vec3, index: number) {
        target[0] = this.unit.model.atomicConformation.x[index]
        target[1] = this.unit.model.atomicConformation.y[index]
        target[2] = this.unit.model.atomicConformation.z[index]
    }

    move() {
        const { residueIt, polymerIt, value } = this
        value.first = false
        value.last = false

        if (this.state === AtomicPolymerTraceIteratorState.nextPolymer) {
            if (polymerIt.hasNext) {
                this.polymerSegment = polymerIt.move();
                residueIt.setSegment(this.polymerSegment);
                this.state = AtomicPolymerTraceIteratorState.nextResidue
                value.first = true
            }
        }

        if (this.state === AtomicPolymerTraceIteratorState.nextResidue) {
            const residueSegment = residueIt.move();
            value.index = setTraceElement(value.center, residueSegment)

            setSegment(this.tmpSegment, residueSegment.index - 2, this.residueSegments, this.polymerSegment)
            this.pos(value.c0, getTraceElement2(this.unit.model, this.tmpSegment))

            setSegment(this.tmpSegment, residueSegment.index - 1, this.residueSegments, this.polymerSegment)
            this.pos(value.c1, getTraceElement2(this.unit.model, this.tmpSegment))

            setSegment(this.tmpSegment, residueSegment.index, this.residueSegments, this.polymerSegment)
            this.pos(value.c2, getTraceElement2(this.unit.model, this.tmpSegment))

            setSegment(this.tmpSegment, residueSegment.index + 1, this.residueSegments, this.polymerSegment)
            this.pos(value.c3, getTraceElement2(this.unit.model, this.tmpSegment))

            setSegment(this.tmpSegment, residueSegment.index + 2, this.residueSegments, this.polymerSegment)
            this.pos(value.c4, getTraceElement2(this.unit.model, this.tmpSegment))

            if (!residueIt.hasNext) {
                this.state = AtomicPolymerTraceIteratorState.nextPolymer
                value.last = true
            }
        }

        this.hasNext = residueIt.hasNext || polymerIt.hasNext

        return this.value;
    }

    constructor(unit: Unit.Atomic) {
        const { polymerSegments, residueSegments } = unit.model.atomicHierarchy
        this.polymerIt = Segmentation.transientSegments(polymerSegments, unit.elements);
        this.residueIt = Segmentation.transientSegments(residueSegments, unit.elements);
        this.residueSegments = residueSegments
        this.value = createPolymerTraceElement(unit)
        this.hasNext = this.residueIt.hasNext || this.polymerIt.hasNext

        this.tmpSegment = { index: 0, start: 0 as Element, end: 0 as Element }

        this.unit = unit
        console.log('model', unit.model)
        console.log('unit', unit)
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