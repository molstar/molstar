/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, StructureProperties, Model, ElementIndex } from 'mol-model/structure';
import { Segmentation, OrderedSet, Interval } from 'mol-data/int';
import { MoleculeType, SecondaryStructureType } from 'mol-model/structure/model/types';
import Iterator from 'mol-data/iterator';
import { Vec3 } from 'mol-math/linear-algebra';
import { SymmetryOperator } from 'mol-math/geometry';
import SortedRanges from 'mol-data/int/sorted-ranges';

export function getPolymerRanges(unit: Unit) {
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
                // TODO add OrderedSet.intersectionSize
                count += OrderedSet.size(OrderedSet.intersect(Interval.ofBounds(elements[start], elements[end - 1]), elements))
            }
            break
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

function setTraceElement(l: StructureElement, residueSegment: Segmentation.Segment) {
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
    } else if (moleculeType === MoleculeType.DNA) {
        // traceName = 'P'
        traceName = 'C3\''
    } else if (moleculeType === MoleculeType.RNA) {
        // traceName = 'P'
        traceName = 'C4\''
    }
    return traceName
}

// TODO fix type
function getTraceElement2(model: Model, residueModelSegment: Segmentation.Segment<number>) {
    const traceName = getTraceName2(model, residueModelSegment.index)

    for (let j = residueModelSegment.start, _j = residueModelSegment.end; j < _j; j++) {
        if (model.atomicHierarchy.atoms.label_atom_id.value(j) === traceName) return j
    }
    console.log(`trace name element "${traceName}" not found`, { ...residueModelSegment })
    return residueModelSegment.start
}

function getDirectionName2(model: Model, residueModelIndex: number) {
    const compId = model.atomicHierarchy.residues.label_comp_id.value(residueModelIndex)
    const chemCompMap = model.properties.chemicalComponentMap
    const cc = chemCompMap.get(compId)
    const moleculeType = cc ? cc.moleculeType : MoleculeType.unknown
    let traceName = ''
    if (moleculeType === MoleculeType.protein) {
        traceName = 'O'
    } else if (moleculeType === MoleculeType.DNA) {
        traceName = 'O4\''
    } else if (moleculeType === MoleculeType.RNA) {
        traceName = 'C3\''
    } else {
        console.log('molecule type unknown', Number.isInteger(residueModelIndex) ? compId : residueModelIndex, chemCompMap)
    }
    return traceName
}

// function residueLabel(model: Model, rI: number) {
//     const { residues, chains, residueSegments, chainSegments } = model.atomicHierarchy
//     const { label_comp_id, label_seq_id } = residues
//     const { label_asym_id } = chains
//     const cI = chainSegments.segmentMap[residueSegments.segments[rI]]
//     return `${label_asym_id.value(cI)} ${label_comp_id.value(rI)} ${label_seq_id.value(rI)}`
// }

// TODO fix type
function getDirectionElement2(model: Model, residueModelSegment: Segmentation.Segment<number>) {
    const traceName = getDirectionName2(model, residueModelSegment.index)

    for (let j = residueModelSegment.start, _j = residueModelSegment.end; j < _j; j++) {
        if (model.atomicHierarchy.atoms.label_atom_id.value(j) === traceName) return j
    }
    // console.log('direction name element not found', { ...residueModelSegment })
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

    private polymerIt: SortedRanges.Iterator<ElementIndex>
    private residueIt: Segmentation.SegmentIterator<number> // TODO specific type
    private polymerSegment: Segmentation.Segment<ElementIndex>
    private state: AtomicPolymerBackboneIteratorState = AtomicPolymerBackboneIteratorState.nextPolymer
    private pos: SymmetryOperator.CoordinateMapper

    hasNext: boolean = false;

    move() {
        const { residueIt, polymerIt, value, pos } = this

        if (this.state === AtomicPolymerBackboneIteratorState.nextPolymer) {
            while (polymerIt.hasNext) {
                this.polymerSegment = polymerIt.move();
                // console.log('polymerSegment', this.polymerSegment)
                residueIt.setSegment(this.polymerSegment);

                const residueSegment = residueIt.move();
                // console.log('first residueSegment', residueSegment, residueIt.hasNext)
                if (residueIt.hasNext) {
                    value.indexB = setTraceElement(value.centerB, residueSegment)
                    pos(value.centerB.element, value.posB)
                    this.state = AtomicPolymerBackboneIteratorState.nextResidue
                    break
                }
            }
        }

        if (this.state === AtomicPolymerBackboneIteratorState.nextResidue) {
            const residueSegment = residueIt.move();
            // console.log('next residueSegment', residueSegment)
            value.centerA.element = value.centerB.element
            value.indexA = value.indexB
            Vec3.copy(value.posA, value.posB)
            value.indexB = setTraceElement(value.centerB, residueSegment)
            pos(value.centerB.element, value.posB)

            if (!residueIt.hasNext) {
                // TODO need to advance to a polymer that has two or more residues (can't assume it has)
                this.state = AtomicPolymerBackboneIteratorState.nextPolymer
            }
        }

        this.hasNext = residueIt.hasNext || polymerIt.hasNext
        
        // console.log('hasNext', this.hasNext)
        // console.log('value', this.value)

        return this.value;
    }

    constructor(unit: Unit.Atomic) {
        const { residueAtomSegments } = unit.model.atomicHierarchy
        // console.log('unit.elements', OrderedSet.toArray(unit.elements))
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements)
        this.residueIt = Segmentation.transientSegments(residueAtomSegments, unit.elements)
        this.pos = unit.conformation.invariantPosition
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
    private pos: SymmetryOperator.CoordinateMapper
    private elementIndex: number

    hasNext: boolean = false;

    move() {
        const { polymerIt, value, pos } = this

        if (this.state === CoarsePolymerBackboneIteratorState.nextPolymer) {
            if (polymerIt.hasNext) {
                this.polymerSegment = polymerIt.move();
                console.log('polymer', this.polymerSegment)
                this.elementIndex = this.polymerSegment.start
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
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);

        this.pos = unit.conformation.invariantPosition
        this.value = createPolymerBackbonePair(unit)

        console.log('CoarsePolymerBackboneIterator', this.polymerIt.hasNext)

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
    index: number
    first: boolean
    last: boolean
    secStrucType: SecondaryStructureType

    t0: Vec3
    t1: Vec3
    t2: Vec3
    t3: Vec3
    t4: Vec3

    d12: Vec3
    d23: Vec3
}

function createPolymerTraceElement (unit: Unit): PolymerTraceElement {
    return {
        center: StructureElement.create(unit),
        index: 0,
        first: false,
        last: false,
        secStrucType: SecondaryStructureType.create(SecondaryStructureType.Flag.NA),

        t0: Vec3.zero(),
        t1: Vec3.zero(),
        t2: Vec3.zero(),
        t3: Vec3.zero(),
        t4: Vec3.zero(),

        d12: Vec3.zero(),
        d23: Vec3.zero(),
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
    private residueIt: Segmentation.SegmentIterator<number> // TODO specialize type
    private residueSegmentMin: number
    private residueSegmentMax: number
    private state: AtomicPolymerTraceIteratorState = AtomicPolymerTraceIteratorState.nextPolymer
    private residueSegments: Segmentation<ElementIndex>

    private tmpSegment: Segmentation.Segment<number>

    private unit: Unit.Atomic

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
        this.residueSegmentMin = residueAtomSegments.index[sMin]
        this.residueSegmentMax = residueAtomSegments.index[sMax]
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
            const { tmpSegment, residueSegments, residueSegmentMin, residueSegmentMax } = this
            const residueSegment = residueIt.move();
            const resSegIdx = residueSegment.index
            // console.log(residueLabel(this.unit.model, resSegIdx), resSegIdx, this.unit.model.properties.secondaryStructure.type[resSegIdx])
            value.index = setTraceElement(value.center, residueSegment)

            setSegment(tmpSegment, resSegIdx - 2, residueSegments, residueSegmentMin, residueSegmentMax)
            this.pos(value.t0, getTraceElement2(this.unit.model, tmpSegment))

            setSegment(tmpSegment, resSegIdx - 1, residueSegments, residueSegmentMin, residueSegmentMax)
            this.pos(value.t1, getTraceElement2(this.unit.model, tmpSegment))
            this.pos(value.d12, getDirectionElement2(this.unit.model, tmpSegment))

            setSegment(tmpSegment, resSegIdx, residueSegments, residueSegmentMin, residueSegmentMax)
            value.secStrucType = this.unit.model.properties.secondaryStructure.type[resSegIdx]
            this.pos(value.t2, getTraceElement2(this.unit.model, tmpSegment))
            this.pos(value.d23, getDirectionElement2(this.unit.model, tmpSegment))

            setSegment(tmpSegment, resSegIdx + 1, residueSegments, residueSegmentMin, residueSegmentMax)
            this.pos(value.t3, getTraceElement2(this.unit.model, tmpSegment))

            setSegment(tmpSegment, resSegIdx + 2, residueSegments, residueSegmentMin, residueSegmentMax)
            this.pos(value.t4, getTraceElement2(this.unit.model, tmpSegment))

            value.first = resSegIdx === residueSegmentMin
            value.last = resSegIdx === residueSegmentMax

            if (!residueIt.hasNext) {
                this.state = AtomicPolymerTraceIteratorState.nextPolymer
            }
        }

        this.hasNext = residueIt.hasNext || polymerIt.hasNext

        return this.value;
    }

    constructor(unit: Unit.Atomic) {
        const { residueAtomSegments } = unit.model.atomicHierarchy
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements)
        this.residueIt = Segmentation.transientSegments(residueAtomSegments, unit.elements);
        this.residueSegments = residueAtomSegments
        this.value = createPolymerTraceElement(unit)
        this.hasNext = this.residueIt.hasNext && this.polymerIt.hasNext

        this.tmpSegment = { index: 0, start: 0 as ElementIndex, end: 0 as ElementIndex }

        this.unit = unit
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