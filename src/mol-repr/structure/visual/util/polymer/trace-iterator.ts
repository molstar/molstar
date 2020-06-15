/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, StructureElement, ElementIndex, ResidueIndex, Structure } from '../../../../../mol-model/structure';
import { Segmentation, SortedArray } from '../../../../../mol-data/int';
import { MoleculeType, SecondaryStructureType } from '../../../../../mol-model/structure/model/types';
import Iterator from '../../../../../mol-data/iterator';
import { Vec3 } from '../../../../../mol-math/linear-algebra';
import SortedRanges from '../../../../../mol-data/int/sorted-ranges';
import { CoarseSphereConformation, CoarseGaussianConformation } from '../../../../../mol-model/structure/model/properties/coarse';
import { getPolymerRanges } from '../polymer';
import { AtomicConformation } from '../../../../../mol-model/structure/model/properties/atomic';
import { SecondaryStructureProvider } from '../../../../../mol-model-props/computed/secondary-structure';
import { SecondaryStructure } from '../../../../../mol-model/structure/model/properties/seconday-structure';

/**
 * Iterates over individual residues/coarse elements in polymers of a unit while
 * providing information about the neighbourhood in the underlying model for drawing splines
 */
export function PolymerTraceIterator(unit: Unit, structure: Structure): Iterator<PolymerTraceElement> {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return new AtomicPolymerTraceIterator(unit, structure);
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return new CoarsePolymerTraceIterator(unit, structure);
    }
}

interface PolymerTraceElement {
    center: StructureElement.Location
    centerPrev: StructureElement.Location
    centerNext: StructureElement.Location
    first: boolean, last: boolean
    initial: boolean, final: boolean
    secStrucFirst: boolean, secStrucLast: boolean
    secStrucType: SecondaryStructureType
    moleculeType: MoleculeType
    isCoarseBackbone: boolean
    coarseBackboneFirst: boolean, coarseBackboneLast: boolean

    p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3
    d12: Vec3, d23: Vec3
}

const SecStrucTypeNA = SecondaryStructureType.create(SecondaryStructureType.Flag.NA);

function createPolymerTraceElement (structure: Structure, unit: Unit): PolymerTraceElement {
    return {
        center: StructureElement.Location.create(structure, unit),
        centerPrev: StructureElement.Location.create(structure, unit),
        centerNext: StructureElement.Location.create(structure, unit),
        first: false, last: false,
        initial: false, final: false,
        secStrucFirst: false, secStrucLast: false,
        secStrucType: SecStrucTypeNA,
        moleculeType: MoleculeType.Unknown,
        coarseBackboneFirst: false, coarseBackboneLast: false,
        isCoarseBackbone: false,
        p0: Vec3(), p1: Vec3(), p2: Vec3(), p3: Vec3(), p4: Vec3(),
        d12: Vec3(), d23: Vec3()
    };
}

const enum AtomicPolymerTraceIteratorState { nextPolymer, nextResidue }

const tmpDir = Vec3();
const tmpVecA = Vec3();
const tmpVecB = Vec3();

export class AtomicPolymerTraceIterator implements Iterator<PolymerTraceElement> {
    private value: PolymerTraceElement
    private polymerIt: SortedRanges.Iterator<ElementIndex, number>
    private residueIt: Segmentation.SegmentIterator<ResidueIndex>
    private polymerSegment: Segmentation.Segment<number>
    private cyclicPolymerMap: Map<ResidueIndex, ResidueIndex>
    private secondaryStructureType: SecondaryStructure['type']
    private secondaryStructureGetIndex: SecondaryStructure['getIndex']
    private residueSegmentMin: ResidueIndex
    private residueSegmentMax: ResidueIndex
    private prevSecStrucType: SecondaryStructureType
    private currSecStrucType: SecondaryStructureType
    private nextSecStrucType: SecondaryStructureType
    private prevCoarseBackbone: boolean
    private currCoarseBackbone: boolean
    private nextCoarseBackbone: boolean
    private state: AtomicPolymerTraceIteratorState = AtomicPolymerTraceIteratorState.nextPolymer
    private polymerRanges: SortedArray<ElementIndex>
    private residueAtomSegments: Segmentation<ElementIndex, ResidueIndex>
    private traceElementIndex: ArrayLike<ElementIndex>
    private directionFromElementIndex: ArrayLike<ElementIndex | -1>
    private directionToElementIndex: ArrayLike<ElementIndex | -1>
    private moleculeType: ArrayLike<MoleculeType>
    private atomicConformation: AtomicConformation

    private p0 = Vec3()
    private p1 = Vec3()
    private p2 = Vec3()
    private p3 = Vec3()
    private p4 = Vec3()
    private p5 = Vec3()
    private p6 = Vec3()

    private d01 = Vec3()
    private d12 = Vec3()
    private d23 = Vec3()
    private d34 = Vec3()

    hasNext: boolean = false;

    private pos(target: Vec3, index: number) {
        if (index !== -1) {
            target[0] = this.atomicConformation.x[index];
            target[1] = this.atomicConformation.y[index];
            target[2] = this.atomicConformation.z[index];
        }
    }

    private updateResidueSegmentRange(polymerSegment: Segmentation.Segment<number>) {
        const { index } = this.residueAtomSegments;
        this.residueSegmentMin = index[this.polymerRanges[polymerSegment.index * 2]];
        this.residueSegmentMax = index[this.polymerRanges[polymerSegment.index * 2 + 1]];
    }

    private getResidueIndex(residueIndex: number) {
        if (residueIndex < this.residueSegmentMin) {
            const cyclicIndex = this.cyclicPolymerMap.get(this.residueSegmentMin);
            if (cyclicIndex !== undefined) {
                residueIndex = cyclicIndex - (this.residueSegmentMin - residueIndex - 1);
            } else {
                residueIndex = this.residueSegmentMin;
            }
        } else if (residueIndex > this.residueSegmentMax) {
            const cyclicIndex = this.cyclicPolymerMap.get(this.residueSegmentMax);
            if (cyclicIndex !== undefined) {
                residueIndex = cyclicIndex + (residueIndex - this.residueSegmentMax - 1);
            } else {
                residueIndex = this.residueSegmentMax;
            }
        }
        return residueIndex as ResidueIndex;
    }

    private getSecStruc(residueIndex: number) {
        return this.secondaryStructureType[this.secondaryStructureGetIndex(residueIndex as ResidueIndex)];
    }

    private setControlPoint(out: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, residueIndex: ResidueIndex) {
        const ss =  this.getSecStruc(residueIndex);
        if (SecondaryStructureType.is(ss, SecondaryStructureType.Flag.Beta)) {
            Vec3.scale(out, Vec3.add(out, p1, Vec3.add(out, p3, Vec3.add(out, p2, p2))), 1 / 4);
        } else {
            Vec3.copy(out, p2);
        }
    }

    private setFromToVector(out: Vec3, residueIndex: ResidueIndex) {
        if (this.value.isCoarseBackbone) {
            Vec3.set(out, 1, 0, 0);
        } else {
            this.pos(tmpVecA, this.directionFromElementIndex[residueIndex]);
            this.pos(tmpVecB, this.directionToElementIndex[residueIndex]);
            Vec3.sub(out, tmpVecB, tmpVecA);
        }
    }

    private setDirection(out: Vec3, v1: Vec3, v2: Vec3, v3: Vec3) {
        Vec3.matchDirection(tmpVecA, v1, v2);
        Vec3.matchDirection(tmpVecB, v3, v2);
        Vec3.scale(out, Vec3.add(out, tmpVecA, Vec3.add(out, tmpVecB, Vec3.add(out, v2, v2))), 1 / 4);
    }

    move() {
        const { residueIt, polymerIt, value } = this;

        if (this.state === AtomicPolymerTraceIteratorState.nextPolymer) {
            while (polymerIt.hasNext) {
                this.polymerSegment = polymerIt.move();
                residueIt.setSegment(this.polymerSegment);
                this.updateResidueSegmentRange(this.polymerSegment);
                if (residueIt.hasNext) {
                    this.state = AtomicPolymerTraceIteratorState.nextResidue;
                    const residueIndexBeg = this.residueAtomSegments.index[this.unit.elements[this.polymerSegment.start]];
                    const residueIndexBegPrev = this.getResidueIndex(residueIndexBeg - 1);
                    this.currSecStrucType = residueIndexBeg === residueIndexBegPrev ? SecStrucTypeNA : this.getSecStruc(residueIndexBegPrev);
                    this.nextSecStrucType = this.getSecStruc(residueIndexBeg);
                    this.currCoarseBackbone = this.directionFromElementIndex[residueIndexBegPrev] === -1 || this.directionToElementIndex[residueIndexBegPrev] === -1;
                    this.nextCoarseBackbone = this.directionFromElementIndex[residueIndexBeg] === -1 || this.directionToElementIndex[residueIndexBeg] === -1;
                    break;
                }
            }
        }

        if (this.state === AtomicPolymerTraceIteratorState.nextResidue) {
            const { index: residueIndex } = residueIt.move();
            const residueIndexPrev3 = this.getResidueIndex(residueIndex - 3);
            const residueIndexPrev2 = this.getResidueIndex(residueIndex - 2);
            const residueIndexPrev1 = this.getResidueIndex(residueIndex - 1);
            const residueIndexNext1 = this.getResidueIndex(residueIndex + 1);
            const residueIndexNext2 = this.getResidueIndex(residueIndex + 2);
            const residueIndexNext3 = this.getResidueIndex(residueIndex + 3);

            this.prevSecStrucType = this.currSecStrucType;
            this.currSecStrucType = this.nextSecStrucType;
            this.nextSecStrucType = residueIndex === residueIndexNext1 ? SecStrucTypeNA : this.getSecStruc(residueIndexNext1);

            this.prevCoarseBackbone = this.currCoarseBackbone;
            this.currCoarseBackbone = this.nextCoarseBackbone;
            this.nextCoarseBackbone = residueIndex === residueIndexNext1 ? false : (this.directionFromElementIndex[residueIndexNext1] === -1 || this.directionToElementIndex[residueIndexNext1] === -1);

            value.secStrucType = this.currSecStrucType;
            value.secStrucFirst = this.prevSecStrucType !== this.currSecStrucType;
            value.secStrucLast = this.currSecStrucType !== this.nextSecStrucType;
            value.isCoarseBackbone = this.currCoarseBackbone;
            value.coarseBackboneFirst = this.prevCoarseBackbone !== this.currCoarseBackbone;
            value.coarseBackboneLast = this.currCoarseBackbone !== this.nextCoarseBackbone;
            value.first = residueIndex === this.residueSegmentMin;
            value.last = residueIndex === this.residueSegmentMax;
            value.moleculeType = this.moleculeType[residueIndex];
            value.isCoarseBackbone = this.directionFromElementIndex[residueIndex] === -1 || this.directionToElementIndex[residueIndex] === -1;

            value.initial = residueIndex === residueIndexPrev1;
            value.final = residueIndex === residueIndexNext1;

            value.centerPrev.element = this.traceElementIndex[residueIndexPrev1];
            value.center.element = this.traceElementIndex[residueIndex];
            value.centerNext.element = this.traceElementIndex[residueIndexNext1];

            this.pos(this.p0, this.traceElementIndex[residueIndexPrev3]);
            this.pos(this.p1, this.traceElementIndex[residueIndexPrev2]);
            this.pos(this.p2, this.traceElementIndex[residueIndexPrev1]);
            this.pos(this.p3, this.traceElementIndex[residueIndex]);
            this.pos(this.p4, this.traceElementIndex[residueIndexNext1]);
            this.pos(this.p5, this.traceElementIndex[residueIndexNext2]);
            this.pos(this.p6, this.traceElementIndex[residueIndexNext3]);

            this.setFromToVector(this.d01, residueIndexPrev1);
            this.setFromToVector(this.d12, residueIndex);
            this.setFromToVector(this.d23, residueIndexNext1);
            this.setFromToVector(this.d34, residueIndexNext2);

            // extend termini
            const f = 0.5;
            if (residueIndex === residueIndexPrev1) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.p3, this.p4), f);
                Vec3.add(this.p2, this.p3, tmpDir);
                Vec3.add(this.p1, this.p2, tmpDir);
                Vec3.add(this.p0, this.p1, tmpDir);
            } else if (residueIndexPrev1 === residueIndexPrev2) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.p2, this.p3), f);
                Vec3.add(this.p1, this.p2, tmpDir);
                Vec3.add(this.p0, this.p1, tmpDir);
            } else if (residueIndexPrev2 === residueIndexPrev3) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.p1, this.p2), f);
                Vec3.add(this.p0, this.p1, tmpDir);
            }
            if (residueIndex === residueIndexNext1) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.p3, this.p2), f);
                Vec3.add(this.p4, this.p3, tmpDir);
                Vec3.add(this.p5, this.p4, tmpDir);
                Vec3.add(this.p6, this.p5, tmpDir);
            } else if (residueIndexNext1 === residueIndexNext2) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.p4, this.p3), f);
                Vec3.add(this.p5, this.p4, tmpDir);
                Vec3.add(this.p6, this.p5, tmpDir);
            } else if (residueIndexNext2 === residueIndexNext3) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.p5, this.p4), f);
                Vec3.add(this.p6, this.p5, tmpDir);
            }

            this.setControlPoint(value.p0, this.p0, this.p1, this.p2, residueIndexPrev2);
            this.setControlPoint(value.p1, this.p1, this.p2, this.p3, residueIndexPrev1);
            this.setControlPoint(value.p2, this.p2, this.p3, this.p4, residueIndex);
            this.setControlPoint(value.p3, this.p3, this.p4, this.p5, residueIndexNext1);
            this.setControlPoint(value.p4, this.p4, this.p5, this.p6, residueIndexNext2);

            this.setDirection(value.d12, this.d01, this.d12, this.d23);
            this.setDirection(value.d23, this.d12, this.d23, this.d34);

            if (!residueIt.hasNext) {
                this.state = AtomicPolymerTraceIteratorState.nextPolymer;
            }
        }

        this.hasNext = residueIt.hasNext || polymerIt.hasNext;

        return this.value;
    }

    constructor(private unit: Unit.Atomic, structure: Structure) {
        this.atomicConformation = unit.model.atomicConformation;
        this.residueAtomSegments = unit.model.atomicHierarchy.residueAtomSegments;
        this.polymerRanges = unit.model.atomicRanges.polymerRanges;
        this.traceElementIndex = unit.model.atomicHierarchy.derived.residue.traceElementIndex as ArrayLike<ElementIndex>; // can assume it won't be -1 for polymer residues
        this.directionFromElementIndex = unit.model.atomicHierarchy.derived.residue.directionFromElementIndex;
        this.directionToElementIndex = unit.model.atomicHierarchy.derived.residue.directionToElementIndex;
        this.moleculeType = unit.model.atomicHierarchy.derived.residue.moleculeType;
        this.cyclicPolymerMap = unit.model.atomicRanges.cyclicPolymerMap;
        this.polymerIt = SortedRanges.transientSegments(this.polymerRanges, unit.elements);
        this.residueIt = Segmentation.transientSegments(this.residueAtomSegments, unit.elements);
        this.value = createPolymerTraceElement(structure, unit);
        this.hasNext = this.residueIt.hasNext && this.polymerIt.hasNext;

        const secondaryStructure = SecondaryStructureProvider.get(structure).value?.get(unit.invariantId);
        if (!secondaryStructure) throw new Error('missing secondary structure');
        this.secondaryStructureType = secondaryStructure.type;
        this.secondaryStructureGetIndex = secondaryStructure.getIndex;
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

    private getElementIndex(elementIndex: number) {
        return Math.min(Math.max(this.polymerSegment.start, elementIndex), this.polymerSegment.end - 1) as ElementIndex;
    }

    private pos(target: Vec3, elementIndex: number) {
        const index = this.unit.elements[elementIndex];
        target[0] = this.conformation.x[index];
        target[1] = this.conformation.y[index];
        target[2] = this.conformation.z[index];
    }

    move() {
        if (this.state === CoarsePolymerTraceIteratorState.nextPolymer) {
            while (this.polymerIt.hasNext) {
                this.polymerSegment = this.polymerIt.move();
                this.elementIndex = this.polymerSegment.start;

                if (this.elementIndex < this.polymerSegment.end) {
                    this.state = CoarsePolymerTraceIteratorState.nextElement;
                    break;
                }
            }
        }

        if (this.state === CoarsePolymerTraceIteratorState.nextElement) {
            const elementIndexPrev2 = this.getElementIndex(this.elementIndex - 2);
            const elementIndexPrev1 = this.getElementIndex(this.elementIndex - 1);
            const elementIndexNext1 = this.getElementIndex(this.elementIndex + 1);
            const elementIndexNext2 = this.getElementIndex(this.elementIndex + 2);

            this.value.centerPrev.element = this.value.center.unit.elements[elementIndexPrev1];
            this.value.center.element = this.value.center.unit.elements[this.elementIndex];
            this.value.centerNext.element = this.value.center.unit.elements[elementIndexNext1];

            this.pos(this.value.p0, elementIndexPrev2);
            this.pos(this.value.p1, elementIndexPrev1);
            this.pos(this.value.p2, this.elementIndex);
            this.pos(this.value.p3, elementIndexNext1);
            this.pos(this.value.p4, elementIndexNext2);

            // extend termini
            const f = 0.5;
            if (this.elementIndex === elementIndexPrev1) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.value.p2, this.value.p3), f);
                Vec3.add(this.value.p1, this.value.p2, tmpDir);
                Vec3.add(this.value.p0, this.value.p1, tmpDir);
            } else if (elementIndexPrev1 === elementIndexPrev2) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.value.p1, this.value.p2), f);
                Vec3.add(this.value.p0, this.value.p1, tmpDir);
            }
            if (this.elementIndex === elementIndexNext1) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.value.p2, this.value.p1), f);
                Vec3.add(this.value.p3, this.value.p2, tmpDir);
                Vec3.add(this.value.p4, this.value.p3, tmpDir);
            } else if (elementIndexNext1 === elementIndexNext2) {
                Vec3.scale(tmpDir, Vec3.sub(tmpDir, this.value.p3, this.value.p2), f);
                Vec3.add(this.value.p4, this.value.p3, tmpDir);
            }

            this.value.first = this.elementIndex === this.polymerSegment.start;
            this.value.last = this.elementIndex === this.polymerSegment.end - 1;

            if (this.elementIndex + 1 >= this.polymerSegment.end) {
                this.state = CoarsePolymerTraceIteratorState.nextPolymer;
            }
        }

        this.hasNext = this.elementIndex + 1 < this.polymerSegment.end || this.polymerIt.hasNext;
        this.elementIndex += 1;
        return this.value;
    }

    constructor(private unit: Unit.Spheres | Unit.Gaussians, structure: Structure) {
        this.polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);
        this.value = createPolymerTraceElement(structure, unit);
        Vec3.set(this.value.d12, 1, 0, 0);
        Vec3.set(this.value.d23, 1, 0, 0);
        switch (unit.kind) {
            case Unit.Kind.Spheres: this.conformation = unit.model.coarseConformation.spheres; break;
            case Unit.Kind.Gaussians: this.conformation = unit.model.coarseConformation.gaussians; break;
        }
        this.hasNext = this.polymerIt.hasNext;
    }
}