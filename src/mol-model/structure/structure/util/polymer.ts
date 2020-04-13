/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, ElementIndex } from '../../../../mol-model/structure';
import { Segmentation, OrderedSet, Interval, SortedArray } from '../../../../mol-data/int';
import SortedRanges from '../../../../mol-data/int/sorted-ranges';
import { isNucleic, isProtein } from '../../../../mol-model/structure/model/types';

export function getAtomicPolymerElements(unit: Unit.Atomic) {
    const indices: ElementIndex[] = [];
    const { elements, model } = unit;
    const { residueAtomSegments } = unit.model.atomicHierarchy;
    const { traceElementIndex } = model.atomicHierarchy.derived.residue;
    const polymerIt = SortedRanges.transientSegments(unit.model.atomicRanges.polymerRanges, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);
    while (polymerIt.hasNext) {
        const polymerSegment = polymerIt.move();
        residueIt.setSegment(polymerSegment);
        while (residueIt.hasNext) {
            const residueSegment = residueIt.move();
            const { start, end, index } = residueSegment;
            if (OrderedSet.areIntersecting(Interval.ofRange(elements[start], elements[end - 1]), elements)) {
                const elementIndex = traceElementIndex[index];
                indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex);
            }
        }
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices);
}

export function getCoarsePolymerElements(unit: Unit.Spheres | Unit.Gaussians) {
    const indices: ElementIndex[] = [];
    const { elements, model } = unit;
    const { spheres, gaussians } = model.coarseHierarchy;
    const polymerRanges = Unit.isSpheres(unit) ? spheres.polymerRanges : gaussians.polymerRanges;
    const polymerIt = SortedRanges.transientSegments(polymerRanges, elements);
    while (polymerIt.hasNext) {
        const { start, end } = polymerIt.move();
        for (let i = start; i < end; ++i) { indices.push(elements[i]); }
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices);
}

//

export function getAtomicGapElements(unit: Unit.Atomic) {
    const indices: ElementIndex[] = [];
    const { elements, model, residueIndex } = unit;
    const { residueAtomSegments } = unit.model.atomicHierarchy;
    const { traceElementIndex } = model.atomicHierarchy.derived.residue;
    const gapIt = SortedRanges.transientSegments(unit.model.atomicRanges.gapRanges, unit.elements);
    while (gapIt.hasNext) {
        const gapSegment = gapIt.move();
        const indexStart = residueIndex[elements[gapSegment.start]];
        const indexEnd = residueIndex[elements[gapSegment.end - 1]];
        const elementIndexStart = traceElementIndex[indexStart];
        const elementIndexEnd = traceElementIndex[indexEnd];
        indices.push(elementIndexStart === -1 ? residueAtomSegments.offsets[indexStart] : elementIndexStart);
        indices.push(elementIndexEnd === -1 ? residueAtomSegments.offsets[indexEnd] : elementIndexEnd);

    }
    return SortedArray.ofSortedArray<ElementIndex>(indices);
}

export function getCoarseGapElements(unit: Unit.Spheres | Unit.Gaussians) {
    const indices: ElementIndex[] = [];
    const { elements, model } = unit;
    const { spheres, gaussians } = model.coarseHierarchy;
    const gapRanges = Unit.isSpheres(unit) ? spheres.gapRanges : gaussians.gapRanges;
    const gapIt = SortedRanges.transientSegments(gapRanges, elements);
    while (gapIt.hasNext) {
        const { start, end } = gapIt.move();
        indices.push(elements[start], elements[end - 1]);
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices);
}

//

export function getNucleotideElements(unit: Unit.Atomic) {
    const indices: ElementIndex[] = [];
    const { elements, model } = unit;
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
    const { moleculeType, traceElementIndex } = model.atomicHierarchy.derived.residue;
    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index } = residueIt.move();

            if (isNucleic(moleculeType[index])) {
                const elementIndex = traceElementIndex[index];
                indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex);
            }
        }
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices);
}

export function getProteinElements(unit: Unit.Atomic) {
    const indices: ElementIndex[] = [];
    const { elements, model } = unit;
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
    const { moleculeType, traceElementIndex } = model.atomicHierarchy.derived.residue;
    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index } = residueIt.move();

            if (isProtein(moleculeType[index])) {
                const elementIndex = traceElementIndex[index];
                indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex);
            }
        }
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices);
}