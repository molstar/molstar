/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, ElementIndex } from 'mol-model/structure';
import { Segmentation, OrderedSet, Interval, SortedArray } from 'mol-data/int';
import SortedRanges from 'mol-data/int/sorted-ranges';
import { getElementIndexForAtomRole } from 'mol-model/structure/util';

export function getAtomicPolymerElements(unit: Unit.Atomic) {
    const indices: ElementIndex[] = []
    const { elements, model } = unit
    const { residueAtomSegments } = unit.model.atomicHierarchy
    const polymerIt = SortedRanges.transientSegments(unit.model.atomicHierarchy.polymerRanges, elements)
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements)
    while (polymerIt.hasNext) {
        const polymerSegment = polymerIt.move()
        residueIt.setSegment(polymerSegment)
        while (residueIt.hasNext) {
            const residueSegment = residueIt.move()
            const { start, end, index } = residueSegment
            if (OrderedSet.areIntersecting(Interval.ofBounds(elements[start], elements[end - 1]), elements)) {
                const elementIndex = getElementIndexForAtomRole(model, index, 'trace')
                indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex)
            }
        }
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices)
}

export function getCoarsePolymerElements(unit: Unit.Spheres | Unit.Gaussians) {
    const indices: ElementIndex[] = []
    const { elements, model } = unit
    const { spheres, gaussians } = model.coarseHierarchy
    const polymerRanges = Unit.isSpheres(unit) ? spheres.polymerRanges : gaussians.polymerRanges
    const polymerIt = SortedRanges.transientSegments(polymerRanges, elements)
    while (polymerIt.hasNext) {
        const { start, end } = polymerIt.move()
        for (let i = start; i < end; ++i) { indices.push(elements[i]) }
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices)
}

export function getAtomicGapElements(unit: Unit.Atomic) {
    const indices: ElementIndex[] = []
    const { elements, model, residueIndex } = unit
    const { residueAtomSegments } = unit.model.atomicHierarchy
    const gapIt = SortedRanges.transientSegments(unit.model.atomicHierarchy.gapRanges, unit.elements);
    while (gapIt.hasNext) {
        const gapSegment = gapIt.move();
        const indexStart = residueIndex[elements[gapSegment.start]]
        const indexEnd = residueIndex[elements[gapSegment.end - 1]]
        const elementIndexStart = getElementIndexForAtomRole(model, indexStart, 'trace')
        const elementIndexEnd = getElementIndexForAtomRole(model, indexEnd, 'trace')
        indices.push(elementIndexStart === -1 ? residueAtomSegments.offsets[indexStart] : elementIndexStart)
        indices.push(elementIndexEnd === -1 ? residueAtomSegments.offsets[indexEnd] : elementIndexEnd)

    }
    return SortedArray.ofSortedArray<ElementIndex>(indices)
}

export function getCoarseGapElements(unit: Unit.Spheres | Unit.Gaussians) {
    const indices: ElementIndex[] = []
    const { elements, model } = unit
    const { spheres, gaussians } = model.coarseHierarchy
    const gapRanges = Unit.isSpheres(unit) ? spheres.gapRanges : gaussians.gapRanges
    const gapIt = SortedRanges.transientSegments(gapRanges, elements)
    while (gapIt.hasNext) {
        const { start, end } = gapIt.move()
        indices.push(elements[start], elements[end - 1])
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices)
}