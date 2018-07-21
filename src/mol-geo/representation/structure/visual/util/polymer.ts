/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, ElementIndex } from 'mol-model/structure';
import { Segmentation, OrderedSet, Interval } from 'mol-data/int';
import SortedRanges from 'mol-data/int/sorted-ranges';

export * from './polymer/backbone-iterator'
export * from './polymer/gap-iterator'
export * from './polymer/trace-iterator'
export * from './polymer/interpolate'

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
