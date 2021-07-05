/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
import { Segmentation } from '../../../../../mol-data/int';
import { SortedRanges } from '../../../../../mol-data/int/sorted-ranges';
import { ElementIndex, ResidueIndex, Unit } from '../../../../../mol-model/structure';
import { MoleculeType } from '../../../../../mol-model/structure/model/types';
import { getPolymerRanges } from '../polymer';

export type PolymerBackboneLinkCallback = (indexA: ElementIndex, indexB: ElementIndex, groupA: number, groupB: number, moleculeType: MoleculeType) => void

export function eachPolymerBackboneLink(unit: Unit, callback: PolymerBackboneLinkCallback) {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return eachAtomicPolymerBackboneLink(unit, callback);
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return eachCoarsePolymerBackboneLink(unit, callback);
    }
}

function eachAtomicPolymerBackboneLink(unit: Unit.Atomic, callback: PolymerBackboneLinkCallback) {
    const cyclicPolymerMap = unit.model.atomicRanges.cyclicPolymerMap;
    const polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);
    const residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
    const traceElementIndex = unit.model.atomicHierarchy.derived.residue.traceElementIndex as ArrayLike<ElementIndex>; // can assume it won't be -1 for polymer residues
    const { moleculeType } = unit.model.atomicHierarchy.derived.residue;

    let indexA = -1 as ResidueIndex;
    let indexB = -1 as ResidueIndex;
    let isFirst = true;
    let firstGroup = -1;
    let i = 0;
    while (polymerIt.hasNext) {
        isFirst = true;
        firstGroup = i;
        residueIt.setSegment(polymerIt.move());
        while (residueIt.hasNext) {
            if (isFirst) {
                const index_1 = residueIt.move().index;
                ++i;
                if (!residueIt.hasNext)
                    continue;
                isFirst = false;
                indexB = index_1;
            }
            const index = residueIt.move().index;
            indexA = indexB;
            indexB = index;
            callback(traceElementIndex[indexA], traceElementIndex[indexB], i - 1, i, moleculeType[indexA]);
            ++i;
        }
        if (cyclicPolymerMap.has(indexB)) {
            indexA = indexB;
            indexB = cyclicPolymerMap.get(indexA)!;
            callback(traceElementIndex[indexA], traceElementIndex[indexB], i - 1, firstGroup, moleculeType[indexA]);
        }
    }
}

function eachCoarsePolymerBackboneLink(unit: Unit.Spheres | Unit.Gaussians, callback: PolymerBackboneLinkCallback) {
    const polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);
    const { elements } = unit;

    let isFirst = true;
    let i = 0;
    while (polymerIt.hasNext) {
        isFirst = true;
        const _a = polymerIt.move(), start = _a.start, end = _a.end;
        for (let j = start, jl = end; j < jl; ++j) {
            if (isFirst) {
                ++j;
                ++i;
                if (j > jl)
                    continue;
                isFirst = false;
            }
            callback(elements[j - 1], elements[j], i - 1, i, 0 /* Unknown */);
            ++i;
        }
    }
}

//

export type PolymerBackboneElementCallback = (index: ElementIndex, group: number) => void

export function eachPolymerBackboneElement(unit: Unit, callback: PolymerBackboneElementCallback) {
    switch (unit.kind) {
        case Unit.Kind.Atomic: return eachAtomicPolymerBackboneElement(unit, callback);
        case Unit.Kind.Spheres:
        case Unit.Kind.Gaussians:
            return eachCoarsePolymerBackboneElement(unit, callback);
    }
}

export function eachAtomicPolymerBackboneElement(unit: Unit.Atomic, callback: PolymerBackboneElementCallback) {
    const polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);
    const residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
    const traceElementIndex = unit.model.atomicHierarchy.derived.residue.traceElementIndex as ArrayLike<ElementIndex>; // can assume it won't be -1 for polymer residues

    let i = 0;
    while (polymerIt.hasNext) {
        residueIt.setSegment(polymerIt.move());
        while (residueIt.hasNext) {
            const index = residueIt.move().index;
            callback(traceElementIndex[index], i);
            ++i;
        }
    }
}

function eachCoarsePolymerBackboneElement(unit: Unit.Spheres | Unit.Gaussians, callback: PolymerBackboneElementCallback) {
    const polymerIt = SortedRanges.transientSegments(getPolymerRanges(unit), unit.elements);
    const { elements } = unit;

    let i = 0;
    while (polymerIt.hasNext) {
        const _a = polymerIt.move(), start = _a.start, end = _a.end;
        for (let j = start, jl = end; j < jl; ++j) {
            callback(elements[j], i);
            ++i;
        }
    }
}