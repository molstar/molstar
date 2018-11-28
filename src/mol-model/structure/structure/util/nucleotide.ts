/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, ElementIndex } from 'mol-model/structure';
import { Segmentation, SortedArray } from 'mol-data/int';
import { isNucleic } from 'mol-model/structure/model/types';

export function getNucleotideElements(unit: Unit.Atomic) {
    const indices: ElementIndex[] = []
    const { elements, model } = unit
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy
    const { moleculeType, traceElementIndex } = model.atomicHierarchy.derived.residue
    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements)
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements)
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index } = residueIt.move();

            if (isNucleic(moleculeType[index])) {
                const elementIndex = traceElementIndex[index]
                indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex)
            }
        }
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices)
}