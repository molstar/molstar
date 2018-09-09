/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, ElementIndex } from 'mol-model/structure';
import { Segmentation, SortedArray } from 'mol-data/int';
import { isNucleic, MoleculeType } from 'mol-model/structure/model/types';
import { getElementIndexForAtomRole } from 'mol-model/structure/util';

export function getNucleotideElements(unit: Unit.Atomic) {
    const indices: ElementIndex[] = []
    const { elements, model } = unit
    const { chemicalComponentMap } = model.properties
    const { chainAtomSegments, residueAtomSegments, residues } = model.atomicHierarchy
    const { label_comp_id } = residues
    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements)
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements)
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index } = residueIt.move();
            const cc = chemicalComponentMap.get(label_comp_id.value(index))
            const moleculeType = cc ? cc.moleculeType : MoleculeType.unknown

            if (isNucleic(moleculeType)) {
                const elementIndex = getElementIndexForAtomRole(model, index, 'trace')
                indices.push(elementIndex === -1 ? residueAtomSegments.offsets[index] : elementIndex)
            }
        }
    }
    return SortedArray.ofSortedArray<ElementIndex>(indices)
}