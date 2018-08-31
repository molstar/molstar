/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */
import { Unit, ElementIndex, StructureElement } from 'mol-model/structure';
import { LocationIterator } from '../../../../util/location-iterator';
import { Segmentation } from 'mol-data/int';
import { isNucleic, MoleculeType } from 'mol-model/structure/model/types';
import { getElementIndexForAtomRole } from 'mol-model/structure/util';

export function getNucleotideElementIndices(unit: Unit) {
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
    return indices
}

export namespace NucleotideLocationIterator {
    export function fromGroup(group: Unit.SymmetryGroup): LocationIterator {
        const nucleotideElementIndices = getNucleotideElementIndices(group.units[0])
        const groupCount = nucleotideElementIndices.length
        const instanceCount = group.units.length
        const location = StructureElement.create()
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex]
            location.unit = unit
            location.element = nucleotideElementIndices[groupIndex]
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }
}