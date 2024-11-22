/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../../mol-data/int/ordered-set';
import { ElementIndex } from '../model/indexing';
import { Structure } from './structure';
import { Unit } from './unit';

/** Serial index of an element in the structure across all units */
export type SerialIndex = { readonly '@type': 'serial-index' } & number

export interface SerialMapping {
    /** Cumulative count of preceding elements for each unit */
    cumulativeUnitElementCount: ArrayLike<number>
    /** Unit index for each serial element in the structure */
    unitIndices: ArrayLike<number>
    /** Element index for each serial element in the structure */
    elementIndices: ArrayLike<ElementIndex>
    /** Get serial index of element in the structure */
    getSerialIndex: (unit: Unit, element: ElementIndex) => SerialIndex
}

export function getSerialMapping(structure: Structure): SerialMapping {
    const { units, elementCount, unitIndexMap } = structure;
    const cumulativeUnitElementCount = new Uint32Array(units.length);
    const unitIndices = new Uint32Array(elementCount);
    const elementIndices = new Uint32Array(elementCount) as unknown as ElementIndex[];
    for (let i = 0, m = 0, il = units.length; i < il; ++i) {
        cumulativeUnitElementCount[i] = m;
        const { elements } = units[i];
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            const mj = m + j;
            unitIndices[mj] = i;
            elementIndices[mj] = elements[j];
        }
        m += elements.length;
    }
    return {
        cumulativeUnitElementCount,
        unitIndices,
        elementIndices,

        getSerialIndex: (unit, element) => cumulativeUnitElementCount[unitIndexMap.get(unit.id)] + OrderedSet.indexOf(unit.elements, element) as SerialIndex
    };
}

//

export interface IntraUnitBondMapping {
    bondCount: number
    unitIndex: ArrayLike<number>
    unitEdgeIndex: ArrayLike<number>
    unitGroupIndex: ArrayLike<number>
    unitGroupOffset: ArrayLike<number>
}

export function getIntraUnitBondMapping(structure: Structure): IntraUnitBondMapping {
    let bondCount = 0;

    const unitGroupOffset: number[] = [];
    for (const ug of structure.unitSymmetryGroups) {
        unitGroupOffset.push(bondCount);
        const unit = ug.units[0];
        if (Unit.isAtomic(unit)) {
            bondCount += unit.bonds.edgeCount * 2 * ug.units.length;
        }
    }

    const unitIndex = new Uint32Array(bondCount);
    const unitEdgeIndex = new Uint32Array(bondCount);
    const unitGroupIndex = new Uint32Array(bondCount);

    let idx = 0;
    let unitIdx = 0;
    for (const ug of structure.unitSymmetryGroups) {
        const unit = ug.units[0];
        if (Unit.isAtomic(unit)) {
            const edgeCount = unit.bonds.edgeCount * 2;
            for (let i = 0, il = ug.units.length; i < il; ++i) {
                for (let j = 0, jl = edgeCount; j < jl; ++j) {
                    unitIndex[idx] = unitIdx;
                    unitEdgeIndex[idx] = j;
                    unitGroupIndex[idx] = i;
                    idx += 1;
                }
            }
        }
        unitIdx += 1;
    }

    return { bondCount, unitIndex, unitEdgeIndex, unitGroupIndex, unitGroupOffset };
}
