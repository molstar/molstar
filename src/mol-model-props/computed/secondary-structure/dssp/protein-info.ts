/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Unit, ResidueIndex, ElementIndex } from '../../../../mol-model/structure';
import { SortedArray } from '../../../../mol-data/int';

export interface ProteinInfo {
    readonly residueIndices: SortedArray<ResidueIndex>
    readonly cIndices: ArrayLike<ElementIndex | -1>
    readonly hIndices: ArrayLike<ElementIndex | -1>
    readonly oIndices: ArrayLike<ElementIndex | -1>
    readonly nIndices: ArrayLike<ElementIndex | -1>
}

export function getUnitProteinInfo(unit: Unit.Atomic): ProteinInfo {
    const { index } = unit.model.atomicHierarchy;
    const { proteinElements, residueIndex } = unit;
    const residueCount = proteinElements.length;

    const unitProteinResidues = new Uint32Array(residueCount);
    const c = new Int32Array(residueCount);
    const h = new Int32Array(residueCount);
    const o = new Int32Array(residueCount);
    const n = new Int32Array(residueCount);

    for (let i = 0; i < residueCount; ++i) {
        const rI = residueIndex[proteinElements[i]];
        unitProteinResidues[i] = rI;
        c[i] = index.findAtomOnResidue(rI, 'C');
        h[i] = index.findAtomOnResidue(rI, 'H');
        o[i] = index.findAtomOnResidue(rI, 'O');
        n[i] = index.findAtomOnResidue(rI, 'N');
    }

    return {
        residueIndices: SortedArray.ofSortedArray(unitProteinResidues),
        cIndices: c as unknown as ArrayLike<ElementIndex | -1>,
        hIndices: h as unknown as ArrayLike<ElementIndex | -1>,
        oIndices: o as unknown as ArrayLike<ElementIndex | -1>,
        nIndices: n as unknown as ArrayLike<ElementIndex | -1>,
    };
}