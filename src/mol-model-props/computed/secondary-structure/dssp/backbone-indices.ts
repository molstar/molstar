/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { AtomicHierarchy } from 'mol-model/structure/model/properties/atomic';
import { SortedArray } from 'mol-data/int';
import { ResidueIndex, ElementIndex } from 'mol-model/structure';
import { BackboneAtomIndices } from './common';

export function calcBackboneAtomIndices(hierarchy: AtomicHierarchy, proteinResidues: SortedArray<ResidueIndex>): BackboneAtomIndices {
    const residueCount = proteinResidues.length
    const { index } = hierarchy

    const c = new Int32Array(residueCount)
    const h = new Int32Array(residueCount)
    const o = new Int32Array(residueCount)
    const n = new Int32Array(residueCount)

    for (let i = 0, il = residueCount; i < il; ++i) {
        const rI = proteinResidues[i]
        c[i] = index.findAtomOnResidue(rI, 'C')
        h[i] = index.findAtomOnResidue(rI, 'H')
        o[i] = index.findAtomOnResidue(rI, 'O')
        n[i] = index.findAtomOnResidue(rI, 'N')
    }

    return {
        cIndices: c as unknown as ArrayLike<ElementIndex | -1>,
        hIndices: h as unknown as ArrayLike<ElementIndex | -1>,
        oIndices: o as unknown as ArrayLike<ElementIndex | -1>,
        nIndices: n as unknown as ArrayLike<ElementIndex | -1>,
    }
}