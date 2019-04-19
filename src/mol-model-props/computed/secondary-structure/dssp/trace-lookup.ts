/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { AtomicHierarchy, AtomicConformation } from 'mol-model/structure/model/properties/atomic';
import { MoleculeType } from 'mol-model/structure/model/types';
import { GridLookup3D } from 'mol-math/geometry';
import { SortedArray } from 'mol-data/int';
import { ResidueIndex } from 'mol-model/structure';

export function calcAtomicTraceLookup3D(hierarchy: AtomicHierarchy, conformation: AtomicConformation) {
    const { x, y, z } = conformation;
    const { moleculeType, traceElementIndex } = hierarchy.derived.residue
    const indices: number[] = []
    const _proteinResidues: number[] = []
    for (let i = 0, il = moleculeType.length; i < il; ++i) {
        if (moleculeType[i] === MoleculeType.protein) {
            indices[indices.length] = traceElementIndex[i]
            _proteinResidues[_proteinResidues.length] = i
        }
    }
    const lookup3d = GridLookup3D({ x, y, z, indices: SortedArray.ofSortedArray(indices) }, 4);
    const proteinResidues = SortedArray.ofSortedArray<ResidueIndex>(_proteinResidues)
    return { lookup3d, proteinResidues }
}