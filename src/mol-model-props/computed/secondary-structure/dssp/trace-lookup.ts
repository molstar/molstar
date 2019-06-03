/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GridLookup3D } from '../../../../mol-math/geometry';
import { SortedArray } from '../../../../mol-data/int';
import { Unit } from '../../../../mol-model/structure/structure';
import { ResidueIndex } from '../../../../mol-model/structure';

export function calcUnitProteinTraceLookup3D(unit: Unit.Atomic, unitProteinResidues: SortedArray<ResidueIndex>): GridLookup3D {
    const { x, y, z } = unit.model.atomicConformation;
    const { traceElementIndex } = unit.model.atomicHierarchy.derived.residue
    const indices = new Uint32Array(unitProteinResidues.length)
    for (let i = 0, il = unitProteinResidues.length; i < il; ++i) {
        indices[i] = traceElementIndex[unitProteinResidues[i]]
    }
    return GridLookup3D({ x, y, z, indices: SortedArray.ofSortedArray(indices) });
}