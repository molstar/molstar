/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GridLookup3D } from '@molstar/core/math/geometry';
import { SortedArray } from '@molstar/core/data/int';
import type { Unit } from '@molstar/model/model/structure/structure';
import type { ResidueIndex } from '@molstar/model/model/structure';
import { getBoundary } from '@molstar/core/math/geometry/boundary';

export function calcUnitProteinTraceLookup3D(
  unit: Unit.Atomic,
  unitProteinResidues: SortedArray<ResidueIndex>,
): GridLookup3D {
  const { x, y, z } = unit.model.atomicConformation;
  const { traceElementIndex } = unit.model.atomicHierarchy.derived.residue;
  const indices = new Uint32Array(unitProteinResidues.length);
  for (let i = 0, il = unitProteinResidues.length; i < il; ++i) {
    indices[i] = traceElementIndex[unitProteinResidues[i]];
  }
  const position = { x, y, z, indices: SortedArray.ofSortedArray(indices) };
  return GridLookup3D(position, getBoundary(position));
}
