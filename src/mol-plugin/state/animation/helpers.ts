/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */


import { SymmetryOperator } from 'mol-math/geometry';
import { Mat4 } from 'mol-math/linear-algebra';
import { Structure } from 'mol-model/structure';
import { StructureUnitTransforms } from 'mol-model/structure/structure/util/unit-transforms';

const _unwindMatrix = Mat4.zero();
export function unwindStructureAssebmly(structure: Structure, unitTransforms: StructureUnitTransforms, t: number) {
    for (let i = 0, _i = structure.units.length; i < _i; i++) {
        const u = structure.units[i];
        SymmetryOperator.lerpFromIdentity(_unwindMatrix, u.conformation.operator, t);
        unitTransforms.setTransform(_unwindMatrix, u);
    }
}