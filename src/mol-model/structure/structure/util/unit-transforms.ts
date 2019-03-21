/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Unit } from 'mol-model/structure';
import { Mat4 } from 'mol-math/linear-algebra';
import { IntMap } from 'mol-data/int';
import { fillIdentityTransform } from 'mol-geo/geometry/transform-data';

export class StructureUnitTransforms {
    private unitTransforms: Float32Array
    private groupUnitTransforms: Float32Array[] = []
    /** maps unit.id to offset of transform in unitTransforms */
    private unitOffsetMap = IntMap.Mutable<number>();
    private groupIndexMap = IntMap.Mutable<number>();
    private size: number;

    constructor(readonly structure: Structure) {
        this.unitTransforms = new Float32Array(structure.units.length * 16)
        this.size = structure.units.length
        fillIdentityTransform(this.unitTransforms, structure.units.length)
        let groupOffset = 0
        for (let i = 0, il = structure.unitSymmetryGroups.length; i <il; ++i) {
            const g = structure.unitSymmetryGroups[i]
            this.groupIndexMap.set(g.hashCode, i)
            const groupTransforms = this.unitTransforms.subarray(groupOffset, groupOffset + g.units.length * 16)
            this.groupUnitTransforms.push(groupTransforms)
            for (let j = 0, jl = g.units.length; j < jl; ++j) {
                this.unitOffsetMap.set(g.units[j].id, groupOffset + j * 16)
            }
            groupOffset += g.units.length * 16
        }
    }

    reset() {
        fillIdentityTransform(this.unitTransforms, this.size);
    }

    setTransform(matrix: Mat4, unit: Unit) {
        Mat4.toArray(matrix, this.unitTransforms, this.unitOffsetMap.get(unit.id))
    }

    getTransform(out: Mat4, unit: Unit) {
        return Mat4.fromArray(out, this.unitTransforms, this.unitOffsetMap.get(unit.id))
    }

    getSymmetryGroupTransforms(group: Unit.SymmetryGroup): Float32Array {
        return this.groupUnitTransforms[this.groupIndexMap.get(group.hashCode)]
    }
}