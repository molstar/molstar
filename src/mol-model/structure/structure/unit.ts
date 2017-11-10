/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import SymmetryOperator from 'mol-math/geometry/symmetry-operator'
import { Model } from '../model'

interface Unit extends SymmetryOperator.ArrayMapping {
    // Structure-level unique identifier of the unit.
    readonly id: number,

    // Provides access to the underlying data.
    readonly model: Model,

    // Determines the operation applied to this unit.
    // The transform and and inverse are baked into the "getPosition" function
    readonly operator: SymmetryOperator,

    // Reference some commonly accessed things for faster access.
    readonly residueIndex: ArrayLike<number>,
    readonly chainIndex: ArrayLike<number>,
    readonly hierarchy: Model['hierarchy'],
    readonly conformation: Model['conformation']

    // TODO: add velocity?
}

namespace Unit {
    export function create(id: number, model: Model, operator: SymmetryOperator): Unit {
        const h = model.hierarchy;
        const { invariantPosition, position, x, y, z } = SymmetryOperator.createMapping(operator, model.conformation);

        return {
            id,
            model,
            operator,
            residueIndex: h.residueSegments.segmentMap,
            chainIndex: h.chainSegments.segmentMap,
            hierarchy: model.hierarchy,
            conformation: model.conformation,
            invariantPosition,
            position,
            x, y, z
        };
    }
}

export default Unit;