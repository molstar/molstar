/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import SymmetryOperator from 'mol-math/geometry/symmetry-operator'
import AtomGroup from './atom/group'
import { Model } from '../model'

// A bulding block of a structure that corresponds
// to a "natural group of atoms" (most often a "chain")
// together with a tranformation (rotation and translation)
// that is dynamically applied to the underlying atom set.
//
// An atom set can be referenced by multiple diffrent units which
// makes construction of assemblies and spacegroups very efficient.
interface Unit extends SymmetryOperator.ArrayMapping {
    // Provides access to the underlying data.
    readonly model: Model,

    // The "full" atom group corresponding to this unit.
    // Every selection is a subset of this atoms group.
    // Things like inter-unit bonds or spatial lookups
    // can be be implemented efficiently as "views" of the
    // full group.
    readonly fullGroup: AtomGroup,

    // Reference some commonly accessed things for faster access.
    readonly residueIndex: ArrayLike<number>,
    readonly chainIndex: ArrayLike<number>,
    readonly hierarchy: Model['hierarchy'],
    readonly conformation: Model['conformation']

    // TODO: add velocity?
}

namespace Unit {
    export function create(model: Model, operator: SymmetryOperator, fullGroup: AtomGroup): Unit {
        const h = model.hierarchy;
        const { invariantPosition, position, x, y, z } = SymmetryOperator.createMapping(operator, model.conformation);

        return {
            model,
            operator,
            fullGroup,
            residueIndex: h.residueSegments.segmentMap,
            chainIndex: h.chainSegments.segmentMap,
            hierarchy: model.hierarchy,
            conformation: model.conformation,
            invariantPosition,
            position,
            x, y, z
        };
    }

    export function withOperator(unit: Unit, operator: SymmetryOperator) {
        return create(unit.model, SymmetryOperator.compose(unit.operator, operator), unit.fullGroup);
    }
}

export default Unit;