/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Model from '../model'
import Operator from './operator'

interface Unit extends Readonly<{
    // Structure-level unique identifier of the unit.
    id: number,

    // Provides access to the underlying data.
    model: Model,

    // Determines the operation applied to this unit.
    // The transform and and inverse are baked into the "getPosition" function
    operator: Operator,

    // Cache residue and chain indices for fast access.
    residueIndex: ArrayLike<number>,
    chainIndex: ArrayLike<number>,
    hierarchy: Model['hierarchy'],
    conformation: Model['conformation']
}> {
    // // returns the untransformed position. Used for spatial queries.
    // getInvariantPosition(atom: number, slot: Vec3): Vec3

    // // gets the transformed position of the specified atom
    // getPosition(atom: number, slot: Vec3): Vec3
}

namespace Unit {
    export function create(model: Model, operator: Operator): Unit {
        const h = model.hierarchy;
        return {
            id: nextUnitId(),
            model,
            operator,
            residueIndex: h.residueSegments.segmentMap,
            chainIndex: h.chainSegments.segmentMap,
            hierarchy: model.hierarchy,
            conformation: model.conformation
        };
    }
}

export default Unit;

let _id = 0;
function nextUnitId() {
    const ret = _id;
    _id = (_id + 1) % 0x3fffffff;
    return ret;
}