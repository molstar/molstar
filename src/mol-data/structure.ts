/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3, Mat4 } from '../mol-base/math/linear-algebra'
import AtomSet from './atom-set'
import Model from './model'
import Conformation from './conformation'

export interface Operator extends Readonly<{
    name: string,
    hkl: number[], // defaults to [0, 0, 0] for non symmetry entries
    transform: Mat4,
    // cache the inverse of the transform
    inverse: Mat4,
    // optimize the identity case
    isIdentity: boolean
}> { }

export interface Unit extends Readonly<{
    // Structure-level unique identifier of the unit.
    id: number,

    // Provides access to the underlying data.
    model: Model,

    // Separate the conformation from the data for faster/more straightforward dynamics
    conformation: Conformation,

    // Determines the operation applied to this unit.
    // The transform and and inverse are baked into the "getPosition" function
    operator: Operator
}> {
    // returns the untransformed position. Used for spatial queries.
    getInvariantPosition(atom: number, slot: Vec3): Vec3

    // gets the transformed position of the specified atom
    getPosition(atom: number, slot: Vec3): Vec3
}

// TODO: do "single model" version of the structure?
export interface Structure extends Readonly<{
    units: Readonly<{ [id: number]: Unit }>,
    atoms: AtomSet
}> { }

export namespace Structure {
    export const Empty: Structure = { units: {}, atoms: AtomSet.Empty };

    // TODO: "lift" atom set operators
    // TODO: "diff"
}

export default Structure