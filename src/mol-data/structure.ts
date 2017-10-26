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
    hkl: number[], // defaults to [0, 0, 0] where not appropriate
    transform: Mat4,
    inverse: Mat4,
    isIdentity: boolean
}> { }

export interface Unit extends Readonly<{
    // Structure-level unique identifier of the unit.
    id: number,

    // Each unit can only contain atoms from a single "chain"
    // the reason for this is to make certain basic data transforms fast
    // without having to look at the actual contents of the unit
    // multiple units can point to the same chain
    chainIndex: number,

    // Provides access to the underlying data.
    model: Model,

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

export interface Structure { units: { [id: number]: Unit }, atoms: AtomSet }

export namespace Structure {
    export const Empty: Structure = { units: {}, atoms: AtomSet.Empty };
}

export default Structure
// export interface Selection { structure: Structure, sets: AtomSet[] }
// type SelectionImpl = Structure | Structure[]