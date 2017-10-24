/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3, Mat4 } from '../mol-base/math/linear-algebra'
import AtomSet from './atom-set'
import Model from './model'

export type Operator =
    | { kind: Operator.Kind.Identity }
    | { kind: Operator.Kind.Symmetry, hkl: number[], index: number, name: string, transform: Mat4, inverse: Mat4 }
    | { kind: Operator.Kind.Assembly, assemblyName: string, index: number, transform: Mat4, inverse: Mat4 }
    | { kind: Operator.Kind.Custom, name: string, transform: Mat4, inverse: Mat4 }

export namespace Operator {
    export enum Kind { Identity, Symmetry, Assembly, Custom }
}

export interface Unit extends Readonly<{
    // Structure-level unique identifier of the unit.
    id: number,

    // Each unit can only contain atoms from a single "chain"
    // the reason for this is to make symmetry and assembly transforms fast
    // without having to look at the actual contents of the unit
    // multiple units can point to the same chain
    chainIndex: number,

    // Provides access to the underlying data.
    model: Model,

    // Determines the operation applied to this unit.
    // The transform and and inverse a baked into the "getPosition" function
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

    export enum Algebra { AddUnit, RemoveUnit, UpdateConformation /* specify which units map to which */ }
}

// export interface Selection { structure: Structure, sets: AtomSet[] }
// type SelectionImpl = Structure | Structure[]