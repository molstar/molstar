/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SymmetryOperator } from 'mol-math/geometry/symmetry-operator'
import ElementGroup from './element/group'
import { Model } from '../model'
import { GridLookup3D } from 'mol-math/geometry'

// A building block of a structure that corresponds to an atomic or a coarse grained representation
// 'conveniently grouped together'.
type Unit = Unit.Atomic | Unit.Coarse

namespace Unit {
    export const enum Kind { Atomic, Coarse }

    export function isAtomic(u: Unit): u is Atomic { return u.kind === Kind.Atomic; }
    export function isCoarse(u: Unit): u is Coarse { return u.kind === Kind.Coarse; }

    export interface Base extends SymmetryOperator.ArrayMapping {
        // Provides access to the underlying data.
        readonly model: Model,

        // The "full" atom group corresponding to this unit.
        // Every selection is a subset of this atoms group.
        // Things like inter-unit bonds or spatial lookups
        // can be be implemented efficiently as "views" of the
        // full group.
        readonly fullGroup: ElementGroup,

        readonly hierarchy: Model['hierarchy'],
    }

    // A bulding block of a structure that corresponds
    // to a "natural group of atoms" (most often a "chain")
    // together with a tranformation (rotation and translation)
    // that is dynamically applied to the underlying atom set.
    //
    // An atom set can be referenced by multiple diffrent units which
    // makes construction of assemblies and spacegroups very efficient.
    export interface Atomic extends Base {
        readonly kind: Unit.Kind.Atomic,

        // Reference some commonly accessed things for faster access.
        readonly residueIndex: ArrayLike<number>,
        readonly chainIndex: ArrayLike<number>,
        readonly conformation: Model['conformation']
    }

    // Coarse grained representations.
    // TODO: can we use the ArrayMapping here?
    export interface Coarse extends Base  {
        readonly kind: Unit.Kind.Coarse
    }

    export function createAtomic(model: Model, operator: SymmetryOperator, fullGroup: ElementGroup): Unit {
        const h = model.hierarchy;
        const { invariantPosition, position, x, y, z } = SymmetryOperator.createMapping(operator, model.conformation);

        return {
            model,
            kind: Kind.Atomic,
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

    export function createCoarse(model: Model, operator: SymmetryOperator, fullGroup: ElementGroup): Unit {
        throw 'not implemented'
    }

    export function withOperator(unit: Unit, operator: SymmetryOperator): Unit {
        switch (unit.kind) {
            case Kind.Atomic: return createAtomic(unit.model, SymmetryOperator.compose(unit.operator, operator), unit.fullGroup);
            case Kind.Coarse: return createCoarse(unit.model, SymmetryOperator.compose(unit.operator, operator), unit.fullGroup);
        }
    }

    export function getLookup3d(unit: Unit, group: ElementGroup) {
        if (group.__lookup3d__)  return group.__lookup3d__;
        if (Unit.isAtomic(unit)) {
            const { x, y, z } = unit.model.conformation;
            group.__lookup3d__ = GridLookup3D({ x, y, z, indices: group.elements });
            return group.__lookup3d__;
        }

        throw 'not implemented';
    }
}

export default Unit;