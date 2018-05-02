/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SymmetryOperator } from 'mol-math/geometry/symmetry-operator'
import ElementGroup from './element/group'
import { Model } from '../model'
import { GridLookup3D } from 'mol-math/geometry'
import { computeUnitBonds } from './element/properties/bonds/group-compute';
import CoarseGrained from '../model/properties/coarse-grained';

// A building block of a structure that corresponds to an atomic or a coarse grained representation
// 'conveniently grouped together'.
type Unit = Unit.Atomic | Unit.CoarseSpheres | Unit.CoarseGaussians

namespace Unit {
    export const enum Kind { Atomic, CoarseSpheres, CoarseGaussians }

    export function isAtomic(u: Unit): u is Atomic { return u.kind === Kind.Atomic; }
    export function isCoarse(u: Unit): u is CoarseSpheres | CoarseGaussians { return u.kind === Kind.CoarseSpheres || u.kind === Kind.CoarseGaussians; }
    export function isCoarseSpheres(u: Unit): u is CoarseSpheres { return u.kind === Kind.CoarseSpheres; }
    export function isCoarseGaussians(u: Unit): u is CoarseGaussians { return u.kind === Kind.CoarseGaussians; }

    export interface Base extends SymmetryOperator.ArrayMapping {
        // Provides access to the underlying data.
        readonly model: Model,

        // The "full" atom group corresponding to this unit.
        // Every selection is a subset of this atoms group.
        // Things like inter-unit bonds or spatial lookups
        // can be be implemented efficiently as "views" of the
        // full group.
        readonly fullGroup: ElementGroup
    }

    // A bulding block of a structure that corresponds
    // to a "natural group of atoms" (most often a "chain")
    // together with a tranformation (rotation and translation)
    // that is dynamically applied to the underlying atom set.
    //
    // An atom set can be referenced by multiple diffrent units which
    // makes construction of assemblies and spacegroups very efficient.
    export interface Atomic extends Base {
        readonly kind: Kind.Atomic,

        // Reference some commonly accessed things for faster access.
        readonly residueIndex: ArrayLike<number>,
        readonly chainIndex: ArrayLike<number>,
        readonly conformation: Model['atomSiteConformation'],
        readonly hierarchy: Model['hierarchy']
    }

    // Coarse grained representations.
    export interface CoarseBase<K extends Kind, S extends CoarseGrained.SitesBase> extends Base  {
        readonly kind: K,
        readonly sites: S
    }

    export interface CoarseSpheres extends CoarseBase<Kind.CoarseSpheres, CoarseGrained.Spheres> { }
    export interface CoarseGaussians extends CoarseBase<Kind.CoarseGaussians, CoarseGrained.Gaussians> { }

    export function create(kind: Kind, model: Model, operator: SymmetryOperator, fullGroup: ElementGroup): Unit {
        switch (kind) {
            case Kind.Atomic: return createAtomic(model, operator, fullGroup);
            case Kind.CoarseSpheres: return createCoarseSpheres(model, operator, fullGroup);
            case Kind.CoarseGaussians: return createCoarseGaussians(model, operator, fullGroup);
        }
    }

    function createAtomic(model: Model, operator: SymmetryOperator, fullGroup: ElementGroup): Unit.Atomic {
        const h = model.hierarchy;
        const { invariantPosition, position, x, y, z } = SymmetryOperator.createMapping(operator, model.atomSiteConformation);

        return {
            model,
            kind: Kind.Atomic,
            operator,
            fullGroup,
            residueIndex: h.residueSegments.segmentMap,
            chainIndex: h.chainSegments.segmentMap,
            hierarchy: model.hierarchy,
            conformation: model.atomSiteConformation,
            invariantPosition,
            position,
            x, y, z
        };
    }

    function createCoarseSpheres(model: Model, operator: SymmetryOperator, fullGroup: ElementGroup): Unit.CoarseSpheres {
        const { invariantPosition, position, x, y, z } = SymmetryOperator.createMapping(operator, model.coarseGrained.spheres);

        return {
            model,
            kind: Kind.CoarseSpheres,
            sites: model.coarseGrained.spheres,
            operator,
            fullGroup,
            invariantPosition,
            position,
            x, y, z
        };
    }

    function createCoarseGaussians(model: Model, operator: SymmetryOperator, fullGroup: ElementGroup): Unit.CoarseGaussians {
        const { invariantPosition, position, x, y, z } = SymmetryOperator.createMapping(operator, model.coarseGrained.gaussians);

        return {
            model,
            kind: Kind.CoarseGaussians,
            sites: model.coarseGrained.gaussians,
            operator,
            fullGroup,
            invariantPosition,
            position,
            x, y, z
        };
    }

    export function withOperator(unit: Unit, operator: SymmetryOperator): Unit {
        return create(unit.kind, unit.model, SymmetryOperator.compose(unit.operator, operator), unit.fullGroup);
    }

    export function getLookup3d(unit: Unit, group: ElementGroup) {
        if (group.__lookup3d__)  return group.__lookup3d__;
        if (Unit.isAtomic(unit)) {
            const { x, y, z } = unit.model.atomSiteConformation;
            group.__lookup3d__ = GridLookup3D({ x, y, z, indices: group.elements });
            return group.__lookup3d__;
        }

        throw 'not implemented';
    }

    export function getGroupBonds(unit: Unit, group: ElementGroup) {
        if (group.__bonds__) return group.__bonds__;
        if (Unit.isAtomic(unit)) {
            group.__bonds__ = computeUnitBonds(unit, group);
            return group.__bonds__;
        }

        throw 'not implemented';
    }
}

export default Unit;