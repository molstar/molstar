/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SymmetryOperator } from 'mol-math/geometry/symmetry-operator'
import { Model } from '../model'
import { GridLookup3D, Lookup3D } from 'mol-math/geometry'
import { CoarseGrainedHierarchy } from '../model/properties/coarse-grained/hierarchy';
import { SortedArray } from 'mol-data/int';
import { idFactory } from 'mol-util/id-factory';
import { IntraUnitBonds, computeIntraUnitBonds } from './unit/bonds'

// A building block of a structure that corresponds to an atomic or a coarse grained representation
// 'conveniently grouped together'.
type Unit = Unit.Atomic | Unit.Spheres | Unit.Gaussians

namespace Unit {
    export const enum Kind { Atomic, Spheres, Gaussians }

    export function isAtomic(u: Unit): u is Atomic { return u.kind === Kind.Atomic; }
    export function isCoarse(u: Unit): u is Spheres | Gaussians { return u.kind === Kind.Spheres || u.kind === Kind.Gaussians; }
    export function isSpheres(u: Unit): u is Spheres { return u.kind === Kind.Spheres; }
    export function isGaussians(u: Unit): u is Gaussians { return u.kind === Kind.Gaussians; }

    export function create(id: number, kind: Kind, model: Model, operator: SymmetryOperator, elements: SortedArray): Unit {
        switch (kind) {
            case Kind.Atomic: return new Atomic(id, unitIdFactory(), model, elements, SymmetryOperator.createMapping(operator, model.atomSiteConformation));
            case Kind.Spheres: return createCoarse(id, unitIdFactory(), model, Kind.Spheres, model.coarseGrained.spheres, elements, SymmetryOperator.createMapping(operator, model.atomSiteConformation));
            case Kind.Gaussians: return createCoarse(id, unitIdFactory(), model, Kind.Gaussians, model.coarseGrained.gaussians, elements, SymmetryOperator.createMapping(operator, model.atomSiteConformation));
        }
    }

    // A group of units that differ only by symmetry operators.
    export type SymmetryGroup = { readonly elements: SortedArray, readonly units: ReadonlyArray<Unit> }

    export interface Base {
        readonly id: number,
        // invariant ID stays the same even if the Operator/conformation changes.
        readonly invariantId: number,
        readonly elements: SortedArray,
        readonly model: Model,
        readonly conformation: SymmetryOperator.ArrayMapping,

        getChild(elements: SortedArray): Unit,
        applyOperator(id: number, operator: SymmetryOperator, dontCompose?: boolean /* = false */): Unit,

        readonly lookup3d: Lookup3D
    }

    const unitIdFactory = idFactory();

    // A bulding block of a structure that corresponds
    // to a "natural group of atoms" (most often a "chain")
    // together with a tranformation (rotation and translation)
    // that is dynamically applied to the underlying atom set.
    //
    // An atom set can be referenced by multiple diffrent units which
    // makes construction of assemblies and spacegroups very efficient.
    export class Atomic implements Base {
        readonly kind = Kind.Atomic;

        readonly id: number;
        readonly invariantId: number;
        readonly elements: SortedArray;
        readonly model: Model;
        readonly conformation: SymmetryOperator.ArrayMapping;

        // Reference some commonly accessed things for faster access.
        readonly residueIndex: ArrayLike<number>;
        readonly chainIndex: ArrayLike<number>;

        getChild(elements: SortedArray): Unit {
            if (elements.length === this.elements.length) return this;
            return new Atomic(this.id, this.invariantId, this.model, elements, this.conformation);
        }

        applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
            const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
            return new Atomic(id, this.invariantId, this.model, this.elements, SymmetryOperator.createMapping(op, this.model.atomSiteConformation));
        }

        private _lookup3d?: Lookup3D = void 0;
        get lookup3d() {
            if (this._lookup3d) return this._lookup3d;
            const { x, y, z } = this.model.atomSiteConformation;
            this._lookup3d = GridLookup3D({ x, y, z, indices: this.elements });
            return this._lookup3d;
        }

        private _bonds?: IntraUnitBonds = void 0;
        get bonds() {
            if (this._bonds) return this._bonds;
            this._bonds = computeIntraUnitBonds(this);
            return this._bonds;
        }

        constructor(id: number, invariantId: number, model: Model, elements: SortedArray, conformation: SymmetryOperator.ArrayMapping) {
            this.id = id;
            this.invariantId = invariantId;
            this.model = model;
            this.elements = elements;
            this.conformation = conformation;

            this.residueIndex = model.hierarchy.residueSegments.segmentMap;
            this.chainIndex = model.hierarchy.chainSegments.segmentMap;
        }
    }

    // Coarse grained representations.
    export interface CoarseBase<S extends CoarseGrainedHierarchy.SitesBase> extends Base  {
        readonly sites: S
    }

    class Coarse<S extends CoarseGrainedHierarchy.SitesBase> implements CoarseBase<S> {
        readonly kind: Kind;

        readonly id: number;
        readonly invariantId: number;
        readonly elements: SortedArray;
        readonly model: Model;
        readonly conformation: SymmetryOperator.ArrayMapping;
        readonly sites: S;

        getChild(elements: SortedArray): Unit {
            if (elements.length === this.elements.length) return this as any as Unit /** lets call this an ugly temporary hack */;
            return createCoarse(this.id, this.invariantId, this.model, this.kind, this.sites, elements, this.conformation);
        }

        applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
            const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
            return createCoarse(id, this.invariantId, this.model, this.kind, this.sites, this.elements, SymmetryOperator.createMapping(op, this.sites));
        }

        private _lookup3d?: Lookup3D = void 0;
        get lookup3d() {
            if (this._lookup3d) return this._lookup3d;
            const { x, y, z } = this.sites;
            // TODO: support sphere radius
            this._lookup3d = GridLookup3D({ x, y, z, indices: this.elements });
            return this._lookup3d;
        }

        constructor(id: number, invariantId: number, model: Model, kind: Kind, sites: S, elements: SortedArray, conformation: SymmetryOperator.ArrayMapping) {
            this.kind = kind;
            this.id = id;
            this.invariantId = invariantId;
            this.model = model;
            this.elements = elements;
            this.conformation = conformation;
            this.sites = sites;
        }
    }

    function createCoarse<S extends CoarseGrainedHierarchy.SitesBase>(id: number, invariantId: number, model: Model, kind: Kind, sites: S, elements: SortedArray, conformation: SymmetryOperator.ArrayMapping): Unit {
        return new Coarse(id, invariantId, model, kind, sites, elements, conformation) as any as Unit /** lets call this an ugly temporary hack */;
    }

    export interface Spheres extends CoarseBase<CoarseGrainedHierarchy.Spheres> { kind: Kind.Spheres }
    export interface Gaussians extends CoarseBase<CoarseGrainedHierarchy.Gaussians> { kind: Kind.Gaussians }
}

export default Unit;