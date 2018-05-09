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
import { SortedArray } from 'mol-data/int';
import { idFactory } from 'mol-util/id-factory';

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

    export function applyOperator(id: number, unit: Unit, operator: SymmetryOperator): Unit {
        return create(id, unit.kind, unit.model, SymmetryOperator.compose(unit.conformation.operator, operator), unit.elements);
    }

    export interface Base {
        readonly id: number,
        // invariant ID stays the same even if the Operator/conformation changes.
        readonly invariantId: number,
        readonly elements: SortedArray,
        readonly model: Model,
        readonly conformation: SymmetryOperator.ArrayMapping,

        getChild(elements: SortedArray): Unit,
        applyOperator(id: number, operator: SymmetryOperator): Unit
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
            return new Atomic(this.id, this.invariantId, this.model, elements, this.conformation);
        }

        applyOperator(id: number, operator: SymmetryOperator): Unit {
            const op = SymmetryOperator.compose(this.conformation.operator, operator);
            return new Atomic(id, this.invariantId, this.model, this.elements, SymmetryOperator.createMapping(op, this.model.atomSiteConformation));
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
    export interface CoarseBase<S extends CoarseGrained.SitesBase> extends Base  {
        readonly sites: S
    }

    class Coarse<S extends CoarseGrained.SitesBase> implements CoarseBase<S> {
        readonly kind: Kind;

        readonly id: number;
        readonly invariantId: number;
        readonly elements: SortedArray;
        readonly model: Model;
        readonly conformation: SymmetryOperator.ArrayMapping;
        readonly sites: S;

        getChild(elements: SortedArray): Unit {
            return createCoarse(this.id, this.invariantId, this.model, this.kind, this.sites, elements, this.conformation);
        }

        applyOperator(id: number, operator: SymmetryOperator): Unit {
            const op = SymmetryOperator.compose(this.conformation.operator, operator);
            return createCoarse(id, this.invariantId, this.model, this.kind, this.sites, this.elements, SymmetryOperator.createMapping(op, this.sites));
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

    function createCoarse<S extends CoarseGrained.SitesBase>(id: number, invariantId: number, model: Model, kind: Kind, sites: S, elements: SortedArray, conformation: SymmetryOperator.ArrayMapping): Unit {
        return new Coarse(id, invariantId, model, kind, sites, elements, conformation) as any as Unit /** lets call this an ugly temporary hack */;
    }

    export interface Spheres extends CoarseBase<CoarseGrained.Spheres> { kind: Kind.Spheres }
    export interface Gaussians extends CoarseBase<CoarseGrained.Gaussians> { kind: Kind.Gaussians }

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