/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SymmetryOperator } from 'mol-math/geometry/symmetry-operator'
import { Model } from '../model'
import { GridLookup3D, Lookup3D } from 'mol-math/geometry'
import { IntraUnitLinks, computeIntraUnitBonds } from './unit/links'
import { CoarseElements, CoarseSphereConformation, CoarseGaussianConformation } from '../model/properties/coarse';
import { ValueRef } from 'mol-util';
import { UnitRings } from './unit/rings';
import StructureElement from './element'
import { ChainIndex, ResidueIndex } from '../model/indexing';
import { IntMap, SortedArray } from 'mol-data/int';
import { hash2 } from 'mol-data/util';

// A building block of a structure that corresponds to an atomic or a coarse grained representation
// 'conveniently grouped together'.
type Unit = Unit.Atomic | Unit.Spheres | Unit.Gaussians

namespace Unit {
    export const enum Kind { Atomic, Spheres, Gaussians }

    export function isAtomic(u: Unit): u is Atomic { return u.kind === Kind.Atomic; }
    export function isCoarse(u: Unit): u is Spheres | Gaussians { return u.kind === Kind.Spheres || u.kind === Kind.Gaussians; }
    export function isSpheres(u: Unit): u is Spheres { return u.kind === Kind.Spheres; }
    export function isGaussians(u: Unit): u is Gaussians { return u.kind === Kind.Gaussians; }

    export function create(id: number, invariantId: number, kind: Kind, model: Model, operator: SymmetryOperator, elements: StructureElement.Set): Unit {
        switch (kind) {
            case Kind.Atomic: return new Atomic(id, invariantId, model, elements, SymmetryOperator.createMapping(operator, model.atomicConformation, void 0), AtomicProperties());
            case Kind.Spheres: return createCoarse(id, invariantId, model, Kind.Spheres, elements, SymmetryOperator.createMapping(operator, model.coarseConformation.spheres, getSphereRadiusFunc(model)));
            case Kind.Gaussians: return createCoarse(id, invariantId, model, Kind.Gaussians, elements, SymmetryOperator.createMapping(operator, model.coarseConformation.gaussians, getGaussianRadiusFunc(model)));
        }
    }

    /** A group of units that differ only by symmetry operators. */
    export type SymmetryGroup = {
        readonly elements: StructureElement.Set
        readonly units: ReadonlyArray<Unit>
        /** Maps unit.id to index of unit in units array */
        readonly unitIndexMap: IntMap<number>
        readonly hashCode: number
    }

    function getUnitIndexMap(units: Unit[]) {
        const unitIndexMap = IntMap.Mutable<number>();
        for (let i = 0, _i = units.length; i < _i; i++) {
            unitIndexMap.set(units[i].id, i);
        }
        return unitIndexMap
    }

    export function SymmetryGroup(units: Unit[]) {
        const props: {
            unitIndexMap?: IntMap<number>
        } = {}

        return {
            elements: units[0].elements,
            units,
            get unitIndexMap () {
                if (props.unitIndexMap) return props.unitIndexMap
                props.unitIndexMap = getUnitIndexMap(units)
                return props.unitIndexMap
            },
            hashCode: hashUnit(units[0])
        }
    }

    export function conformationId (unit: Unit) {
        return Unit.isAtomic(unit) ? unit.model.atomicConformation.id : unit.model.coarseConformation.id
    }

    export function hashUnit(u: Unit) {
        return hash2(u.invariantId, SortedArray.hashCode(u.elements));
    }

    export interface Base {
        readonly id: number,
        // invariant ID stays the same even if the Operator/conformation changes.
        readonly invariantId: number,
        readonly elements: StructureElement.Set,
        readonly model: Model,
        readonly conformation: SymmetryOperator.ArrayMapping,

        getChild(elements: StructureElement.Set): Unit,
        applyOperator(id: number, operator: SymmetryOperator, dontCompose?: boolean /* = false */): Unit,

        readonly lookup3d: Lookup3D
    }

    function getSphereRadiusFunc(model: Model) {
        const r = model.coarseConformation.spheres.radius;
        return (i: number) => r[i];
    }

    function getGaussianRadiusFunc(model: Model) {
        // TODO: compute radius for gaussians
        return (i: number) => 0;
    }

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
        readonly elements: StructureElement.Set;
        readonly model: Model;
        readonly conformation: SymmetryOperator.ArrayMapping;

        // Reference some commonly accessed things for faster access.
        readonly residueIndex: ArrayLike<ResidueIndex>;
        readonly chainIndex: ArrayLike<ChainIndex>;

        private props: AtomicProperties;

        getChild(elements: StructureElement.Set): Unit {
            if (elements.length === this.elements.length) return this;
            return new Atomic(this.id, this.invariantId, this.model, elements, this.conformation, AtomicProperties());
        }

        applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
            const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
            return new Atomic(id, this.invariantId, this.model, this.elements, SymmetryOperator.createMapping(op, this.model.atomicConformation, this.conformation.r), this.props);
        }

        get lookup3d() {
            if (this.props.lookup3d.ref) return this.props.lookup3d.ref;
            const { x, y, z } = this.model.atomicConformation;
            this.props.lookup3d.ref = GridLookup3D({ x, y, z, indices: this.elements });
            return this.props.lookup3d.ref;
        }

        get links() {
            if (this.props.links.ref) return this.props.links.ref;
            this.props.links.ref = computeIntraUnitBonds(this);
            return this.props.links.ref;
        }

        get rings() {
            if (this.props.rings.ref) return this.props.rings.ref;
            this.props.rings.ref = UnitRings.create(this);
            return this.props.rings.ref;
        }

        getResidueIndex(elementIndex: StructureElement.UnitIndex) {
            return this.model.atomicHierarchy.residueAtomSegments.index[this.elements[elementIndex]];
        }

        constructor(id: number, invariantId: number, model: Model, elements: StructureElement.Set, conformation: SymmetryOperator.ArrayMapping, props: AtomicProperties) {
            this.id = id;
            this.invariantId = invariantId;
            this.model = model;
            this.elements = elements;
            this.conformation = conformation;

            this.residueIndex = model.atomicHierarchy.residueAtomSegments.index;
            this.chainIndex = model.atomicHierarchy.chainAtomSegments.index;
            this.props = props;
        }
    }

    interface AtomicProperties {
        lookup3d: ValueRef<Lookup3D | undefined>,
        links: ValueRef<IntraUnitLinks | undefined>,
        rings: ValueRef<UnitRings | undefined>
    }

    function AtomicProperties(): AtomicProperties {
        return { lookup3d: ValueRef.create(void 0), links: ValueRef.create(void 0), rings: ValueRef.create(void 0) };
    }

    class Coarse<K extends Kind.Gaussians | Kind.Spheres, C extends CoarseSphereConformation | CoarseGaussianConformation> implements Base {
        readonly kind: K;

        readonly id: number;
        readonly invariantId: number;
        readonly elements: StructureElement.Set;
        readonly model: Model;
        readonly conformation: SymmetryOperator.ArrayMapping;

        readonly coarseElements: CoarseElements;
        readonly coarseConformation: C;

        getChild(elements: StructureElement.Set): Unit {
            if (elements.length === this.elements.length) return this as any as Unit /** lets call this an ugly temporary hack */;
            return createCoarse(this.id, this.invariantId, this.model, this.kind, elements, this.conformation);
        }

        applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
            const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
            const ret = createCoarse(id, this.invariantId, this.model, this.kind, this.elements, SymmetryOperator.createMapping(op, this.getCoarseElements(), this.conformation.r));
            (ret as Coarse<K, C>)._lookup3d = this._lookup3d;
            return ret;
        }

        private _lookup3d: ValueRef<Lookup3D | undefined> = ValueRef.create(void 0);
        get lookup3d() {
            if (this._lookup3d.ref) return this._lookup3d.ref;
            // TODO: support sphere radius?
            const { x, y, z } = this.getCoarseElements();
            this._lookup3d.ref = GridLookup3D({ x, y, z, indices: this.elements });
            return this._lookup3d.ref;
        }

        private getCoarseElements() {
            return this.kind === Kind.Spheres ? this.model.coarseConformation.spheres : this.model.coarseConformation.gaussians;
        }

        constructor(id: number, invariantId: number, model: Model, kind: K, elements: StructureElement.Set, conformation: SymmetryOperator.ArrayMapping) {
            this.kind = kind;
            this.id = id;
            this.invariantId = invariantId;
            this.model = model;
            this.elements = elements;
            this.conformation = conformation;
            this.coarseElements = kind === Kind.Spheres ? model.coarseHierarchy.spheres : model.coarseHierarchy.gaussians;
            this.coarseConformation = (kind === Kind.Spheres ? model.coarseConformation.spheres : model.coarseConformation.gaussians) as C;
        }
    }

    function createCoarse<K extends Kind.Gaussians | Kind.Spheres>(id: number, invariantId: number, model: Model, kind: K, elements: StructureElement.Set, conformation: SymmetryOperator.ArrayMapping): Unit {
        return new Coarse(id, invariantId, model, kind, elements, conformation) as any as Unit /** lets call this an ugly temporary hack */;
    }

    export class Spheres extends Coarse<Kind.Spheres, CoarseSphereConformation> { }
    export class Gaussians extends Coarse<Kind.Gaussians, CoarseGaussianConformation> { }
}

export default Unit;