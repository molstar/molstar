/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SymmetryOperator } from '../../../mol-math/geometry/symmetry-operator';
import { Model } from '../model';
import { GridLookup3D, Lookup3D } from '../../../mol-math/geometry';
import { IntraUnitBonds, computeIntraUnitBonds } from './unit/bonds';
import { CoarseElements, CoarseSphereConformation, CoarseGaussianConformation } from '../model/properties/coarse';
import { ValueRef, BitFlags } from '../../../mol-util';
import { UnitRings } from './unit/rings';
import StructureElement from './element';
import { ChainIndex, ResidueIndex, ElementIndex } from '../model/indexing';
import { IntMap, SortedArray, Segmentation } from '../../../mol-data/int';
import { hash2, hashFnv32a } from '../../../mol-data/util';
import { getAtomicPolymerElements, getCoarsePolymerElements, getAtomicGapElements, getCoarseGapElements, getNucleotideElements, getProteinElements } from './util/polymer';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { getPrincipalAxes } from './util/principal-axes';
import { Boundary, getBoundary } from '../../../mol-math/geometry/boundary';

/**
 * A building block of a structure that corresponds to an atomic or
 * a coarse grained representation 'conveniently grouped together'.
 */
type Unit = Unit.Atomic | Unit.Spheres | Unit.Gaussians

namespace Unit {
    export const enum Kind { Atomic, Spheres, Gaussians }

    export function isAtomic(u: Unit): u is Atomic { return u.kind === Kind.Atomic; }
    export function isCoarse(u: Unit): u is Spheres | Gaussians { return u.kind === Kind.Spheres || u.kind === Kind.Gaussians; }
    export function isSpheres(u: Unit): u is Spheres { return u.kind === Kind.Spheres; }
    export function isGaussians(u: Unit): u is Gaussians { return u.kind === Kind.Gaussians; }

    export function create(id: number, invariantId: number, chainGroupId: number, traits: Traits, kind: Kind, model: Model, operator: SymmetryOperator, elements: StructureElement.Set): Unit {
        switch (kind) {
            case Kind.Atomic: return new Atomic(id, invariantId, chainGroupId, traits, model, elements, SymmetryOperator.createMapping(operator, model.atomicConformation, void 0), AtomicProperties());
            case Kind.Spheres: return createCoarse(id, invariantId, chainGroupId, traits, model, Kind.Spheres, elements, SymmetryOperator.createMapping(operator, model.coarseConformation.spheres, getSphereRadiusFunc(model)), CoarseProperties());
            case Kind.Gaussians: return createCoarse(id, invariantId, chainGroupId, traits, model, Kind.Gaussians, elements, SymmetryOperator.createMapping(operator, model.coarseConformation.gaussians, getGaussianRadiusFunc(model)), CoarseProperties());
        }
    }

    /** A group of units that differ only by symmetry operators. */
    export type SymmetryGroup = {
        readonly elements: StructureElement.Set
        readonly units: ReadonlyArray<Unit>
        /** Maps unit.id to index of unit in units array */
        readonly unitIndexMap: IntMap<number>
        /** Hash based on unit.invariantId which is the same for all units in the group */
        readonly hashCode: number
        /** Hash based on all unit.id values in the group, reflecting the units transformation*/
        readonly transformHash: number
    }

    function getUnitIndexMap(units: Unit[]) {
        const unitIndexMap = IntMap.Mutable<number>();
        for (let i = 0, _i = units.length; i < _i; i++) {
            unitIndexMap.set(units[i].id, i);
        }
        return unitIndexMap;
    }

    export function SymmetryGroup(units: Unit[]) {
        const props: {
            unitIndexMap?: IntMap<number>
        } = {};

        return {
            elements: units[0].elements,
            units,
            get unitIndexMap () {
                if (props.unitIndexMap) return props.unitIndexMap;
                props.unitIndexMap = getUnitIndexMap(units);
                return props.unitIndexMap;
            },
            hashCode: hashUnit(units[0]),
            transformHash: hashFnv32a(units.map(u => u.id))
        };
    }

    export namespace SymmetryGroup {
        export function areInvariantElementsEqual(a: SymmetryGroup, b: SymmetryGroup) {
            if (a.hashCode !== b.hashCode) return false;
            return SortedArray.areEqual(a.elements, b.elements);
        }

        export function getUnitSymmetryGroupsIndexMap(symmetryGroups: ReadonlyArray<Unit.SymmetryGroup>): IntMap<number> {
            const unitSymmetryGroupsIndexMap = IntMap.Mutable<number>();
            for (let i = 0, _i = symmetryGroups.length; i < _i; i++) {
                unitSymmetryGroupsIndexMap.set(symmetryGroups[i].units[0].invariantId, i);
            }
            return unitSymmetryGroupsIndexMap;
        }
    }

    export function conformationId (unit: Unit) {
        return Unit.isAtomic(unit) ? unit.model.atomicConformation.id : unit.model.coarseConformation.id;
    }

    export function hashUnit(u: Unit) {
        return hash2(u.invariantId, SortedArray.hashCode(u.elements));
    }

    export type Traits = BitFlags<Trait>
    export const enum Trait {
        None = 0x0,
        MultiChain = 0x1,
        Partitioned = 0x2
    }
    export namespace Traits {
        export const is: (t: Traits, f: Trait) => boolean = BitFlags.has;
        export const create: (f: Trait) => Traits = BitFlags.create;
    }

    export interface Base {
        readonly id: number,
        /** invariant ID stays the same even if the Operator/conformation changes. */
        readonly invariantId: number,
        readonly chainGroupId: number,
        readonly traits: Traits,
        readonly elements: StructureElement.Set,
        readonly model: Model,
        readonly conformation: SymmetryOperator.ArrayMapping<ElementIndex>,

        getChild(elements: StructureElement.Set): Unit,
        applyOperator(id: number, operator: SymmetryOperator, dontCompose?: boolean /* = false */): Unit,

        readonly boundary: Boundary
        readonly lookup3d: Lookup3D<StructureElement.UnitIndex>
        readonly polymerElements: SortedArray<ElementIndex>
        readonly gapElements: SortedArray<ElementIndex>
        /**
         * From mmCIF/IHM schema: `_ihm_model_representation_details.model_object_primitive`.
         */
        readonly objectPrimitive: mmCIF_Schema['ihm_model_representation_details']['model_object_primitive']['T']
    }

    interface BaseProperties {
        boundary: ValueRef<Boundary | undefined>,
        lookup3d: ValueRef<Lookup3D<StructureElement.UnitIndex> | undefined>,
        principalAxes: ValueRef<PrincipalAxes | undefined>,
        polymerElements: ValueRef<SortedArray<ElementIndex> | undefined>
        gapElements: ValueRef<SortedArray<ElementIndex> | undefined>
    }


    function BaseProperties(): BaseProperties {
        return {
            boundary: ValueRef.create(void 0),
            lookup3d: ValueRef.create(void 0),
            principalAxes: ValueRef.create(void 0),
            polymerElements: ValueRef.create(void 0),
            gapElements: ValueRef.create(void 0),
        };
    }

    function getSphereRadiusFunc(model: Model) {
        const r = model.coarseConformation.spheres.radius;
        return (i: number) => r[i];
    }

    function getGaussianRadiusFunc(model: Model) {
        // TODO: compute radius for gaussians
        return (i: number) => 0;
    }

    /**
     * A bulding block of a structure that corresponds
     * to a "natural group of atoms" (most often a "chain")
     * together with a transformation (rotation and translation)
     * that is dynamically applied to the underlying atom set.
     *
     * An atom set can be referenced by multiple different units which
     * makes construction of assemblies and spacegroups very efficient.
     */
    export class Atomic implements Base {
        readonly kind = Kind.Atomic;
        readonly objectPrimitive = 'atomistic';

        readonly id: number;
        readonly invariantId: number;
        /** Used to identify a single chain split into multiple units. */
        readonly chainGroupId: number;
        readonly traits: Traits;
        readonly elements: StructureElement.Set;
        readonly model: Model;
        readonly conformation: SymmetryOperator.ArrayMapping<ElementIndex>;

        /** Reference `residueIndex` from `model` for faster access. */
        readonly residueIndex: ArrayLike<ResidueIndex>;
        /** Reference `chainIndex` from `model` for faster access. */
        readonly chainIndex: ArrayLike<ChainIndex>;

        private props: AtomicProperties;

        getChild(elements: StructureElement.Set): Unit {
            if (elements.length === this.elements.length) return this;
            return new Atomic(this.id, this.invariantId, this.chainGroupId, this.traits, this.model, elements, this.conformation, AtomicProperties());
        }

        applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
            const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
            return new Atomic(id, this.invariantId, this.chainGroupId, this.traits, this.model, this.elements, SymmetryOperator.createMapping(op, this.model.atomicConformation, this.conformation.r), this.props);
        }

        get boundary() {
            if (this.props.boundary.ref) return this.props.boundary.ref;
            const { x, y, z } = this.model.atomicConformation;
            this.props.boundary.ref = getBoundary({ x, y, z, indices: this.elements });
            return this.props.boundary.ref;
        }

        get lookup3d() {
            if (this.props.lookup3d.ref) return this.props.lookup3d.ref;
            const { x, y, z } = this.model.atomicConformation;
            this.props.lookup3d.ref = GridLookup3D({ x, y, z, indices: this.elements }, this.boundary);
            return this.props.lookup3d.ref;
        }

        get principalAxes() {
            if (this.props.principalAxes.ref) return this.props.principalAxes.ref;
            this.props.principalAxes.ref = getPrincipalAxes(this);
            return this.props.principalAxes.ref;
        }

        get bonds() {
            if (this.props.bonds.ref) return this.props.bonds.ref;
            this.props.bonds.ref = computeIntraUnitBonds(this);
            return this.props.bonds.ref;
        }

        get rings() {
            if (this.props.rings.ref) return this.props.rings.ref;
            this.props.rings.ref = UnitRings.create(this);
            return this.props.rings.ref;
        }

        get polymerElements() {
            if (this.props.polymerElements.ref) return this.props.polymerElements.ref;
            this.props.polymerElements.ref = getAtomicPolymerElements(this);
            return this.props.polymerElements.ref;
        }

        get gapElements() {
            if (this.props.gapElements.ref) return this.props.gapElements.ref;
            this.props.gapElements.ref = getAtomicGapElements(this);
            return this.props.gapElements.ref;
        }

        get nucleotideElements() {
            if (this.props.nucleotideElements.ref) return this.props.nucleotideElements.ref;
            this.props.nucleotideElements.ref = getNucleotideElements(this);
            return this.props.nucleotideElements.ref;
        }

        get proteinElements() {
            if (this.props.proteinElements.ref) return this.props.proteinElements.ref;
            this.props.proteinElements.ref = getProteinElements(this);
            return this.props.proteinElements.ref;
        }

        get residueCount(): number {
            if (this.props.residueCount.ref !== undefined) return this.props.residueCount.ref;

            let residueCount = 0;
            const residueIt = Segmentation.transientSegments(this.model.atomicHierarchy.residueAtomSegments, this.elements);
            while (residueIt.hasNext) {
                residueIt.move();
                residueCount += 1;
            }

            this.props.residueCount.ref = residueCount;
            return this.props.residueCount.ref!;
        }

        getResidueIndex(elementIndex: StructureElement.UnitIndex) {
            return this.residueIndex[this.elements[elementIndex]];
        }

        constructor(id: number, invariantId: number, chainGroupId: number, traits: Traits, model: Model, elements: StructureElement.Set, conformation: SymmetryOperator.ArrayMapping<ElementIndex>, props: AtomicProperties) {
            this.id = id;
            this.invariantId = invariantId;
            this.chainGroupId = chainGroupId;
            this.traits = traits;
            this.model = model;
            this.elements = elements;
            this.conformation = conformation;

            this.residueIndex = model.atomicHierarchy.residueAtomSegments.index;
            this.chainIndex = model.atomicHierarchy.chainAtomSegments.index;
            this.props = props;
        }
    }

    interface AtomicProperties extends BaseProperties {
        bonds: ValueRef<IntraUnitBonds | undefined>,
        rings: ValueRef<UnitRings | undefined>
        nucleotideElements: ValueRef<SortedArray<ElementIndex> | undefined>
        proteinElements: ValueRef<SortedArray<ElementIndex> | undefined>
        residueCount: ValueRef<number | undefined>
    }

    function AtomicProperties(): AtomicProperties {
        return {
            ...BaseProperties(),
            bonds: ValueRef.create(void 0),
            rings: ValueRef.create(void 0),
            nucleotideElements: ValueRef.create(void 0),
            proteinElements: ValueRef.create(void 0),
            residueCount: ValueRef.create(void 0),
        };
    }

    class Coarse<K extends Kind.Gaussians | Kind.Spheres, C extends CoarseSphereConformation | CoarseGaussianConformation> implements Base {
        readonly kind: K;
        readonly objectPrimitive: 'sphere' | 'gaussian';

        readonly id: number;
        readonly invariantId: number;
        readonly chainGroupId: number;
        readonly traits: Traits;
        readonly elements: StructureElement.Set;
        readonly model: Model;
        readonly conformation: SymmetryOperator.ArrayMapping<ElementIndex>;

        readonly coarseElements: CoarseElements;
        readonly coarseConformation: C;

        private props: CoarseProperties;

        getChild(elements: StructureElement.Set): Unit {
            if (elements.length === this.elements.length) return this as any as Unit /** lets call this an ugly temporary hack */;
            return createCoarse(this.id, this.invariantId, this.chainGroupId, this.traits, this.model, this.kind, elements, this.conformation, CoarseProperties());
        }

        applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
            const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
            const ret = createCoarse(id, this.invariantId, this.chainGroupId, this.traits, this.model, this.kind, this.elements, SymmetryOperator.createMapping(op, this.getCoarseConformation(), this.conformation.r), this.props);
            // (ret as Coarse<K, C>)._lookup3d = this._lookup3d;
            return ret;
        }

        get boundary() {
            if (this.props.boundary.ref) return this.props.boundary.ref;
            // TODO: support sphere radius?
            const { x, y, z } = this.getCoarseConformation();
            this.props.boundary.ref = getBoundary({ x, y, z, indices: this.elements });
            return this.props.boundary.ref;
        }

        get lookup3d() {
            if (this.props.lookup3d.ref) return this.props.lookup3d.ref;
            // TODO: support sphere radius?
            const { x, y, z } = this.getCoarseConformation();
            this.props.lookup3d.ref = GridLookup3D({ x, y, z, indices: this.elements }, this.boundary);
            return this.props.lookup3d.ref;
        }

        get principalAxes() {
            if (this.props.principalAxes.ref) return this.props.principalAxes.ref;
            this.props.principalAxes.ref = getPrincipalAxes(this as Unit.Spheres | Unit.Gaussians); // TODO get rid of casting
            return this.props.principalAxes.ref;
        }

        get polymerElements() {
            if (this.props.polymerElements.ref) return this.props.polymerElements.ref;
            this.props.polymerElements.ref = getCoarsePolymerElements(this as Unit.Spheres | Unit.Gaussians); // TODO get rid of casting
            return this.props.polymerElements.ref;
        }

        get gapElements() {
            if (this.props.gapElements.ref) return this.props.gapElements.ref;
            this.props.gapElements.ref = getCoarseGapElements(this as Unit.Spheres | Unit.Gaussians); // TODO get rid of casting
            return this.props.gapElements.ref;
        }

        private getCoarseConformation() {
            return this.kind === Kind.Spheres ? this.model.coarseConformation.spheres : this.model.coarseConformation.gaussians;
        }

        constructor(id: number, invariantId: number, chainGroupId: number, traits: Traits, model: Model, kind: K, elements: StructureElement.Set, conformation: SymmetryOperator.ArrayMapping<ElementIndex>, props: CoarseProperties) {
            this.kind = kind;
            this.objectPrimitive = kind === Kind.Spheres ? 'sphere' : 'gaussian';
            this.id = id;
            this.invariantId = invariantId;
            this.chainGroupId = chainGroupId;
            this.traits = traits;
            this.model = model;
            this.elements = elements;
            this.conformation = conformation;
            this.coarseElements = kind === Kind.Spheres ? model.coarseHierarchy.spheres : model.coarseHierarchy.gaussians;
            this.coarseConformation = (kind === Kind.Spheres ? model.coarseConformation.spheres : model.coarseConformation.gaussians) as C;
            this.props = props;
        }
    }

    interface CoarseProperties extends BaseProperties { }

    function CoarseProperties(): CoarseProperties {
        return BaseProperties();
    }

    function createCoarse<K extends Kind.Gaussians | Kind.Spheres>(id: number, invariantId: number, chainGroupId: number, traits: Traits, model: Model, kind: K, elements: StructureElement.Set, conformation: SymmetryOperator.ArrayMapping<ElementIndex>, props: CoarseProperties): Unit {
        return new Coarse(id, invariantId, chainGroupId, traits, model, kind, elements, conformation, props) as any as Unit /** lets call this an ugly temporary hack */;
    }

    export class Spheres extends Coarse<Kind.Spheres, CoarseSphereConformation> { }
    export class Gaussians extends Coarse<Kind.Gaussians, CoarseGaussianConformation> { }

    export function areSameChainOperatorGroup(a: Unit, b: Unit) {
        return a.chainGroupId === b.chainGroupId && a.conformation.operator.name === b.conformation.operator.name;
    }
}

export default Unit;