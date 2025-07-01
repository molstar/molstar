/**
 * Copyright (c) 2017-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SymmetryOperator } from '../../../mol-math/geometry/symmetry-operator';
import { Model } from '../model';
import { GridLookup3D, Lookup3D, Spacegroup } from '../../../mol-math/geometry';
import { IntraUnitBonds, computeIntraUnitBonds } from './unit/bonds';
import { CoarseElements, CoarseSphereConformation, CoarseGaussianConformation } from '../model/properties/coarse';
import { BitFlags } from '../../../mol-util';
import { UnitRings } from './unit/rings';
import { StructureElement } from './element';
import { ChainIndex, ResidueIndex, ElementIndex } from '../model/indexing';
import { IntMap, SortedArray, Segmentation } from '../../../mol-data/int';
import { hash2, hashFnv32a } from '../../../mol-data/util';
import { getAtomicPolymerElements, getCoarsePolymerElements, getAtomicGapElements, getCoarseGapElements, getNucleotideElements, getProteinElements } from './util/polymer';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { getPrincipalAxes } from './util/principal-axes';
import { Boundary, getBoundary, getFastBoundary } from '../../../mol-math/geometry/boundary';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { IndexPairBonds } from '../../../mol-model-formats/structure/property/bonds/index-pair';
import { ElementSetIntraBondCache } from './unit/bonds/element-set-intra-bond-cache';
import { ModelSymmetry } from '../../../mol-model-formats/structure/property/symmetry';
import { getResonance, UnitResonance } from './unit/resonance';

/**
 * A building block of a structure that corresponds to an atomic or
 * a coarse grained representation 'conveniently grouped together'.
 */
type Unit = Unit.Atomic | Unit.Spheres | Unit.Gaussians

namespace Unit {
    export const enum Kind { Atomic, Spheres, Gaussians }

    // To use with isolatedModules
    export enum Kinds { Atomic = Kind.Atomic, Spheres = Kind.Spheres, Gaussians = Kind.Gaussians }

    export function isAtomic(u: Unit): u is Atomic { return u.kind === Kind.Atomic; }
    export function isCoarse(u: Unit): u is Spheres | Gaussians { return u.kind === Kind.Spheres || u.kind === Kind.Gaussians; }
    export function isSpheres(u: Unit): u is Spheres { return u.kind === Kind.Spheres; }
    export function isGaussians(u: Unit): u is Gaussians { return u.kind === Kind.Gaussians; }

    export function create<K extends Kind>(id: number, invariantId: number, chainGroupId: number, traits: Traits, kind: Kind, model: Model, operator: SymmetryOperator, elements: StructureElement.Set, props?: K extends Kind.Atomic ? AtomicProperties : CoarseProperties): Unit {
        switch (kind) {
            case Kind.Atomic: return new Atomic(id, invariantId, chainGroupId, traits, model, elements, SymmetryOperator.createMapping(operator, model.atomicConformation), props ?? AtomicProperties());
            case Kind.Spheres: return createCoarse(id, invariantId, chainGroupId, traits, model, Kind.Spheres, elements, SymmetryOperator.createMapping(operator, model.coarseConformation.spheres, getSphereRadiusFunc(model)), props ?? CoarseProperties());
            case Kind.Gaussians: return createCoarse(id, invariantId, chainGroupId, traits, model, Kind.Gaussians, elements, SymmetryOperator.createMapping(operator, model.coarseConformation.gaussians, getGaussianRadiusFunc(model)), props ?? CoarseProperties());
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

    function getTransformHash(units: Unit[]) {
        const ids: number[] = [];
        for (let i = 0, _i = units.length; i < _i; i++) {
            ids.push(units[i].id);
        }
        return hashFnv32a(ids);
    }

    export function SymmetryGroup(units: Unit[]) {
        const props: {
            unitIndexMap?: IntMap<number>
        } = {};

        return {
            elements: units[0].elements,
            units,
            get unitIndexMap() {
                if (props.unitIndexMap) return props.unitIndexMap;
                props.unitIndexMap = getUnitIndexMap(units);
                return props.unitIndexMap;
            },
            hashCode: hashUnit(units[0]),
            transformHash: getTransformHash(units)
        };
    }

    export namespace SymmetryGroup {
        export function areInvariantElementsEqual(a: SymmetryGroup, b: SymmetryGroup) {
            if (a.hashCode !== b.hashCode) return false;
            return SortedArray.areEqual(a.elements, b.elements);
        }

        export function getUnitSymmetryGroupsIndexMap(symmetryGroups: ReadonlyArray<Unit.SymmetryGroup>): IntMap<number> {
            const unitSymmetryGroupsIndexMap = IntMap.Mutable<number>();
            for (let i = 0, il = symmetryGroups.length; i < il; ++i) {
                const sg = symmetryGroups[i];
                for (let j = 0, jl = sg.units.length; j < jl; ++j) {
                    unitSymmetryGroupsIndexMap.set(sg.units[j].id, i);
                }
            }
            return unitSymmetryGroupsIndexMap;
        }
    }

    export function conformationId(unit: Unit) {
        return Unit.isAtomic(unit) ? unit.model.atomicConformation.id : unit.model.coarseConformation.id;
    }

    export function hashUnit(u: Unit) {
        return hash2(u.invariantId, SortedArray.hashCode(u.elements));
    }

    export type Traits = BitFlags<Trait>
    export enum Trait {
        None = 0x0,
        MultiChain = 0x1,
        Partitioned = 0x2,
        FastBoundary = 0x4,
        Water = 0x8,
        CoarseGrained = 0x10,
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
        readonly props: BaseProperties,
        readonly transientCache: Map<any, any>,

        getChild(elements: StructureElement.Set): Unit,
        getCopy(id: number, invariantId: number, chainGroupId: number, options?: GetCopyOptions): Unit,
        applyOperator(id: number, operator: SymmetryOperator, dontCompose?: boolean /* = false */): Unit,
        remapModel(model: Model, dynamicBonds: boolean): Unit,

        readonly boundary: Boundary
        readonly lookup3d: Lookup3D<StructureElement.UnitIndex>
        readonly polymerElements: SortedArray<ElementIndex>
        readonly gapElements: SortedArray<ElementIndex>
        /**
         * From mmCIF/IHM schema: `_ihm_model_representation_details.model_object_primitive`.
         */
        readonly objectPrimitive: mmCIF_Schema['ihm_model_representation_details']['model_object_primitive']['T']
    }

    export interface GetCopyOptions {
        propagateTransientCache?: boolean;
    }

    interface BaseProperties {
        boundary?: Boundary
        lookup3d?: Lookup3D<StructureElement.UnitIndex>
        principalAxes?: PrincipalAxes
        polymerElements?: SortedArray<ElementIndex>
        gapElements?: SortedArray<ElementIndex>
    }


    function BaseProperties(): BaseProperties {
        return {};
    }

    function getSphereRadiusFunc(model: Model) {
        const r = model.coarseConformation.spheres.radius;
        return (i: ElementIndex) => r[i];
    }

    function getGaussianRadiusFunc(_model: Model) {
        // TODO: compute radius for gaussians
        return (i: ElementIndex) => 0;
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

        readonly props: AtomicProperties;

        private _transientCache: Map<any, any> | undefined = undefined;
        get transientCache() {
            if (this._transientCache === void 0) this._transientCache = new Map<any, any>();
            return this._transientCache;
        }

        getChild(elements: StructureElement.Set): Unit {
            if (elements.length === this.elements.length) return this;
            return new Atomic(this.id, this.invariantId, this.chainGroupId, this.traits, this.model, elements, this.conformation, AtomicProperties());
        }

        getCopy(id: number, invariantId: number, chainGroupId: number, options?: GetCopyOptions): Unit {
            const unit = new Atomic(id, invariantId, chainGroupId, this.traits, this.model, this.elements, this.conformation, this.props);
            if (options?.propagateTransientCache) {
                unit._transientCache = this._transientCache;
            }
            return unit;
        }

        applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
            const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
            return new Atomic(id, this.invariantId, this.chainGroupId, this.traits, this.model, this.elements, SymmetryOperator.createMapping(op, this.model.atomicConformation, this.conformation.r), this.props);
        }

        remapModel(model: Model, dynamicBonds: boolean, props?: AtomicProperties) {
            if (!props) {
                props = {
                    ...this.props,
                    bonds: dynamicBonds && !this.props.bonds?.props?.canRemap
                        ? undefined
                        : tryRemapBonds(this, this.props.bonds, model, dynamicBonds)
                };
                if (!Unit.isSameConformation(this, model)) {
                    props.boundary = undefined;
                    props.lookup3d = undefined;
                    props.principalAxes = undefined;
                }
            }

            let operator = this.conformation.operator;
            const symmetry = ModelSymmetry.Provider.get(model);
            if (operator.spgrOp !== -1 && symmetry && symmetry !== ModelSymmetry.Provider.get(this.model)) {
                const [i, j, k] = operator.hkl;
                const { toFractional } = symmetry.spacegroup.cell;
                const ref = Vec3.transformMat4(Vec3(), Model.getCenter(model), toFractional);
                operator = Spacegroup.getSymmetryOperatorRef(symmetry.spacegroup, operator.spgrOp, i, j, k, ref);
            }

            const conformation = (this.model.atomicConformation !== model.atomicConformation || operator !== this.conformation.operator)
                ? SymmetryOperator.createMapping<ElementIndex>(operator, model.atomicConformation)
                : this.conformation;
            return new Atomic(this.id, this.invariantId, this.chainGroupId, this.traits, model, this.elements, conformation, props);
        }

        get boundary() {
            if (this.props.boundary) return this.props.boundary;
            const { x, y, z } = this.model.atomicConformation;
            this.props.boundary = Traits.is(this.traits, Trait.FastBoundary)
                ? getFastBoundary({ x, y, z, indices: this.elements })
                : getBoundary({ x, y, z, indices: this.elements });
            return this.props.boundary;
        }

        get lookup3d() {
            if (this.props.lookup3d) return this.props.lookup3d;
            const { x, y, z } = this.model.atomicConformation;
            this.props.lookup3d = GridLookup3D({ x, y, z, indices: this.elements }, this.boundary);
            return this.props.lookup3d;
        }

        get principalAxes() {
            if (this.props.principalAxes) return this.props.principalAxes;
            this.props.principalAxes = getPrincipalAxes(this);
            return this.props.principalAxes;
        }

        get bonds() {
            if (this.props.bonds) return this.props.bonds;

            const cache = ElementSetIntraBondCache.get(this.model);
            let bonds = cache.get(this.elements);
            if (!bonds) {
                bonds = computeIntraUnitBonds(this);
                if (bonds.props?.cacheable) {
                    cache.set(this.elements, bonds);
                }
            }
            this.props.bonds = bonds;
            return this.props.bonds;
        }

        get rings() {
            if (this.props.rings) return this.props.rings;
            this.props.rings = UnitRings.create(this);
            return this.props.rings;
        }

        get resonance() {
            if (this.props.resonance) return this.props.resonance;
            this.props.resonance = getResonance(this);
            return this.props.resonance;
        }

        get polymerElements() {
            if (this.props.polymerElements) return this.props.polymerElements;
            this.props.polymerElements = getAtomicPolymerElements(this);
            return this.props.polymerElements;
        }

        get gapElements() {
            if (this.props.gapElements) return this.props.gapElements;
            this.props.gapElements = getAtomicGapElements(this);
            return this.props.gapElements;
        }

        get nucleotideElements() {
            if (this.props.nucleotideElements) return this.props.nucleotideElements;
            this.props.nucleotideElements = getNucleotideElements(this);
            return this.props.nucleotideElements;
        }

        get proteinElements() {
            if (this.props.proteinElements) return this.props.proteinElements;
            this.props.proteinElements = getProteinElements(this);
            return this.props.proteinElements;
        }

        get residueCount(): number {
            if (this.props.residueCount !== undefined) return this.props.residueCount;

            let residueCount = 0;
            const residueIt = Segmentation.transientSegments(this.model.atomicHierarchy.residueAtomSegments, this.elements);
            while (residueIt.hasNext) {
                residueIt.move();
                residueCount += 1;
            }

            this.props.residueCount = residueCount;
            return this.props.residueCount!;
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
        bonds?: IntraUnitBonds
        rings?: UnitRings
        resonance?: UnitResonance
        nucleotideElements?: SortedArray<ElementIndex>
        proteinElements?: SortedArray<ElementIndex>
        residueCount?: number
    }

    function AtomicProperties(): AtomicProperties {
        return BaseProperties();
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

        readonly props: CoarseProperties;

        private _transientCache: Map<any, any> | undefined = undefined;
        get transientCache() {
            if (this._transientCache === void 0) this._transientCache = new Map<any, any>();
            return this._transientCache;
        }

        getChild(elements: StructureElement.Set): Unit {
            if (elements.length === this.elements.length) return this as any as Unit; // lets call this an ugly temporary hack
            return createCoarse(this.id, this.invariantId, this.chainGroupId, this.traits, this.model, this.kind, elements, this.conformation, CoarseProperties());
        }

        getCopy(id: number, invariantId: number, chainGroupId: number, options?: GetCopyOptions): Unit {
            const unit = createCoarse(id, invariantId, chainGroupId, this.traits, this.model, this.kind, this.elements, this.conformation, this.props);
            if (options?.propagateTransientCache) {
                unit._transientCache = this._transientCache;
            }
            return unit;
        }

        applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
            const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
            return createCoarse(id, this.invariantId, this.chainGroupId, this.traits, this.model, this.kind, this.elements, SymmetryOperator.createMapping(op, this.getCoarseConformation(), this.conformation.r), this.props);
        }

        remapModel(model: Model, dynamicBonds: boolean, props?: CoarseProperties): Unit.Spheres | Unit.Gaussians {
            const coarseConformation = this.getCoarseConformation();
            const modelCoarseConformation = getCoarseConformation(this.kind, model);

            if (!props) {
                props = { ...this.props };
                if (!Unit.isSameConformation(this as Unit.Spheres | Unit.Gaussians, model)) { // TODO get rid of casting
                    props.boundary = undefined;
                    props.lookup3d = undefined;
                    props.principalAxes = undefined;
                }
            }

            const conformation = coarseConformation !== modelCoarseConformation
                ? SymmetryOperator.createMapping(this.conformation.operator, modelCoarseConformation, this.kind === Unit.Kind.Spheres ? getSphereRadiusFunc(model) : getGaussianRadiusFunc(model))
                : this.conformation;
            return new Coarse(this.id, this.invariantId, this.chainGroupId, this.traits, model, this.kind, this.elements, conformation, props) as Unit.Spheres | Unit.Gaussians; // TODO get rid of casting
        }

        get boundary() {
            if (this.props.boundary) return this.props.boundary;
            // TODO: support sphere radius?
            const { x, y, z } = this.getCoarseConformation();
            this.props.boundary = Traits.is(this.traits, Trait.FastBoundary)
                ? getFastBoundary({ x, y, z, indices: this.elements })
                : getBoundary({ x, y, z, indices: this.elements });
            return this.props.boundary;
        }

        get lookup3d() {
            if (this.props.lookup3d) return this.props.lookup3d;
            // TODO: support sphere radius?
            const { x, y, z } = this.getCoarseConformation();
            this.props.lookup3d = GridLookup3D({ x, y, z, indices: this.elements }, this.boundary);
            return this.props.lookup3d;
        }

        get principalAxes() {
            if (this.props.principalAxes) return this.props.principalAxes;
            this.props.principalAxes = getPrincipalAxes(this as Unit.Spheres | Unit.Gaussians); // TODO get rid of casting
            return this.props.principalAxes;
        }

        get polymerElements() {
            if (this.props.polymerElements) return this.props.polymerElements;
            this.props.polymerElements = getCoarsePolymerElements(this as Unit.Spheres | Unit.Gaussians); // TODO get rid of casting
            return this.props.polymerElements;
        }

        get gapElements() {
            if (this.props.gapElements) return this.props.gapElements;
            this.props.gapElements = getCoarseGapElements(this as Unit.Spheres | Unit.Gaussians); // TODO get rid of casting
            return this.props.gapElements;
        }

        private getCoarseConformation() {
            return getCoarseConformation(this.kind, this.model);
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

    function getCoarseConformation(kind: Kind, model: Model) {
        return kind === Kind.Spheres ? model.coarseConformation.spheres : model.coarseConformation.gaussians;
    }

    interface CoarseProperties extends BaseProperties { }

    function CoarseProperties(): CoarseProperties {
        return BaseProperties();
    }

    export class Spheres extends Coarse<Kind.Spheres, CoarseSphereConformation> { }
    export class Gaussians extends Coarse<Kind.Gaussians, CoarseGaussianConformation> { }

    function createCoarse<K extends Kind.Gaussians | Kind.Spheres>(id: number, invariantId: number, chainGroupId: number, traits: Traits, model: Model, kind: K, elements: StructureElement.Set, conformation: SymmetryOperator.ArrayMapping<ElementIndex>, props: CoarseProperties): K extends Kind.Spheres ? Spheres : Gaussians {
        return new Coarse(id, invariantId, chainGroupId, traits, model, kind, elements, conformation, props) as any; // lets call this an ugly temporary hack
    }

    export function areSameChainOperatorGroup(a: Unit, b: Unit) {
        return a.chainGroupId === b.chainGroupId && a.conformation.operator.name === b.conformation.operator.name;
    }

    export function areOperatorsEqual(a: Unit, b: Unit) {
        return Mat4.areEqual(a.conformation.operator.matrix, b.conformation.operator.matrix, 1e-6);
    }

    export function areConformationsEqual(a: Unit, b: Unit) {
        if (a === b) return true;
        if (!SortedArray.areEqual(a.elements, b.elements)) return false;
        return isSameConformation(a, b.model);
    }

    function tryRemapBonds(a: Atomic, old: IntraUnitBonds | undefined, model: Model, dynamicBonds: boolean) {
        // TODO: should include additional checks?

        if (!old) return void 0;
        if (a.model.atomicConformation.id === model.atomicConformation.id) return old;

        const oldIndex = IndexPairBonds.Provider.get(a.model);
        if (oldIndex) {
            const newIndex = IndexPairBonds.Provider.get(model);
            // TODO: check the actual indices instead of just reference equality?
            if (!newIndex || oldIndex === newIndex) return old;
            return void 0;
        }

        if (old.props?.canRemap || !dynamicBonds) {
            return old;
        }
        return isSameConformation(a, model) ? old : void 0;
    }

    export function isSameConformation(u: Unit, model: Model) {
        const coordsHistory = Model.CoordinatesHistory.get(Model.getRoot(model));
        if (coordsHistory) return coordsHistory.areEqual(u.elements, u.kind, model);

        const xs = u.elements;
        const { x: xa, y: ya, z: za } = u.conformation.coordinates;
        const { x: xb, y: yb, z: zb } = getModelConformationOfKind(u.kind, model);

        for (let i = 0, _i = xs.length; i < _i; i++) {
            const u = xs[i];
            if (xa[u] !== xb[u] || ya[u] !== yb[u] || za[u] !== zb[u]) return false;
        }

        return true;
    }

    export function getModelConformationOfKind(kind: Unit.Kind, model: Model) {
        return kind === Kind.Atomic ? model.atomicConformation :
            kind === Kind.Spheres ? model.coarseConformation.spheres :
                model.coarseConformation.gaussians;
    }

    export function getConformation(u: Unit) {
        return getModelConformationOfKind(u.kind, u.model);
    }

    export function getModelHierarchyOfKind(kind: Unit.Kind, model: Model) {
        return kind === Kind.Atomic ? model.atomicHierarchy :
            kind === Kind.Spheres ? model.coarseHierarchy.spheres :
                model.coarseHierarchy.gaussians;
    }

    export function getHierarchy(u: Unit) {
        return getModelHierarchyOfKind(u.kind, u.model);
    }
}

export { Unit };