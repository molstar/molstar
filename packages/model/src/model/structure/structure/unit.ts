/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { Model } from '../model.js';
import { GridLookup3D, type Lookup3D, Spacegroup } from '@molstar/core/math/geometry';
import { type IntraUnitBonds, computeIntraUnitBonds } from './unit/bonds.js';
import { VdwRadius } from '../model/properties/atomic.js';
import type {
  CoarseElements,
  CoarseSphereConformation,
  CoarseGaussianConformation,
} from '../model/properties/coarse.js';
import { BitFlags } from '@molstar/core/util';
import { UnitRings } from './unit/rings.js';
import type { StructureElement } from './element.js';
import type { ChainIndex, ResidueIndex, ElementIndex } from '../model/indexing.js';
import { IntMap, SortedArray, Segmentation } from '@molstar/core/data/int';
import { hash2, hashFnv32a } from '@molstar/core/data/util';
import {
  getAtomicPolymerElements,
  getCoarsePolymerElements,
  getAtomicGapElements,
  getCoarseGapElements,
  getNucleotideElements,
  getProteinElements,
} from './util/polymer.js';
import type { mmCIF_Schema } from '@molstar/io/reader/cif/schema/mmcif';
import type { PrincipalAxes } from '@molstar/core/math/linear-algebra/matrix/principal-axes';
import { getPrincipalAxes } from './util/principal-axes.js';
import { type Boundary, getBoundary, getFastBoundary } from '@molstar/core/math/geometry/boundary';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import { ElementSetIntraBondCache } from './unit/bonds/element-set-intra-bond-cache.js';
import { subsetIntraUnitBonds } from './unit/bonds/subset.js';
import { ModelSymmetry } from '@molstar/model/formats/structure/property/symmetry';
import { getResonance, type UnitResonance } from './unit/resonance.js';
import { BondProvider } from './unit/bonds/provider.js';
import { Structure } from './structure.js';

/**
 * A building block of a structure that corresponds to an atomic or
 * a coarse grained representation 'conveniently grouped together'.
 */
type Unit = Unit.Atomic | Unit.Spheres | Unit.Gaussians;

namespace Unit {
  export const enum Kind {
    Atomic,
    Spheres,
    Gaussians,
  }

  // To use with isolatedModules
  export enum Kinds {
    Atomic = Kind.Atomic,
    Spheres = Kind.Spheres,
    Gaussians = Kind.Gaussians,
  }

  export function isAtomic(u: Unit): u is Atomic {
    return u.kind === Kind.Atomic;
  }
  export function isCoarse(u: Unit): u is Spheres | Gaussians {
    return u.kind === Kind.Spheres || u.kind === Kind.Gaussians;
  }
  export function isSpheres(u: Unit): u is Spheres {
    return u.kind === Kind.Spheres;
  }
  export function isGaussians(u: Unit): u is Gaussians {
    return u.kind === Kind.Gaussians;
  }

  export function create<K extends Kind>(
    id: number,
    invariantId: number,
    chainGroupId: number,
    traits: Traits,
    kind: Kind,
    model: Model,
    operator: SymmetryOperator,
    elements: StructureElement.Set,
    props?: K extends Kind.Atomic ? AtomicProperties : CoarseProperties,
    /** Atomic units only: where its bonds come from; the built-in bonds by default (bond-providers plan §2.1). */
    bondProvider: BondProvider.Instance = BondProvider.DefaultInstance,
  ): Unit {
    switch (kind) {
      case Kind.Atomic:
        return new Atomic(
          id,
          invariantId,
          chainGroupId,
          traits,
          model,
          elements,
          SymmetryOperator.createMapping(operator, model.atomicConformation, getAtomicRadiusFunc(model)),
          props ?? AtomicProperties(),
          undefined,
          bondProvider,
        );
      case Kind.Spheres:
        return createCoarse(
          id,
          invariantId,
          chainGroupId,
          traits,
          model,
          Kind.Spheres,
          elements,
          SymmetryOperator.createMapping(operator, model.coarseConformation.spheres, getSphereRadiusFunc(model)),
          props ?? CoarseProperties(),
        );
      case Kind.Gaussians:
        return createCoarse(
          id,
          invariantId,
          chainGroupId,
          traits,
          model,
          Kind.Gaussians,
          elements,
          SymmetryOperator.createMapping(operator, model.coarseConformation.gaussians, getGaussianRadiusFunc(model)),
          props ?? CoarseProperties(),
        );
    }
  }

  /** A group of units that differ only by symmetry operators. */
  export type SymmetryGroup = {
    readonly elements: StructureElement.Set;
    readonly units: ReadonlyArray<Unit>;
    /** Maps unit.id to index of unit in units array */
    readonly unitIndexMap: IntMap<number>;
    /** Hash based on unit.invariantId which is the same for all units in the group */
    readonly hashCode: number;
    /** Hash based on all unit.id values in the group, reflecting the units transformation*/
    readonly transformHash: number;
  };

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
      unitIndexMap?: IntMap<number>;
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
      transformHash: getTransformHash(units),
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

  export type Traits = BitFlags<Trait>;
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
    readonly id: number;
    /** invariant ID stays the same even if the Operator/conformation changes. */
    readonly invariantId: number;
    readonly chainGroupId: number;
    readonly traits: Traits;
    readonly elements: StructureElement.Set;
    readonly model: Model;
    readonly conformation: SymmetryOperator.ArrayMapping<ElementIndex>;
    /**
     * The unit `getChild` made this one from, flattened to the root as `Structure.parent` is: same kind and model,
     * a superset of the elements, possibly another operator. Only its static props may be used (unit.md §2).
     */
    readonly parent: Unit | undefined;
    /** Static props: shared by every unit with these elements and model (symmetry copies), read-only (unit.md §1). */
    readonly props: Readonly<BaseProperties>;
    /** Dynamic data of this unit instance, never inherited by children (unit.md §1). */
    readonly transientCache: Map<any, any>;

    getChild(elements: StructureElement.Set): Unit;
    getCopy(id: number, invariantId: number, chainGroupId: number, options?: GetCopyOptions): Unit;
    applyOperator(id: number, operator: SymmetryOperator, dontCompose?: boolean /* = false */): Unit;
    remapModel(model: Model, dynamicBonds: boolean): Unit;

    readonly boundary: Boundary;
    readonly lookup3d: Lookup3D<StructureElement.UnitIndex>;
    readonly polymerElements: SortedArray<ElementIndex>;
    readonly gapElements: SortedArray<ElementIndex>;
    /**
     * From mmCIF/IHM schema: `_ihm_model_representation_details.model_object_primitive`.
     */
    readonly objectPrimitive: mmCIF_Schema['ihm_model_representation_details']['model_object_primitive']['T'];
  }

  export interface GetCopyOptions {
    propagateTransientCache?: boolean;
  }

  interface BaseProperties {
    boundary?: Boundary;
    lookup3d?: Lookup3D<StructureElement.UnitIndex>;
    principalAxes?: PrincipalAxes;
    polymerElements?: SortedArray<ElementIndex>;
    gapElements?: SortedArray<ElementIndex>;
  }

  function BaseProperties(): BaseProperties {
    return {};
  }

  function getAtomicRadiusFunc(model: Model) {
    const type_symbol = model.atomicHierarchy.atoms.type_symbol.value;
    return (i: ElementIndex) => VdwRadius(type_symbol(i));
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

    readonly parent: Atomic | undefined;
    /**
     * Where the unit's intra-unit bonds come from: fixed at construction and passed on to children, copies and
     * remapped units (bond-providers plan §2.1). Not a static prop: units with the same props may differ in it only
     * when made by `Unit.create`, which never shares props.
     */
    readonly bondProvider: BondProvider.Instance;
    private readonly _props: AtomicProperties;
    get props(): Readonly<AtomicProperties> {
      return this._props;
    }

    private _transientCache: Map<any, any> | undefined = undefined;
    get transientCache() {
      if (this._transientCache === void 0) this._transientCache = new Map<any, any>();
      return this._transientCache;
    }

    getChild(elements: StructureElement.Set): Unit {
      if (elements.length === this.elements.length) return this;
      return new Atomic(
        this.id,
        this.invariantId,
        this.chainGroupId,
        this.traits,
        this.model,
        elements,
        this.conformation,
        AtomicProperties(),
        this.parent ?? this,
        this.bondProvider,
      );
    }

    getCopy(id: number, invariantId: number, chainGroupId: number, options?: GetCopyOptions): Unit {
      const unit = new Atomic(
        id,
        invariantId,
        chainGroupId,
        this.traits,
        this.model,
        this.elements,
        this.conformation,
        this._props,
        this.parent,
        this.bondProvider,
      );
      if (options?.propagateTransientCache) {
        unit._transientCache = this._transientCache;
      }
      return unit;
    }

    applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
      const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
      return new Atomic(
        id,
        this.invariantId,
        this.chainGroupId,
        this.traits,
        this.model,
        this.elements,
        SymmetryOperator.createMapping(op, this.model.atomicConformation, this.conformation.r),
        this._props,
        this.parent,
        this.bondProvider,
      );
    }

    remapModel(model: Model, dynamicBonds: boolean, props?: AtomicProperties) {
      if (!props) {
        props = {
          ...this._props,
          bonds:
            dynamicBonds && !this._props.bonds?.props?.canRemap
              ? undefined
              : tryRemapBonds(this, this._props.bonds, model, dynamicBonds),
        };
        // rings and resonance are derived from the bonds (unit.md §3)
        if (!props.bonds) {
          props.rings = undefined;
          props.resonance = undefined;
        }
        // a frame made whole has the ranges found on whole coordinates (model periodic.ts)
        if (this.model.atomicRanges !== model.atomicRanges) {
          props.polymerElements = undefined;
          props.gapElements = undefined;
        }
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

      const conformation =
        this.model.atomicConformation !== model.atomicConformation || operator !== this.conformation.operator
          ? SymmetryOperator.createMapping<ElementIndex>(operator, model.atomicConformation, getAtomicRadiusFunc(model))
          : this.conformation;
      return new Atomic(
        this.id,
        this.invariantId,
        this.chainGroupId,
        this.traits,
        model,
        this.elements,
        conformation,
        props,
        // the parent belongs to the old model (unit.md §3)
        model === this.model ? this.parent : undefined,
        // model-independent: a frame keeps the structure's provider (bond-providers plan §4)
        this.bondProvider,
      );
    }

    get boundary() {
      if (this._props.boundary) return this._props.boundary;
      const fast = Traits.is(this.traits, Trait.FastBoundary);
      const { x, y, z } = this.model.atomicConformation;
      const radius = fast ? undefined : Model.getAtomicRadii(this.model);
      this._props.boundary = fast
        ? getFastBoundary({ x, y, z, radius, indices: this.elements })
        : getBoundary({ x, y, z, radius, indices: this.elements });
      return this._props.boundary;
    }

    get lookup3d() {
      if (this._props.lookup3d) return this._props.lookup3d;
      const { x, y, z } = this.model.atomicConformation;
      this._props.lookup3d = GridLookup3D({ x, y, z, indices: this.elements }, this.boundary);
      return this._props.lookup3d;
    }

    get principalAxes() {
      if (this._props.principalAxes) return this._props.principalAxes;
      this._props.principalAxes = getPrincipalAxes(this);
      return this._props.principalAxes;
    }

    get bonds(): IntraUnitBonds {
      if (this._props.bonds) return this._props.bonds;

      // a child restricts its parent's bonds if they exist, and never computes them (unit.md §4), unless its provider
      // asks for the root's bonds (`childBonds: 'root'`, bond-providers plan §2.3); only bonds of the same provider
      const { parent, bondProvider } = this;
      const canRestrict = !!parent && parent.bondProvider === bondProvider && canRestrictParentBonds(parent, this);
      const parentBonds: IntraUnitBonds | undefined = canRestrict
        ? bondProvider.childBonds === 'root'
          ? parent!.bonds
          : parent!.props.bonds
        : undefined;
      let bonds: IntraUnitBonds | undefined =
        parentBonds && subsetIntraUnitBonds(parentBonds, parent!.elements, this.elements);
      if (!bonds) {
        const cacheKey = intraCacheKey(bondProvider);
        const cache = cacheKey !== undefined ? ElementSetIntraBondCache.get(this.model, cacheKey) : undefined;
        bonds = cache?.get(this.elements);
        if (!bonds) {
          bonds = bondProvider.intra ? bondProvider.intra(this, intraContext(this)) : computeIntraUnitBonds(this);
          if (cache && bonds.props?.cacheable) {
            cache.set(this.elements, bonds);
          }
        }
      }
      this._props.bonds = bonds;
      return this._props.bonds;
    }

    get rings() {
      if (this._props.rings) return this._props.rings;
      this._props.rings = UnitRings.create(this);
      return this._props.rings;
    }

    get resonance() {
      if (this._props.resonance) return this._props.resonance;
      this._props.resonance = getResonance(this);
      return this._props.resonance;
    }

    get polymerElements() {
      if (this._props.polymerElements) return this._props.polymerElements;
      this._props.polymerElements = getAtomicPolymerElements(this);
      return this._props.polymerElements;
    }

    get gapElements() {
      if (this._props.gapElements) return this._props.gapElements;
      this._props.gapElements = getAtomicGapElements(this);
      return this._props.gapElements;
    }

    get nucleotideElements() {
      if (this._props.nucleotideElements) return this._props.nucleotideElements;
      this._props.nucleotideElements = getNucleotideElements(this);
      return this._props.nucleotideElements;
    }

    get proteinElements() {
      if (this._props.proteinElements) return this._props.proteinElements;
      this._props.proteinElements = getProteinElements(this);
      return this._props.proteinElements;
    }

    get residueCount(): number {
      if (this._props.residueCount !== undefined) return this._props.residueCount;

      let residueCount = 0;
      const residueIt = Segmentation.transientSegments(this.model.atomicHierarchy.residueAtomSegments, this.elements);
      while (residueIt.hasNext) {
        residueIt.move();
        residueCount += 1;
      }

      this._props.residueCount = residueCount;
      return this._props.residueCount!;
    }

    getResidueIndex(elementIndex: StructureElement.UnitIndex) {
      return this.residueIndex[this.elements[elementIndex]];
    }

    /** Check if any element of the unit belongs to the given residue. */
    hasResidueIndex(residueIndex: ResidueIndex) {
      const { offsets } = this.model.atomicHierarchy.residueAtomSegments;
      const start = offsets[residueIndex];
      const end = offsets[residueIndex + 1];
      return SortedArray.hasRange(this.elements, start, end - 1);
    }

    constructor(
      id: number,
      invariantId: number,
      chainGroupId: number,
      traits: Traits,
      model: Model,
      elements: StructureElement.Set,
      conformation: SymmetryOperator.ArrayMapping<ElementIndex>,
      props: AtomicProperties,
      parent?: Atomic,
      bondProvider: BondProvider.Instance = BondProvider.DefaultInstance,
    ) {
      this.id = id;
      this.invariantId = invariantId;
      this.chainGroupId = chainGroupId;
      this.traits = traits;
      this.model = model;
      this.elements = elements;
      this.conformation = conformation;

      this.residueIndex = model.atomicHierarchy.residueAtomSegments.index;
      this.chainIndex = model.atomicHierarchy.chainAtomSegments.index;
      this._props = props;
      this.parent = parent;
      this.bondProvider = bondProvider;
    }
  }

  /**
   * The provider's cache key under its name, so two providers never share a cache by choosing the same key; the
   * built-in bonds keep the model's original cache (bond-providers plan §2.1).
   */
  function intraCacheKey(instance: BondProvider.Instance): string | undefined {
    if (BondProvider.isDefault(instance)) return BondProvider.DefaultName;
    if (instance.cacheKey === undefined) return undefined;
    return `${instance.provider.name}:${instance.cacheKey}`;
  }

  function intraContext(unit: Atomic): BondProvider.IntraContext {
    return {
      builtIn: (props) => computeIntraUnitBonds(unit, props),
      modelStructure: () => modelStructureOf(unit.model),
    };
  }

  /**
   * A model as a structure with the built-in bonds, for providers (`IntraContext.modelStructure`), kept with the
   * model's other derived data. Its units are the default provider's, so reading its bonds never reaches back into a
   * provider.
   */
  function modelStructureOf(model: Model): Structure {
    const data = model._dynamicPropertyData;
    if (!data.BondProviderModelStructure) {
      data.BondProviderModelStructure = Structure.ofModel(model, { bondProvider: BondProvider.DefaultInstance });
    }
    return data.BondProviderModelStructure;
  }

  /**
   * Index-pair bonds may be listed per operator (`molstar_bond_site` with operator ids), so a unit's bonds then depend
   * on its operator: a child restricts the bonds of a parent with another operator only without index pairs (unit.md §4).
   */
  function canRestrictParentBonds(parent: Atomic, child: Atomic) {
    if (parent.conformation.operator.key === child.conformation.operator.key) return true;
    return !IndexPairBonds.Provider.get(child.model);
  }

  interface AtomicProperties extends BaseProperties {
    bonds?: IntraUnitBonds;
    rings?: UnitRings;
    resonance?: UnitResonance;
    nucleotideElements?: SortedArray<ElementIndex>;
    proteinElements?: SortedArray<ElementIndex>;
    residueCount?: number;
  }

  function AtomicProperties(): AtomicProperties {
    return BaseProperties();
  }

  class Coarse<K extends Kind.Gaussians | Kind.Spheres, C extends CoarseSphereConformation | CoarseGaussianConformation>
    implements Base
  {
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

    readonly parent: Unit.Spheres | Unit.Gaussians | undefined;
    private readonly _props: CoarseProperties;
    get props(): Readonly<CoarseProperties> {
      return this._props;
    }

    private _transientCache: Map<any, any> | undefined = undefined;
    get transientCache() {
      if (this._transientCache === void 0) this._transientCache = new Map<any, any>();
      return this._transientCache;
    }

    getChild(elements: StructureElement.Set): Unit {
      if (elements.length === this.elements.length) return this as any as Unit; // lets call this an ugly temporary hack
      return createCoarse(
        this.id,
        this.invariantId,
        this.chainGroupId,
        this.traits,
        this.model,
        this.kind,
        elements,
        this.conformation,
        CoarseProperties(),
        this.parent ?? (this as Unit.Spheres | Unit.Gaussians), // TODO get rid of casting
      );
    }

    getCopy(id: number, invariantId: number, chainGroupId: number, options?: GetCopyOptions): Unit {
      const unit = createCoarse(
        id,
        invariantId,
        chainGroupId,
        this.traits,
        this.model,
        this.kind,
        this.elements,
        this.conformation,
        this._props,
        this.parent,
      );
      if (options?.propagateTransientCache) {
        unit._transientCache = this._transientCache;
      }
      return unit;
    }

    applyOperator(id: number, operator: SymmetryOperator, dontCompose = false): Unit {
      const op = dontCompose ? operator : SymmetryOperator.compose(this.conformation.operator, operator);
      return createCoarse(
        id,
        this.invariantId,
        this.chainGroupId,
        this.traits,
        this.model,
        this.kind,
        this.elements,
        SymmetryOperator.createMapping(op, this.getCoarseConformation(), this.conformation.r),
        this._props,
        this.parent,
      );
    }

    remapModel(model: Model, dynamicBonds: boolean, props?: CoarseProperties): Unit.Spheres | Unit.Gaussians {
      const coarseConformation = this.getCoarseConformation();
      const modelCoarseConformation = getCoarseConformation(this.kind, model);

      if (!props) {
        props = { ...this._props };
        if (!Unit.isSameConformation(this as Unit.Spheres | Unit.Gaussians, model)) {
          // TODO get rid of casting
          props.boundary = undefined;
          props.lookup3d = undefined;
          props.principalAxes = undefined;
        }
      }

      const conformation =
        coarseConformation !== modelCoarseConformation
          ? SymmetryOperator.createMapping(
              this.conformation.operator,
              modelCoarseConformation,
              this.kind === Unit.Kind.Spheres ? getSphereRadiusFunc(model) : getGaussianRadiusFunc(model),
            )
          : this.conformation;
      return new Coarse(
        this.id,
        this.invariantId,
        this.chainGroupId,
        this.traits,
        model,
        this.kind,
        this.elements,
        conformation,
        props,
        model === this.model ? this.parent : undefined,
      ) as Unit.Spheres | Unit.Gaussians; // TODO get rid of casting
    }

    get boundary() {
      if (this._props.boundary) return this._props.boundary;
      const { x, y, z } = this.getCoarseConformation();
      const radius = this.getCoarseRadii();
      this._props.boundary = Traits.is(this.traits, Trait.FastBoundary)
        ? getFastBoundary({ x, y, z, radius, indices: this.elements })
        : getBoundary({ x, y, z, radius, indices: this.elements });
      return this._props.boundary;
    }

    get lookup3d() {
      if (this._props.lookup3d) return this._props.lookup3d;
      // TODO: support sphere radius?
      const { x, y, z } = this.getCoarseConformation();
      this._props.lookup3d = GridLookup3D({ x, y, z, indices: this.elements }, this.boundary);
      return this._props.lookup3d;
    }

    get principalAxes() {
      if (this._props.principalAxes) return this._props.principalAxes;
      this._props.principalAxes = getPrincipalAxes(this as Unit.Spheres | Unit.Gaussians); // TODO get rid of casting
      return this._props.principalAxes;
    }

    get polymerElements() {
      if (this._props.polymerElements) return this._props.polymerElements;
      this._props.polymerElements = getCoarsePolymerElements(this as Unit.Spheres | Unit.Gaussians); // TODO get rid of casting
      return this._props.polymerElements;
    }

    get gapElements() {
      if (this._props.gapElements) return this._props.gapElements;
      this._props.gapElements = getCoarseGapElements(this as Unit.Spheres | Unit.Gaussians); // TODO get rid of casting
      return this._props.gapElements;
    }

    private getCoarseConformation() {
      return getCoarseConformation(this.kind, this.model);
    }

    private getCoarseRadii() {
      return getCoarseRadii(this.kind, this.model);
    }

    constructor(
      id: number,
      invariantId: number,
      chainGroupId: number,
      traits: Traits,
      model: Model,
      kind: K,
      elements: StructureElement.Set,
      conformation: SymmetryOperator.ArrayMapping<ElementIndex>,
      props: CoarseProperties,
      parent?: Unit.Spheres | Unit.Gaussians,
    ) {
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
      this.coarseConformation = (
        kind === Kind.Spheres ? model.coarseConformation.spheres : model.coarseConformation.gaussians
      ) as C;
      this._props = props;
      this.parent = parent;
    }
  }

  function getCoarseConformation(kind: Kind, model: Model) {
    return kind === Kind.Spheres ? model.coarseConformation.spheres : model.coarseConformation.gaussians;
  }

  function getCoarseRadii(kind: Kind, model: Model) {
    return kind === Kind.Spheres ? model.coarseConformation.spheres.radius : undefined; // Zero radius for gaussians
  }

  interface CoarseProperties extends BaseProperties {}

  function CoarseProperties(): CoarseProperties {
    return BaseProperties();
  }

  export class Spheres extends Coarse<Kind.Spheres, CoarseSphereConformation> {}
  export class Gaussians extends Coarse<Kind.Gaussians, CoarseGaussianConformation> {}

  function createCoarse<K extends Kind.Gaussians | Kind.Spheres>(
    id: number,
    invariantId: number,
    chainGroupId: number,
    traits: Traits,
    model: Model,
    kind: K,
    elements: StructureElement.Set,
    conformation: SymmetryOperator.ArrayMapping<ElementIndex>,
    props: CoarseProperties,
    parent?: Unit.Spheres | Unit.Gaussians,
  ): K extends Kind.Spheres ? Spheres : Gaussians {
    return new Coarse(id, invariantId, chainGroupId, traits, model, kind, elements, conformation, props, parent) as any; // lets call this an ugly temporary hack
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
    return kind === Kind.Atomic
      ? model.atomicConformation
      : kind === Kind.Spheres
        ? model.coarseConformation.spheres
        : model.coarseConformation.gaussians;
  }

  export function getConformation(u: Unit) {
    return getModelConformationOfKind(u.kind, u.model);
  }

  export function getModelHierarchyOfKind(kind: Unit.Kind, model: Model) {
    return kind === Kind.Atomic
      ? model.atomicHierarchy
      : kind === Kind.Spheres
        ? model.coarseHierarchy.spheres
        : model.coarseHierarchy.gaussians;
  }

  export function getHierarchy(u: Unit) {
    return getModelHierarchyOfKind(u.kind, u.model);
  }
}

export { Unit };
