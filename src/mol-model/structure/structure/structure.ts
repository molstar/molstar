/**
 * Copyright (c) 2017-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { IntMap, SortedArray, Iterator, Segmentation, Interval } from '../../../mol-data/int';
import { UniqueArray } from '../../../mol-data/generic';
import { SymmetryOperator } from '../../../mol-math/geometry/symmetry-operator';
import { Model, ElementIndex } from '../model';
import { sort, arraySwap, hash1, sortArray, hashString, hashFnv32a } from '../../../mol-data/util';
import { StructureElement } from './element';
import { Unit } from './unit';
import { StructureLookup3D } from './util/lookup3d';
import { CoarseElements } from '../model/properties/coarse';
import { StructureSubsetBuilder } from './util/subset-builder';
import { InterUnitBonds, computeInterUnitBonds, Bond } from './unit/bonds';
import { StructureSymmetry } from './symmetry';
import { StructureProperties } from './properties';
import { ResidueIndex, ChainIndex, EntityIndex } from '../model/indexing';
import { Carbohydrates } from './carbohydrates/data';
import { computeCarbohydrates } from './carbohydrates/compute';
import { Vec3, Mat4 } from '../../../mol-math/linear-algebra';
import { idFactory } from '../../../mol-util/id-factory';
import { UUID } from '../../../mol-util';
import { CustomProperties } from '../../custom-property';
import { StructureSelection } from '../query/selection';
import { Boundary } from '../../../mol-math/geometry/boundary';
import { ElementSymbol } from '../model/types';
import { CustomStructureProperty } from '../../../mol-model-props/common/custom-structure-property';
import { Trajectory } from '../trajectory';
import { RuntimeContext, Task } from '../../../mol-task';
import { computeStructureBoundary } from './util/boundary';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { IntraUnitBondMapping, getIntraUnitBondMapping, getSerialMapping, SerialMapping } from './mapping';

/** Internal structure state */
type State = {
    parent?: Structure,
    boundary?: Boundary,
    lookup3d?: StructureLookup3D,
    interUnitBonds?: InterUnitBonds,
    dynamicBonds: boolean,
    interBondsValidUnit?: (unit: Unit) => boolean,
    interBondsValidUnitPair?: (structure: Structure, unitA: Unit, unitB: Unit) => boolean,
    unitSymmetryGroups?: ReadonlyArray<Unit.SymmetryGroup>,
    unitSymmetryGroupsIndexMap?: IntMap<number>,
    unitsSortedByVolume?: ReadonlyArray<Unit>;
    carbohydrates?: Carbohydrates,
    models?: ReadonlyArray<Model>,
    model?: Model,
    masterModel?: Model,
    representativeModel?: Model,
    uniqueResidueNames?: Set<string>,
    uniqueElementSymbols?: Set<ElementSymbol>,
    entityIndices?: ReadonlyArray<EntityIndex>,
    uniqueAtomicResidueIndices?: ReadonlyMap<UUID, ReadonlyArray<ResidueIndex>>,
    serialMapping?: SerialMapping,
    intraUnitBondMapping?: IntraUnitBondMapping,
    hashCode: number,
    transformHash: number,
    elementCount: number,
    bondCount: number,
    uniqueElementCount: number,
    atomicResidueCount: number,
    polymerResidueCount: number,
    polymerGapCount: number,
    polymerUnitCount: number,
    coordinateSystem: SymmetryOperator,
    label: string,
    propertyData?: any,
    customProps?: CustomProperties
}

class Structure {
    subsetBuilder(isSorted: boolean) {
        return new StructureSubsetBuilder(this, isSorted);
    }

    /** Count of all elements in the structure, i.e. the sum of the elements in the units */
    get elementCount() {
        return this.state.elementCount;
    }

    /** Count of all bonds (intra- and inter-unit) in the structure */
    get bondCount() {
        if (this.state.bondCount === -1) {
            this.state.bondCount = (this.interUnitBonds.edgeCount / 2) + Bond.getIntraUnitBondCount(this);
        }
        return this.state.bondCount;
    }

    get hasCustomProperties() {
        return !!this.state.customProps && this.state.customProps.all.length > 0;
    }

    get customPropertyDescriptors() {
        if (!this.state.customProps) this.state.customProps = new CustomProperties();
        return this.state.customProps;
    }

    /**
     * Property data unique to this instance of the structure.
     */
    get currentPropertyData() {
        if (!this.state.propertyData) this.state.propertyData = Object.create(null);
        return this.state.propertyData;
    }

    /**
     * Property data of the parent structure if it exists, currentPropertyData otherwise.
     */
    get inheritedPropertyData() {
        return this.parent ? this.parent.currentPropertyData : this.currentPropertyData;
    }

    /** Count of all polymer residues in the structure */
    get polymerResidueCount() {
        if (this.state.polymerResidueCount === -1) {
            this.state.polymerResidueCount = getPolymerResidueCount(this);
        }
        return this.state.polymerResidueCount;
    }

    /** Count of all polymer gaps in the structure */
    get polymerGapCount() {
        if (this.state.polymerGapCount === -1) {
            this.state.polymerGapCount = getPolymerGapCount(this);
        }
        return this.state.polymerGapCount;
    }

    get polymerUnitCount() {
        if (this.state.polymerUnitCount === -1) {
            this.state.polymerUnitCount = getPolymerUnitCount(this);
        }
        return this.state.polymerUnitCount;
    }

    get uniqueElementCount() {
        if (this.state.uniqueElementCount === -1) {
            this.state.uniqueElementCount = getUniqueElementCount(this);
        }
        return this.state.uniqueElementCount;
    }

    get atomicResidueCount() {
        if (this.state.atomicResidueCount === -1) {
            this.state.atomicResidueCount = getAtomicResidueCount(this);
        }
        return this.state.atomicResidueCount;
    }

    /**
     * True if any model the structure is based on is coarse grained.
     * @see Model.isCoarseGrained
     */
    get isCoarseGrained() {
        return this.models.some(m => Model.isCoarseGrained(m));
    }

    get isEmpty() {
        return this.units.length === 0;
    }

    get hashCode() {
        if (this.state.hashCode !== -1) return this.state.hashCode;
        return this.computeHash();
    }

    /** Hash based on all unit.id values in the structure, reflecting the units transformation */
    get transformHash() {
        if (this.state.transformHash !== -1) return this.state.transformHash;
        this.state.transformHash = hashFnv32a(this.units.map(u => u.id));
        return this.state.transformHash;
    }

    private computeHash() {
        let hash = 23;
        for (let i = 0, _i = this.units.length; i < _i; i++) {
            const u = this.units[i];
            hash = (31 * hash + u.id) | 0;
            hash = (31 * hash + SortedArray.hashCode(u.elements)) | 0;
        }
        hash = (31 * hash + this.elementCount) | 0;
        hash = hash1(hash);
        if (hash === -1) hash = 0;
        this.state.hashCode = hash;
        return hash;
    }

    /** Returns a new element location iterator */
    elementLocations(): Iterator<StructureElement.Location> {
        return new Structure.ElementLocationIterator(this);
    }

    /** The parent or itself in case this is the root */
    get root(): Structure {
        return this.state.parent ?? this;
    }

    /** The root/top-most parent or `undefined` in case this is the root */
    get parent() {
        return this.state.parent;
    }

    /**
     * Conformation transformation that was applied to every unit of this structure.
     *
     * Coordinate system applies to the *current* structure only.
     * A parent structure can have a different coordinate system and thefore it has to be composed "manualy"
     * by the consumer.
     */
    get coordinateSystem(): SymmetryOperator {
        return this.state.coordinateSystem;
    }

    get label() {
        return this.state.label;
    }

    get boundary() {
        if (this.state.boundary) return this.state.boundary;
        this.state.boundary = computeStructureBoundary(this);
        return this.state.boundary;
    }

    get lookup3d() {
        if (this.state.lookup3d) return this.state.lookup3d;
        this.state.lookup3d = new StructureLookup3D(this);
        return this.state.lookup3d;
    }

    get interUnitBonds() {
        if (this.state.interUnitBonds) return this.state.interUnitBonds;
        if (this.parent && this.state.dynamicBonds === this.parent.state.dynamicBonds &&
            this.parent.state.interUnitBonds && this.parent.state.interUnitBonds.edgeCount === 0
        ) {
            // no need to compute InterUnitBonds if parent's ones are empty
            this.state.interUnitBonds = new InterUnitBonds(new Map());
        } else {
            this.state.interUnitBonds = computeInterUnitBonds(this, {
                ignoreWater: !this.dynamicBonds,
                ignoreIon: !this.dynamicBonds,
                validUnit: this.state.interBondsValidUnit,
                validUnitPair: this.state.interBondsValidUnitPair,
            });
        }
        return this.state.interUnitBonds;
    }

    get dynamicBonds() {
        return this.state.dynamicBonds;
    }

    get interBondsValidUnit() {
        return this.state.interBondsValidUnit;
    }

    get interBondsValidUnitPair() {
        return this.state.interBondsValidUnitPair;
    }

    get unitSymmetryGroups(): ReadonlyArray<Unit.SymmetryGroup> {
        if (this.state.unitSymmetryGroups) return this.state.unitSymmetryGroups;
        this.state.unitSymmetryGroups = StructureSymmetry.computeTransformGroups(this);
        return this.state.unitSymmetryGroups;
    }

    /** Maps unit.id to index of SymmetryGroup in unitSymmetryGroups array */
    get unitSymmetryGroupsIndexMap(): IntMap<number> {
        if (this.state.unitSymmetryGroupsIndexMap) return this.state.unitSymmetryGroupsIndexMap;
        this.state.unitSymmetryGroupsIndexMap = Unit.SymmetryGroup.getUnitSymmetryGroupsIndexMap(this.unitSymmetryGroups);
        return this.state.unitSymmetryGroupsIndexMap;
    }

    get carbohydrates(): Carbohydrates {
        if (this.state.carbohydrates) return this.state.carbohydrates;
        this.state.carbohydrates = computeCarbohydrates(this);
        return this.state.carbohydrates;
    }

    get models(): ReadonlyArray<Model> {
        if (this.state.models) return this.state.models;
        this.state.models = getModels(this);
        return this.state.models;
    }

    get uniqueResidueNames() {
        return this.state.uniqueResidueNames
            || (this.state.uniqueResidueNames = getUniqueResidueNames(this));
    }

    get uniqueElementSymbols() {
        return this.state.uniqueElementSymbols
            || (this.state.uniqueElementSymbols = getUniqueElementSymbols(this));
    }

    get entityIndices() {
        return this.state.entityIndices
            || (this.state.entityIndices = getEntityIndices(this));
    }

    get uniqueAtomicResidueIndices() {
        return this.state.uniqueAtomicResidueIndices
            || (this.state.uniqueAtomicResidueIndices = getUniqueAtomicResidueIndices(this));
    }

    /** Contains only atomic units */
    get isAtomic() {
        for (const u of this.units) if (!Unit.isAtomic(u)) return false;
        return true;
    }

    /** Contains some atomic units */
    get hasAtomic() {
        for (const u of this.units) if (Unit.isAtomic(u)) return true;
        return false;
    }

    /** Contains only coarse units */
    get isCoarse() {
        for (const u of this.units) if (!Unit.isCoarse(u)) return false;
        return true;
    }

    /** Contains some coarse units */
    get hasCoarse() {
        for (const u of this.units) if (Unit.isCoarse(u)) return true;
        return false;
    }

    /**
     * Provides mapping for serial element indices accross all units.
     *
     * Note that this is especially costly for structures with many units that are grouped
     * into few symmetry groups. Use only when needed and prefer `StructureElement`
     * to address elements in a structure.
     */
    get serialMapping() {
        return this.state.serialMapping || (this.state.serialMapping = getSerialMapping(this));
    }

    get intraUnitBondMapping() {
        return this.state.intraUnitBondMapping || (this.state.intraUnitBondMapping = getIntraUnitBondMapping(this));
    }

    /**
     * If the structure is based on a single model or has a master-/representative-model, return it.
     * Otherwise throw an exception.
     */
    get model(): Model {
        if (this.state.model) return this.state.model;
        if (this.state.representativeModel) return this.state.representativeModel;
        if (this.state.masterModel) return this.state.masterModel;
        const models = this.models;
        if (models.length > 1) {
            throw new Error('The structure is based on multiple models and has neither a master- nor a representative-model.');
        }
        this.state.model = models[0];
        return this.state.model;
    }

    /** The master-model, other models can have bonds to it  */
    get masterModel(): Model | undefined {
        return this.state.masterModel;
    }

    /** A representative model, e.g. the first model of a trajectory */
    get representativeModel(): Model | undefined {
        return this.state.representativeModel;
    }

    hasElement(e: StructureElement.Location) {
        if (!this.unitMap.has(e.unit.id)) return false;
        return SortedArray.has(this.unitMap.get(e.unit.id).elements, e.element);
    }

    getModelIndex(m: Model) {
        return this.models.indexOf(m);
    }

    remapModel(m: Model): Structure {
        const { dynamicBonds, interUnitBonds, parent } = this.state;
        const units: Unit[] = [];
        for (const ug of this.unitSymmetryGroups) {
            const unit = ug.units[0].remapModel(m, dynamicBonds);
            units.push(unit);
            for (let i = 1, il = ug.units.length; i < il; ++i) {
                const u = ug.units[i];
                units.push(u.remapModel(m, dynamicBonds, unit.props));
            }
        }
        return Structure.create(units, {
            parent: parent?.remapModel(m),
            label: this.label,
            interUnitBonds: dynamicBonds ? undefined : interUnitBonds,
            dynamicBonds,
            interBondsValidUnit: this.state.interBondsValidUnit,
            interBondsValidUnitPair: this.state.interBondsValidUnitPair,
            coordinateSystem: this.state.coordinateSystem,
            masterModel: this.state.masterModel,
            representativeModel: this.state.representativeModel,
        });
    }

    private _child: Structure | undefined;
    private _target: Structure | undefined;
    private _proxy: Structure | undefined;

    /**
     * For `structure` with `parent` this returns a proxy that
     * targets `parent` and has `structure` attached as a child.
     */
    asParent(): Structure {
        if (this._proxy) return this._proxy;
        if (this.parent) {
            const p = this.parent.coordinateSystem.isIdentity ? this.parent : Structure.transform(this.parent, this.parent.coordinateSystem.inverse);
            const s = this.coordinateSystem.isIdentity ? p : Structure.transform(p, this.coordinateSystem.matrix);
            this._proxy = new Structure(s.units, s.unitMap, s.unitIndexMap, { ...s.state, dynamicBonds: this.dynamicBonds }, { child: this, target: this.parent });
        } else {
            this._proxy = this;
        }
        return this._proxy;
    }

    get child(): Structure | undefined {
        return this._child;
    }

    /** Get the proxy target. Usefull for equality checks. */
    get target(): Structure {
        return this._target ?? this;
    }

    /**
     * @param units Array of all units in the structure, sorted by unit.id
     * @param unitMap Maps unit.id to index of unit in units array
     * @param unitIndexMap Array of all units in the structure, sorted by unit.id
     */
    constructor(readonly units: ReadonlyArray<Unit>, readonly unitMap: IntMap<Unit>, readonly unitIndexMap: IntMap<number>, private readonly state: State, asParent?: { child: Structure, target: Structure }) {
        // always assign to ensure object shape
        this._child = asParent?.child;
        this._target = asParent?.target;
        this._proxy = undefined;
    }
}

function cmpUnits(units: ArrayLike<Unit>, i: number, j: number) {
    return units[i].id - units[j].id;
}

function getModels(s: Structure) {
    const { units } = s;
    const arr = UniqueArray.create<Model['id'], Model>();
    for (const u of units) {
        UniqueArray.add(arr, u.model.id, u.model);
    }
    return arr.array;
}

function getUniqueResidueNames(s: Structure) {
    const { microheterogeneityCompIds } = StructureProperties.residue;
    const names = new Set<string>();
    const loc = StructureElement.Location.create(s);
    for (const unitGroup of s.unitSymmetryGroups) {
        const unit = unitGroup.units[0];
        // TODO: support coarse unit?
        if (!Unit.isAtomic(unit)) continue;
        const residues = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
        loc.unit = unit;
        while (residues.hasNext) {
            const seg = residues.move();
            loc.element = unit.elements[seg.start];
            const compIds = microheterogeneityCompIds(loc);
            for (const compId of compIds) names.add(compId);
        }
    }
    return names;
}

function getUniqueElementSymbols(s: Structure) {
    const prop = StructureProperties.atom.type_symbol;
    const symbols = new Set<ElementSymbol>();
    const loc = StructureElement.Location.create(s);
    for (const unitGroup of s.unitSymmetryGroups) {
        const unit = unitGroup.units[0];
        if (!Unit.isAtomic(unit)) continue;
        loc.unit = unit;
        for (let i = 0, il = unit.elements.length; i < il; ++i) {
            loc.element = unit.elements[i];
            symbols.add(prop(loc));
        }
    }
    return symbols;
}

function getEntityIndices(structure: Structure): ReadonlyArray<EntityIndex> {
    const { units } = structure;
    const l = StructureElement.Location.create(structure);
    const keys = UniqueArray.create<number, EntityIndex>();

    for (const unit of units) {
        const prop = unit.kind === Unit.Kind.Atomic ? StructureProperties.entity.key : StructureProperties.coarse.entityKey;

        l.unit = unit;
        const elements = unit.elements;

        const chainsIt = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, elements);
        while (chainsIt.hasNext) {
            const chainSegment = chainsIt.move();
            l.element = elements[chainSegment.start];
            const key = prop(l);
            UniqueArray.add(keys, key, key);
        }
    }

    sortArray(keys.array);
    return keys.array;
}

function getUniqueAtomicResidueIndices(structure: Structure): ReadonlyMap<UUID, ReadonlyArray<ResidueIndex>> {
    const map = new Map<UUID, UniqueArray<ResidueIndex, ResidueIndex>>();
    const modelIds: UUID[] = [];

    const unitGroups = structure.unitSymmetryGroups;
    for (const unitGroup of unitGroups) {
        const unit = unitGroup.units[0];
        if (!Unit.isAtomic(unit)) continue;

        let uniqueResidues: UniqueArray<ResidueIndex, ResidueIndex>;
        if (map.has(unit.model.id)) uniqueResidues = map.get(unit.model.id)!;
        else {
            uniqueResidues = UniqueArray.create<ResidueIndex, ResidueIndex>();
            modelIds.push(unit.model.id);
            map.set(unit.model.id, uniqueResidues);
        }

        const residues = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
        while (residues.hasNext) {
            const seg = residues.move();
            UniqueArray.add(uniqueResidues, seg.index, seg.index);
        }
    }

    const ret = new Map<UUID, ReadonlyArray<ResidueIndex>>();
    for (const id of modelIds) {
        const array = map.get(id)!.array;
        sortArray(array);
        ret.set(id, array);
    }
    return ret;
}

function getUniqueElementCount(structure: Structure): number {
    const { unitSymmetryGroups } = structure;
    let uniqueElementCount = 0;
    for (let i = 0, _i = unitSymmetryGroups.length; i < _i; i++) {
        uniqueElementCount += unitSymmetryGroups[i].elements.length;
    }
    return uniqueElementCount;
}

function getPolymerResidueCount(structure: Structure): number {
    const { units } = structure;
    let polymerResidueCount = 0;
    for (let i = 0, _i = units.length; i < _i; i++) {
        polymerResidueCount += units[i].polymerElements.length;
    }
    return polymerResidueCount;
}

function getPolymerGapCount(structure: Structure): number {
    const { units } = structure;
    let polymerGapCount = 0;
    for (let i = 0, _i = units.length; i < _i; i++) {
        polymerGapCount += units[i].gapElements.length / 2;
    }
    return polymerGapCount;
}

function getPolymerUnitCount(structure: Structure): number {
    const { units } = structure;
    let polymerUnitCount = 0;
    for (let i = 0, _i = units.length; i < _i; i++) {
        if (units[i].polymerElements.length > 0) polymerUnitCount += 1;
    }
    return polymerUnitCount;
}

function getAtomicResidueCount(structure: Structure): number {
    const { units } = structure;
    let atomicResidueCount = 0;
    for (let i = 0, _i = units.length; i < _i; i++) {
        const unit = units[i];
        if (!Unit.isAtomic(unit)) continue;
        const { elements, residueIndex } = unit;
        let idx = -1;
        let prevIdx = -1;
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            idx = residueIndex[elements[j]];
            if (idx !== prevIdx) {
                atomicResidueCount += 1;
                prevIdx = idx;
            }
        }
    }
    return atomicResidueCount;
}

namespace Structure {
    export const Empty = create([]);

    export interface Props {
        parent?: Structure
        interUnitBonds?: InterUnitBonds
        /**
         * Ensure bonds are recalculated upon model changes.
         * Also enables calculation of inter-unit bonds in water molecules.
         */
        dynamicBonds?: boolean,
        interBondsValidUnit?: (unit: Unit) => boolean,
        interBondsValidUnitPair?: (structure: Structure, unitA: Unit, unitB: Unit) => boolean,
        coordinateSystem?: SymmetryOperator
        label?: string
        /** Master model for structures of a protein model and multiple ligand models */
        masterModel?: Model
        /** Representative model for structures of a model trajectory */
        representativeModel?: Model
    }

    /** Represents a single structure */
    export interface Loci {
        readonly kind: 'structure-loci',
        readonly structure: Structure,
    }
    export function Loci(structure: Structure): Loci {
        return { kind: 'structure-loci', structure };
    }

    export function toStructureElementLoci(structure: Structure): StructureElement.Loci {
        const elements: StructureElement.Loci['elements'][0][] = [];
        for (const unit of structure.units) {
            elements.push({ unit, indices: Interval.ofBounds(0, unit.elements.length) });
        }
        return StructureElement.Loci(structure, elements);
    }

    export function toSubStructureElementLoci(parent: Structure, structure: Structure): StructureElement.Loci {
        return StructureSelection.toLociWithSourceUnits(StructureSelection.Singletons(parent, structure));
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'structure-loci';
    }

    export function areLociEqual(a: Loci, b: Loci) {
        return a.structure === b.structure;
    }

    export function isLociEmpty(loci: Loci) {
        return loci.structure.isEmpty;
    }

    export function remapLoci(loci: Loci, structure: Structure) {
        if (structure === loci.structure) return loci;
        return Loci(structure);
    }

    export function create(units: ReadonlyArray<Unit>, props: Props = {}): Structure {
        // init units
        const unitMap = IntMap.Mutable<Unit>();
        const unitIndexMap = IntMap.Mutable<number>();
        let elementCount = 0;
        let isSorted = true;
        let lastId = units.length > 0 ? units[0].id : 0;
        for (let i = 0, _i = units.length; i < _i; i++) {
            const u = units[i];
            unitMap.set(u.id, u);
            elementCount += u.elements.length;
            if (u.id < lastId) isSorted = false;
            lastId = u.id;
        }
        if (!isSorted) sort(units, 0, units.length, cmpUnits, arraySwap);
        for (let i = 0, _i = units.length; i < _i; i++) {
            unitIndexMap.set(units[i].id, i);
        }

        // initial state
        const state: State = {
            hashCode: -1,
            transformHash: -1,
            elementCount,
            bondCount: -1,
            uniqueElementCount: -1,
            atomicResidueCount: -1,
            polymerResidueCount: -1,
            polymerGapCount: -1,
            polymerUnitCount: -1,
            dynamicBonds: false,
            coordinateSystem: SymmetryOperator.Default,
            label: ''
        };

        // handle props
        if (props.parent) state.parent = props.parent.parent || props.parent;
        if (props.interUnitBonds) state.interUnitBonds = props.interUnitBonds;

        if (props.interBondsValidUnit) state.interBondsValidUnit = props.interBondsValidUnit;
        else if (props.parent) state.interBondsValidUnit = props.parent.interBondsValidUnit;

        if (props.interBondsValidUnitPair) state.interBondsValidUnitPair = props.interBondsValidUnitPair;
        else if (props.parent) state.interBondsValidUnitPair = props.parent.interBondsValidUnitPair;

        if (props.dynamicBonds) state.dynamicBonds = props.dynamicBonds;
        else if (props.parent) state.dynamicBonds = props.parent.dynamicBonds;

        if (props.coordinateSystem) state.coordinateSystem = props.coordinateSystem;
        else if (props.parent) state.coordinateSystem = props.parent.coordinateSystem;

        if (props.label) state.label = props.label;
        else if (props.parent) state.label = props.parent.label;

        if (props.masterModel) state.masterModel = props.masterModel;
        else if (props.parent) state.masterModel = props.parent.masterModel;

        if (props.representativeModel) state.representativeModel = props.representativeModel;
        else if (props.parent) state.representativeModel = props.parent.representativeModel;

        return new Structure(units, unitMap, unitIndexMap, state);
    }

    export async function ofTrajectory(trajectory: Trajectory, ctx: RuntimeContext): Promise<Structure> {
        if (trajectory.frameCount === 0) return Empty;

        const units: Unit[] = [];

        let first: Model | undefined = void 0;
        let count = 0;
        for (let i = 0, il = trajectory.frameCount; i < il; ++i) {
            const frame = await Task.resolveInContext(trajectory.getFrameAtIndex(i), ctx);
            if (!first) first = frame;
            const structure = ofModel(frame);
            for (let j = 0, jl = structure.units.length; j < jl; ++j) {
                const u = structure.units[j];
                const invariantId = u.invariantId + count;
                const chainGroupId = u.chainGroupId + count;
                const newUnit = Unit.create(units.length, invariantId, chainGroupId, u.traits, u.kind, u.model, u.conformation.operator, u.elements);
                units.push(newUnit);
            }
            count = units.length;
        }

        return create(units, { representativeModel: first!, label: first!.label });
    }

    /**
     * Construct a Structure from a model.
     *
     * Generally, a single unit corresponds to a single chain, with the exception
     * of consecutive "single atom chains" with same entity_id and same auth_asym_id.
     */
    export function ofModel(model: Model, props: Props = {}): Structure {
        const chains = model.atomicHierarchy.chainAtomSegments;
        const { index } = model.atomicHierarchy;
        const { auth_asym_id } = model.atomicHierarchy.chains;
        const { atomicChainOperatorMappinng } = model;
        const builder = new StructureBuilder({ label: model.label, ...props });

        for (let c = 0 as ChainIndex; c < chains.count; c++) {
            const operator = atomicChainOperatorMappinng.get(c) || SymmetryOperator.Default;
            const start = chains.offsets[c];

            let isMultiChain = false;
            let isWater = false;

            if (isWaterChain(model, c)) {
                isWater = true;
                // merge consecutive water chains
                while (c + 1 < chains.count && isWaterChain(model, c + 1 as ChainIndex)) {
                    const op1 = atomicChainOperatorMappinng.get(c);
                    const op2 = atomicChainOperatorMappinng.get(c + 1 as ChainIndex);
                    if (op1 !== op2) break;

                    isMultiChain = true;
                    c++;
                }
            } else {
                // merge consecutive "single atom chains" with same entity_id and auth_asym_id
                while (c + 1 < chains.count
                    && chains.offsets[c + 1] - chains.offsets[c] === 1
                    && chains.offsets[c + 2] - chains.offsets[c + 1] === 1
                ) {
                    const e1 = index.getEntityFromChain(c);
                    const e2 = index.getEntityFromChain(c + 1 as ChainIndex);
                    if (e1 !== e2) break;

                    const a1 = auth_asym_id.value(c);
                    const a2 = auth_asym_id.value(c + 1);
                    if (a1 !== a2) break;

                    const op1 = atomicChainOperatorMappinng.get(c);
                    const op2 = atomicChainOperatorMappinng.get(c + 1 as ChainIndex);
                    if (op1 !== op2) break;

                    isMultiChain = true;
                    c++;
                }
            }

            const elements = SortedArray.ofBounds(start as ElementIndex, chains.offsets[c + 1] as ElementIndex);

            let traits = Unit.Trait.None;
            if (isMultiChain) traits |= Unit.Trait.MultiChain;
            if (isWater) traits |= Unit.Trait.Water;
            if (Model.isCoarseGrained(model)) traits |= Unit.Trait.CoarseGrained;

            builder.addUnit(Unit.Kind.Atomic, model, operator, elements, traits);
        }

        const cs = model.coarseHierarchy;
        if (cs.isDefined) {
            if (cs.spheres.count > 0) {
                addCoarseUnits(builder, model, model.coarseHierarchy.spheres, Unit.Kind.Spheres);
            }
            if (cs.gaussians.count > 0) {
                addCoarseUnits(builder, model, model.coarseHierarchy.gaussians, Unit.Kind.Gaussians);
            }
        }

        return builder.getStructure();
    }

    function isWaterChain(model: Model, chainIndex: ChainIndex) {
        const e = model.atomicHierarchy.index.getEntityFromChain(chainIndex);
        return model.entities.data.type.value(e) === 'water';
    }

    function addCoarseUnits(builder: StructureBuilder, model: Model, elements: CoarseElements, kind: Unit.Kind) {
        const { chainElementSegments } = elements;
        for (let cI = 0; cI < chainElementSegments.count; cI++) {
            const elements = SortedArray.ofBounds<ElementIndex>(chainElementSegments.offsets[cI], chainElementSegments.offsets[cI + 1]);
            builder.addUnit(kind, model, SymmetryOperator.Default, elements, Unit.Trait.None);
        }
    }

    export function transform(s: Structure, transform: Mat4) {
        if (Mat4.isIdentity(transform)) return s;
        if (!Mat4.isRotationAndTranslation(transform, SymmetryOperator.RotationTranslationEpsilon)) throw new Error('Only rotation/translation combination can be applied.');

        const units: Unit[] = [];
        for (const u of s.units) {
            const old = u.conformation.operator;
            const op = SymmetryOperator.create(old.name, transform, old);
            units.push(u.applyOperator(u.id, op));
        }

        const cs = s.coordinateSystem;
        const newCS = SymmetryOperator.compose(SymmetryOperator.create(cs.name, transform, cs), cs);
        return create(units, { parent: s, coordinateSystem: newCS });
    }

     export function instances(s: Structure, transforms: Mat4[]) {
        for (const t of transforms) {
            if (!Mat4.isRotationAndTranslation(t, SymmetryOperator.RotationTranslationEpsilon)) {
                throw new Error('Only rotation/translation combination can be applied.');
            }
        }

        const units: Unit[] = [];
        let id = 0;
        for (const u of s.units) {
            const old = u.conformation.operator;
            for (const transform of transforms) {
                if (Mat4.isIdentity(transform)) {
                    units.push(u);
                    continue;
                }
                const op = SymmetryOperator.create(old.name, transform, old);
                units.push(u.applyOperator(id++, op));
            }
        }

        return create(units, { parent: s });
    }


    export class StructureBuilder {
        private units: Unit[] = [];
        private invariantId = idFactory();

        private chainGroupId = -1;
        private inChainGroup = false;

        private p = Vec3();
        private singleElementUnits = new Map<string, Unit>();

        beginChainGroup() {
            this.chainGroupId++;
            this.inChainGroup = true;
        }

        endChainGroup() {
            this.inChainGroup = false;
        }

        addUnit(kind: Unit.Kind, model: Model, operator: SymmetryOperator, elements: StructureElement.Set, traits: Unit.Traits, invariantId?: number): Unit {
            if (invariantId === undefined) invariantId = this.invariantId();
            const chainGroupId = this.inChainGroup ? this.chainGroupId : ++this.chainGroupId;
            const unit = Unit.create(this.units.length, invariantId, chainGroupId, traits, kind, model, operator, elements);
            return this.add(unit);
        }

        copyUnit(unit: Unit, options?: Unit.GetCopyOptions & { invariantId?: number }): Unit {
            let invariantId = options?.invariantId;
            if (invariantId === undefined) invariantId = this.invariantId();
            const chainGroupId = this.inChainGroup ? this.chainGroupId : ++this.chainGroupId;
            const newUnit = unit.getCopy(this.units.length, invariantId, chainGroupId, options);
            return this.add(newUnit);
        }

        private add(unit: Unit) {
            // this is to avoid adding many identical single atom units,
            // for example, from 'degenerate' spacegroups
            // - Diamond (COD 9008564)
            if (unit.elements.length === 1) {
                unit.conformation.position(unit.elements[0], this.p);
                const hash = [unit.invariantId, this.p[0], this.p[1], this.p[2]].join('|');
                if (this.singleElementUnits.has(hash)) {
                    return this.singleElementUnits.get(hash)!;
                } else {
                    this.singleElementUnits.set(hash, unit);
                }
            }
            this.units.push(unit);
            return unit;
        }

        addWithOperator(unit: Unit, operator: SymmetryOperator, dontCompose = false): Unit {
            return this.add(unit.applyOperator(this.units.length, operator, dontCompose));
        }

        getStructure(): Structure {
            return create(this.units, this.props);
        }

        get isEmpty() {
            return this.units.length === 0;
        }

        constructor(private props: Props = {}) {

        }
    }

    export function Builder(props: Props = {}) {
        return new StructureBuilder(props);
    }

    export function hashCode(s: Structure) {
        return s.hashCode;
    }

    /** Hash based on all unit.model conformation values in the structure */
    export function conformationHash(s: Structure) {
        return hashString(s.units.map(u => Unit.conformationId(u)).join('|'));
    }

    // TODO: there should be a version that properly supports partitioned units
    export function areUnitIdsEqual(a: Structure, b: Structure) {
        if (a === b) return true;

        if (a.elementCount !== b.elementCount) return false;
        const len = a.units.length;
        if (len !== b.units.length) return false;

        for (let i = 0; i < len; i++) {
            if (a.units[i].id !== b.units[i].id) return false;
        }

        return true;
    }

    export function areUnitIdsAndIndicesEqual(a: Structure, b: Structure) {
        if (a === b) return true;
        if (!areUnitIdsEqual(a, b)) return false;

        for (let i = 0, il = a.units.length; i < il; i++) {
            if (!SortedArray.areEqual(a.units[i].elements, b.units[i].elements)) return false;
        }
        return true;
    }

    export function areHierarchiesEqual(a: Structure, b: Structure) {
        if (a.hashCode !== b.hashCode) return false;

        const len = a.models.length;
        if (len !== b.models.length) return false;

        for (let i = 0; i < len; i++) {
            if (!Model.areHierarchiesEqual(a.models[i], b.models[i])) return false;
        }
        return true;
    }

    export function areEquivalent(a: Structure, b: Structure) {
        return a === b || (
            a.hashCode === b.hashCode &&
            StructureSymmetry.areTransformGroupsEquivalent(a.unitSymmetryGroups, b.unitSymmetryGroups)
        );
    }

    /** Check if the structures or their parents are equivalent */
    export function areRootsEquivalent(a: Structure, b: Structure) {
        return areEquivalent(a.root, b.root);
    }

    /** Check if the structures or their parents are equal */
    export function areRootsEqual(a: Structure, b: Structure) {
        return a.root === b.root;
    }

    export class ElementLocationIterator implements Iterator<StructureElement.Location> {
        private current: StructureElement.Location;
        private unitIndex = 0;
        private elements: StructureElement.Set;
        private maxIdx = 0;
        private idx = -1;

        hasNext: boolean;
        move(): StructureElement.Location {
            this.advance();
            this.current.element = this.elements[this.idx];
            return this.current;
        }

        private advance() {
            if (this.idx < this.maxIdx) {
                this.idx++;

                if (this.idx === this.maxIdx) this.hasNext = this.unitIndex + 1 < this.structure.units.length;
                return;
            }

            this.idx = 0;
            this.unitIndex++;
            if (this.unitIndex >= this.structure.units.length) {
                this.hasNext = false;
                return;
            }

            this.current.unit = this.structure.units[this.unitIndex];
            this.elements = this.current.unit.elements;
            this.maxIdx = this.elements.length - 1;
            if (this.maxIdx === 0) {
                this.hasNext = this.unitIndex + 1 < this.structure.units.length;
            }
        }

        constructor(private structure: Structure) {
            this.current = StructureElement.Location.create(structure);
            this.hasNext = structure.elementCount > 0;
            if (this.hasNext) {
                this.elements = structure.units[0].elements;
                this.maxIdx = this.elements.length - 1;
                this.current.unit = structure.units[0];
            }
        }
    }

    const distVec = Vec3();
    function unitElementMinDistance(unit: Unit, p: Vec3, eRadius: number) {
        const { elements, conformation: c } = unit, dV = distVec;
        let minD = Number.MAX_VALUE;
        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            const d = Vec3.distance(p, c.position(e, dV)) - eRadius - c.r(e);
            if (d < minD) minD = d;
        }
        return minD;
    }

    export function minDistanceToPoint(s: Structure, point: Vec3, radius: number) {
        const { units } = s;
        let minD = Number.MAX_VALUE;
        for (let i = 0, _i = units.length; i < _i; i++) {
            const unit = units[i];
            const d = unitElementMinDistance(unit, point, radius);
            if (d < minD) minD = d;
        }
        return minD;
    }

    const distPivot = Vec3();
    export function distance(a: Structure, b: Structure) {
        if (a.elementCount === 0 || b.elementCount === 0) return 0;

        const { units } = a;
        let minD = Number.MAX_VALUE;

        for (let i = 0, _i = units.length; i < _i; i++) {
            const unit = units[i];
            const { elements, conformation: c } = unit;
            for (let i = 0, _i = elements.length; i < _i; i++) {
                const e = elements[i];
                const d = minDistanceToPoint(b, c.position(e, distPivot), c.r(e));
                if (d < minD) minD = d;
            }
        }
        return minD;
    }

    export function elementDescription(s: Structure) {
        return s.elementCount === 1 ? '1 element' : `${s.elementCount} elements`;
    }

    export function validUnitPair(s: Structure, a: Unit, b: Unit) {
        return s.masterModel
            ? a.model === b.model || a.model === s.masterModel || b.model === s.masterModel
            : a.model === b.model;
    }

    export interface EachUnitPairProps {
        maxRadius: number
        validUnit: (unit: Unit) => boolean
        validUnitPair: (unitA: Unit, unitB: Unit) => boolean
    }

    /**
     * Iterate over all unit pairs of a structure and invokes callback for valid units
     * and unit pairs if their boundaries are within a max distance.
     */
    export function eachUnitPair(structure: Structure, callback: (unitA: Unit, unitB: Unit) => void, props: EachUnitPairProps) {
        const { maxRadius, validUnit, validUnitPair } = props;
        if (!structure.units.some(u => validUnit(u))) return;

        const lookup = structure.lookup3d;
        const imageCenter = Vec3();

        for (const unitA of structure.units) {
            if (!validUnit(unitA)) continue;

            const bs = unitA.boundary.sphere;
            Vec3.transformMat4(imageCenter, bs.center, unitA.conformation.operator.matrix);
            const closeUnits = lookup.findUnitIndices(imageCenter[0], imageCenter[1], imageCenter[2], bs.radius + maxRadius);
            for (let i = 0; i < closeUnits.count; i++) {
                const unitB = structure.units[closeUnits.indices[i]];
                if (unitA.id >= unitB.id) continue;
                if (!validUnit(unitB) || !validUnitPair(unitA, unitB)) continue;

                if (unitB.elements.length >= unitA.elements.length) callback(unitA, unitB);
                else callback(unitB, unitA);
            }
        }
    }

    export interface ForEachAtomicHierarchyElementParams {
        // Called for 1st element of each chain
        // Note that chains can be split, meaning each chain would be called multiple times.
        chain?: (e: StructureElement.Location<Unit.Atomic>) => void,
        // Called for 1st element of each residue
        residue?: (e: StructureElement.Location<Unit.Atomic>) => void,
        // Called for each element
        atom?: (e: StructureElement.Location<Unit.Atomic>) => void,
    };

    export function eachAtomicHierarchyElement(structure: Structure, { chain, residue, atom }: ForEachAtomicHierarchyElementParams) {
        const l = StructureElement.Location.create<Unit.Atomic>(structure);
        for (const unit of structure.units) {
            if (unit.kind !== Unit.Kind.Atomic) continue;

            l.unit = unit;

            const { elements } = unit;
            const chainsIt = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, elements);
            const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);

            while (chainsIt.hasNext) {
                const chainSegment = chainsIt.move();

                if (chain) {
                    l.element = elements[chainSegment.start];
                    chain(l);
                }

                if (!residue && !atom) continue;

                residuesIt.setSegment(chainSegment);
                while (residuesIt.hasNext) {
                    const residueSegment = residuesIt.move();

                    if (residue) {
                        l.element = elements[residueSegment.start];
                        residue(l);
                    }

                    if (!atom) continue;

                    for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
                        l.element = elements[j];
                        atom(l);
                    }
                }
            }
        }
    }

    //

    export const DefaultSizeThresholds = {
        /** Must be lower to be small */
        smallResidueCount: 10,
        /** Must be lower to be medium */
        mediumResidueCount: 5000,
        /** Must be lower to be large (big ribosomes like 4UG0 should still be `large`) */
        largeResidueCount: 30000,
        /**
         * Structures above `largeResidueCount` are consider huge when they have
         * a `highSymmetryUnitCount` or gigantic when not
         */
        highSymmetryUnitCount: 10,
        /** Fiber-like structure are consider small when below this */
        fiberResidueCount: 15,
    };
    export type SizeThresholds = typeof DefaultSizeThresholds

    function getPolymerSymmetryGroups(structure: Structure) {
        return structure.unitSymmetryGroups.filter(ug => ug.units[0].polymerElements.length > 0);
    }

    /**
     * Try to match fiber-like structures like 6nk4
     */
    function isFiberLike(structure: Structure, thresholds: SizeThresholds) {
        const polymerSymmetryGroups = getPolymerSymmetryGroups(structure);
        return (
            polymerSymmetryGroups.length === 1 &&
            polymerSymmetryGroups[0].units.length > 2 &&
            polymerSymmetryGroups[0].units[0].polymerElements.length < thresholds.fiberResidueCount
        );
    }

    function hasHighSymmetry(structure: Structure, thresholds: SizeThresholds) {
        const polymerSymmetryGroups = getPolymerSymmetryGroups(structure);
        return (
            polymerSymmetryGroups.length >= 1 &&
            polymerSymmetryGroups[0].units.length > thresholds.highSymmetryUnitCount
        );
    }

    export enum Size { Small, Medium, Large, Huge, Gigantic }

    /**
     * @param residueCountFactor - modifies the threshold counts, useful when estimating
     *                             the size of a structure comprised of multiple models
     */
    export function getSize(structure: Structure, thresholds: Partial<SizeThresholds> = {}, residueCountFactor = 1): Size {
        const t = { ...DefaultSizeThresholds, ...thresholds };
        if (structure.polymerResidueCount >= t.largeResidueCount * residueCountFactor) {
            if (hasHighSymmetry(structure, t)) {
                return Size.Huge;
            } else {
                return Size.Gigantic;
            }
        } else if (isFiberLike(structure, t)) {
            return Size.Small;
        } else if (structure.polymerResidueCount < t.smallResidueCount * residueCountFactor) {
            return Size.Small;
        } else if (structure.polymerResidueCount < t.mediumResidueCount * residueCountFactor) {
            return Size.Medium;
        } else {
            return Size.Large;
        }
    }

    //

    export type Index = number;
    export const Index = CustomStructureProperty.createSimple<Index>('index', 'root');

    export type MaxIndex = number;
    export const MaxIndex = CustomStructureProperty.createSimple<MaxIndex>('max_index', 'root');

    const PrincipalAxesProp = '__PrincipalAxes__';
    export function getPrincipalAxes(structure: Structure): PrincipalAxes {
        if (structure.currentPropertyData[PrincipalAxesProp]) return structure.currentPropertyData[PrincipalAxesProp];
        const principalAxes = StructureElement.Loci.getPrincipalAxes(Structure.toStructureElementLoci(structure));
        structure.currentPropertyData[PrincipalAxesProp] = principalAxes;
        return principalAxes;
    }
}

export { Structure };