/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { IntMap, SortedArray, Iterator, Segmentation, Interval, OrderedSet } from '../../../mol-data/int';
import { UniqueArray } from '../../../mol-data/generic';
import { SymmetryOperator } from '../../../mol-math/geometry/symmetry-operator';
import { Model, ElementIndex } from '../model';
import { sort, arraySwap, hash1, sortArray, hashString, hashFnv32a } from '../../../mol-data/util';
import StructureElement from './element';
import Unit from './unit';
import { StructureLookup3D } from './util/lookup3d';
import { CoarseElements } from '../model/properties/coarse';
import { StructureSubsetBuilder } from './util/subset-builder';
import { InterUnitBonds, computeInterUnitBonds, Bond } from './unit/bonds';
import StructureSymmetry from './symmetry';
import StructureProperties from './properties';
import { ResidueIndex, ChainIndex, EntityIndex } from '../model/indexing';
import { Carbohydrates } from './carbohydrates/data';
import { computeCarbohydrates } from './carbohydrates/compute';
import { Vec3, Mat4 } from '../../../mol-math/linear-algebra';
import { idFactory } from '../../../mol-util/id-factory';
import { GridLookup3D } from '../../../mol-math/geometry';
import { UUID } from '../../../mol-util';
import { CustomProperties } from '../../custom-property';
import { AtomicHierarchy } from '../model/properties/atomic';
import { StructureSelection } from '../query/selection';
import { getBoundary } from '../../../mol-math/geometry/boundary';
import { ElementSymbol } from '../model/types';
import { CustomStructureProperty } from '../../../mol-model-props/common/custom-structure-property';

class Structure {
    /** Maps unit.id to unit */
    readonly unitMap: IntMap<Unit>;
    /** Maps unit.id to index of unit in units array */
    readonly unitIndexMap: IntMap<number>;
    /** Array of all units in the structure, sorted by unit.id */
    readonly units: ReadonlyArray<Unit>;

    private _props: {
        parent?: Structure,
        lookup3d?: StructureLookup3D,
        interUnitBonds?: InterUnitBonds,
        unitSymmetryGroups?: ReadonlyArray<Unit.SymmetryGroup>,
        unitSymmetryGroupsIndexMap?: IntMap<number>,
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
        hashCode: number,
        /** Hash based on all unit.id values in the structure, reflecting the units transformation */
        transformHash: number,
        elementCount: number,
        bondCount: number,
        uniqueElementCount: number,
        atomicResidueCount: number,
        polymerResidueCount: number,
        polymerUnitCount: number,
        coordinateSystem: SymmetryOperator,
        label: string,
        propertyData?: any,
        customProps?: CustomProperties
    } = {
        hashCode: -1,
        transformHash: -1,
        elementCount: -1,
        bondCount: -1,
        uniqueElementCount: -1,
        atomicResidueCount: -1,
        polymerResidueCount: -1,
        polymerUnitCount: -1,
        coordinateSystem: SymmetryOperator.Default,
        label: ''
    };

    subsetBuilder(isSorted: boolean) {
        return new StructureSubsetBuilder(this, isSorted);
    }

    /** Count of all elements in the structure, i.e. the sum of the elements in the units */
    get elementCount() {
        return this._props.elementCount;
    }

    /** Count of all bonds (intra- and inter-unit) in the structure */
    get bondCount() {
        if (!this._props.bondCount) {
            this._props.bondCount = this.interUnitBonds.edgeCount + Bond.getIntraUnitBondCount(this);
        }
        return this._props.bondCount;
    }

    get hasCustomProperties() {
        return !!this._props.customProps && this._props.customProps.all.length > 0;
    }

    get customPropertyDescriptors() {
        if (!this._props.customProps) this._props.customProps = new CustomProperties();
        return this._props.customProps;
    }

    /**
     * Property data unique to this instance of the structure.
     */
    get currentPropertyData() {
        if (!this._props.propertyData) this._props.propertyData = Object.create(null);
        return this._props.propertyData;
    }

    /**
     * Property data of the parent structure if it exists, currentPropertyData otherwise.
     */
    get inheritedPropertyData() {
        return this.parent ? this.parent.currentPropertyData : this.currentPropertyData;
    }

    /** Count of all polymer residues in the structure */
    get polymerResidueCount() {
        if (this._props.polymerResidueCount === -1) {
            this._props.polymerResidueCount = getPolymerResidueCount(this);
        }
        return this._props.polymerResidueCount;
    }

    get polymerUnitCount() {
        if (this._props.polymerUnitCount === -1) {
            this._props.polymerUnitCount = getPolymerUnitCount(this);
        }
        return this._props.polymerUnitCount;
    }

    get uniqueElementCount() {
        if (this._props.uniqueElementCount === -1) {
            this._props.uniqueElementCount = getUniqueElementCount(this);
        }
        return this._props.uniqueElementCount;
    }

    get atomicResidueCount() {
        if (this._props.atomicResidueCount === -1) {
            this._props.atomicResidueCount = getAtomicResidueCount(this);
        }
        return this._props.atomicResidueCount;
    }

    /**
     * Coarse-grained structure, defined as containing less than
     * twice as many elements as polymer residues
     */
    get isCoarseGrained() {
        const ec = this.elementCount;
        const prc = this.polymerResidueCount;
        return prc && ec ? ec / prc < 2 : false;
    }

    get isEmpty() {
        return this.units.length === 0;
    }

    get hashCode() {
        if (this._props.hashCode !== -1) return this._props.hashCode;
        return this.computeHash();
    }

    /** Hash based on all unit.id values in the structure, reflecting the units transformation */
    get transformHash() {
        if (this._props.transformHash !== -1) return this._props.transformHash;
        this._props.transformHash = hashFnv32a(this.units.map(u => u.id));
        return this._props.transformHash;
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
        this._props.hashCode = hash;
        return hash;
    }

    /** Returns a new element location iterator */
    elementLocations(): Iterator<StructureElement.Location> {
        return new Structure.ElementLocationIterator(this);
    }

    /** The parent or itself in case this is the root */
    get root() {
        return this._props.parent || this;
    }

    /** The root/top-most parent or `undefined` in case this is the root */
    get parent() {
        return this._props.parent;
    }

    /**
     * Conformation transformation that was applied to every unit of this structure.
     *
     * Coordinate system applies to the *current* structure only.
     * A parent structure can have a different coordinate system and thefore it has to be composed "manualy"
     * by the consumer.
     */
    get coordinateSystem(): SymmetryOperator {
        return this._props.coordinateSystem;
    }

    get label() {
        return this._props.label;
    }

    get boundary() {
        return this.lookup3d.boundary;
    }

    get lookup3d() {
        if (this._props.lookup3d) return this._props.lookup3d;
        this._props.lookup3d = new StructureLookup3D(this);
        return this._props.lookup3d;
    }

    get interUnitBonds() {
        if (this._props.interUnitBonds) return this._props.interUnitBonds;
        this._props.interUnitBonds = computeInterUnitBonds(this);
        return this._props.interUnitBonds;
    }

    get unitSymmetryGroups(): ReadonlyArray<Unit.SymmetryGroup> {
        if (this._props.unitSymmetryGroups) return this._props.unitSymmetryGroups;
        this._props.unitSymmetryGroups = StructureSymmetry.computeTransformGroups(this);
        return this._props.unitSymmetryGroups;
    }

    /** Maps unit.invariantId to index of SymmetryGroup in unitSymmetryGroups array */
    get unitSymmetryGroupsIndexMap(): IntMap<number> {
        if (this._props.unitSymmetryGroupsIndexMap) return this._props.unitSymmetryGroupsIndexMap;
        this._props.unitSymmetryGroupsIndexMap = Unit.SymmetryGroup.getUnitSymmetryGroupsIndexMap(this.unitSymmetryGroups);
        return this._props.unitSymmetryGroupsIndexMap;
    }

    get carbohydrates(): Carbohydrates {
        if (this._props.carbohydrates) return this._props.carbohydrates;
        this._props.carbohydrates = computeCarbohydrates(this);
        return this._props.carbohydrates;
    }

    get models(): ReadonlyArray<Model> {
        if (this._props.models) return this._props.models;
        this._props.models = getModels(this);
        return this._props.models;
    }

    get uniqueResidueNames() {
        return this._props.uniqueResidueNames
            || (this._props.uniqueResidueNames = getUniqueResidueNames(this));
    }

    get uniqueElementSymbols() {
        return this._props.uniqueElementSymbols
            || (this._props.uniqueElementSymbols = getUniqueElementSymbols(this));
    }

    get entityIndices() {
        return this._props.entityIndices
            || (this._props.entityIndices = getEntityIndices(this));
    }

    get uniqueAtomicResidueIndices() {
        return this._props.uniqueAtomicResidueIndices
            || (this._props.uniqueAtomicResidueIndices = getUniqueAtomicResidueIndices(this));
    }

    /** Contains only atomic units */
    get isAtomic() {
        for (const u of this.units) if (Unit.isAtomic(u)) return false;
        return true;
    }

    /** Contains some atomic units */
    get hasAtomic() {
        for (const u of this.units) if (Unit.isAtomic(u)) return true;
        return false;
    }

    /** Contains only coarse units */
    get isCoarse() {
        for (const u of this.units) if (Unit.isCoarse(u)) return false;
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
        return this._props.serialMapping || (this._props.serialMapping = getSerialMapping(this));
    }

    /**
     * If the structure is based on a single model or has a master-/representative-model, return it.
     * Otherwise throw an exception.
     */
    get model(): Model {
        if (this._props.model) return this._props.model;
        if (this._props.representativeModel) return this._props.representativeModel;
        if (this._props.masterModel) return this._props.masterModel;
        const models = this.models;
        if (models.length > 1) {
            throw new Error('The structure is based on multiple models and has neither a master- nor a representative-model.');
        }
        this._props.model = models[0];
        return this._props.model;
    }

    /** The master-model, other models can have bonds to it  */
    get masterModel(): Model | undefined {
        return this._props.masterModel;
    }

    /** A representative model, e.g. the first model of a trajectory */
    get representativeModel(): Model | undefined {
        return this._props.representativeModel;
    }

    hasElement(e: StructureElement.Location) {
        if (!this.unitMap.has(e.unit.id)) return false;
        return SortedArray.has(this.unitMap.get(e.unit.id).elements, e.element);
    }

    getModelIndex(m: Model) {
        return this.models.indexOf(m);
    }

    private initUnits(units: ArrayLike<Unit>) {
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
        this._props.elementCount = elementCount;
        return { unitMap, unitIndexMap };
    }

    constructor(units: ArrayLike<Unit>, props: Structure.Props = {}) {
        const { unitMap, unitIndexMap } = this.initUnits(units);
        this.unitMap = unitMap;
        this.unitIndexMap = unitIndexMap;
        this.units = units as ReadonlyArray<Unit>;

        if (props.parent) this._props.parent = props.parent.parent || props.parent;

        if (props.coordinateSystem) this._props.coordinateSystem = props.coordinateSystem;
        else if (props.parent) this._props.coordinateSystem = props.parent.coordinateSystem;

        if (props.label) this._props.label = props.label;
        else if (props.parent) this._props.label = props.parent.label;

        if (props.masterModel) this._props.masterModel = props.masterModel;
        else if (props.parent) this._props.masterModel = props.parent.masterModel;

        if (props.representativeModel) this._props.representativeModel = props.representativeModel;
        else if (props.parent) this._props.representativeModel = props.parent.representativeModel;
    }
}

function cmpUnits(units: ArrayLike<Unit>, i: number, j: number) { return units[i].id - units[j].id; }

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

interface SerialMapping {
    /** Cumulative count of preceding elements for each unit */
    cumulativeUnitElementCount: ArrayLike<number>
    /** Unit index for each serial element in the structure */
    unitIndices: ArrayLike<number>
    /** Element index for each serial element in the structure */
    elementIndices: ArrayLike<ElementIndex>
    /** Get serial index of element in the structure */
    getSerialIndex: (unit: Unit, element: ElementIndex) => Structure.SerialIndex
}
function getSerialMapping(structure: Structure): SerialMapping {
    const { units, elementCount, unitIndexMap } = structure;
    const cumulativeUnitElementCount = new Uint32Array(units.length);
    const unitIndices = new Uint32Array(elementCount);
    const elementIndices = new Uint32Array(elementCount) as unknown as ElementIndex[];
    for (let i = 0, m = 0, il = units.length; i < il; ++i) {
        cumulativeUnitElementCount[i] = m;
        const { elements } = units[i];
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            const mj = m + j;
            unitIndices[mj] = i;
            elementIndices[mj] = elements[j];
        }
        m += elements.length;
    }
    return {
        cumulativeUnitElementCount,
        unitIndices,
        elementIndices,

        getSerialIndex: (unit, element) => cumulativeUnitElementCount[unitIndexMap.get(unit.id)] + OrderedSet.indexOf(unit.elements, element) as Structure.SerialIndex
    };
}

namespace Structure {
    export const Empty = new Structure([]);

    export interface Props {
        parent?: Structure
        coordinateSystem?: SymmetryOperator
        label?: string
        /** Master model for structures of a protein model and multiple ligand models */
        masterModel?: Model
        /** Representative model for structures of a model trajectory */
        representativeModel?: Model
    }

    /** Serial index of an element in the structure accross all units */
    export type SerialIndex = { readonly '@type': 'serial-index' } & number

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

    export function create(units: ReadonlyArray<Unit>, props?: Props): Structure {
        return new Structure(units, props);
    }

    export function ofTrajectory(trajectory: ReadonlyArray<Model>): Structure {
        if (trajectory.length === 0) return Empty;

        const units: Unit[] = [];

        let count = 0;
        for (let i = 0, il = trajectory.length; i < il; ++i) {
            const structure = ofModel(trajectory[i]);
            for (let j = 0, jl = structure.units.length; j < jl; ++j) {
                const u = structure.units[j];
                const invariantId = u.invariantId + count;
                const chainGroupId = u.chainGroupId + count;
                const newUnit = Unit.create(units.length, invariantId, chainGroupId, u.traits, u.kind, u.model, u.conformation.operator, u.elements);
                units.push(newUnit);
            }
            count = units.length;
        }

        return create(units, { representativeModel: trajectory[0], label: trajectory[0].label });
    }

    /**
     * Construct a Structure from a model.
     *
     * Generally, a single unit corresponds to a single chain, with the exception
     * of consecutive "single atom chains" with same entity_id and same auth_asym_id.
     */
    export function ofModel(model: Model): Structure {
        const chains = model.atomicHierarchy.chainAtomSegments;
        const { index } = model.atomicHierarchy;
        const { auth_asym_id } = model.atomicHierarchy.chains;
        const { atomicChainOperatorMappinng } = model;
        const builder = new StructureBuilder({ label: model.label });

        for (let c = 0 as ChainIndex; c < chains.count; c++) {
            const operator = atomicChainOperatorMappinng.get(c) || SymmetryOperator.Default;
            const start = chains.offsets[c];

            // set to true for chains that consist of "single atom residues",
            // note that it assumes there are no "zero atom residues"
            let singleAtomResidues = AtomicHierarchy.chainResidueCount(model.atomicHierarchy, c) === chains.offsets[c + 1] - chains.offsets[c];

            // merge all consecutive "single atom chains" with same entity_id and same auth_asym_id
            let multiChain = false;
            while (c + 1 < chains.count
                && chains.offsets[c + 1] - chains.offsets[c] === 1
                && chains.offsets[c + 2] - chains.offsets[c + 1] === 1
            ) {
                singleAtomResidues = true;
                const e1 = index.getEntityFromChain(c);
                const e2 = index.getEntityFromChain(c + 1 as ChainIndex);
                if (e1 !== e2) break;

                const a1 = auth_asym_id.value(c);
                const a2 = auth_asym_id.value(c + 1);
                if (a1 !== a2) break;

                const op1 = atomicChainOperatorMappinng.get(c);
                const op2 = atomicChainOperatorMappinng.get(c + 1 as ChainIndex);
                if (op1 !== op2) break;

                multiChain = true;
                c++;
            }

            const elements = SortedArray.ofBounds(start as ElementIndex, chains.offsets[c + 1] as ElementIndex);

            // check for polymer to exclude CA/P-only models
            if (singleAtomResidues && !isPolymerChain(model, c)) {
                partitionAtomicUnitByAtom(model, elements, builder, multiChain, operator);
            } else if (elements.length > 200000 || isWaterChain(model, c)) {
                // split up very large chains e.g. lipid bilayers, micelles or water with explicit H
                partitionAtomicUnitByResidue(model, elements, builder, multiChain, operator);
            } else {
                builder.addUnit(Unit.Kind.Atomic, model, operator, elements, multiChain ? Unit.Trait.MultiChain : Unit.Trait.None);
            }
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

    function isPolymerChain(model: Model, chainIndex: ChainIndex) {
        const e = model.atomicHierarchy.index.getEntityFromChain(chainIndex);
        return model.entities.data.type.value(e) === 'polymer';
    }

    function partitionAtomicUnitByAtom(model: Model, indices: SortedArray, builder: StructureBuilder, multiChain: boolean, operator: SymmetryOperator) {
        const { x, y, z } = model.atomicConformation;
        const position = { x, y, z, indices };
        const lookup = GridLookup3D(position, getBoundary(position), 8192);
        const { offset, count, array } = lookup.buckets;

        const traits = (multiChain ? Unit.Trait.MultiChain : Unit.Trait.None) | (offset.length > 1 ? Unit.Trait.Partitioned : Unit.Trait.None);

        builder.beginChainGroup();
        for (let i = 0, _i = offset.length; i < _i; i++) {
            const start = offset[i];
            const set = new Int32Array(count[i]);
            for (let j = 0, _j = count[i]; j < _j; j++) {
                set[j] = indices[array[start + j]];
            }
            builder.addUnit(Unit.Kind.Atomic, model, operator, SortedArray.ofSortedArray(set), traits);
        }
        builder.endChainGroup();
    }

    // keeps atoms of residues together
    function partitionAtomicUnitByResidue(model: Model, indices: SortedArray, builder: StructureBuilder, multiChain: boolean, operator: SymmetryOperator) {
        const { residueAtomSegments } = model.atomicHierarchy;

        const startIndices: number[] = [];
        const endIndices: number[] = [];

        const residueIt = Segmentation.transientSegments(residueAtomSegments, indices);
        while (residueIt.hasNext) {
            const residueSegment = residueIt.move();
            startIndices[startIndices.length] = indices[residueSegment.start];
            endIndices[endIndices.length] = indices[residueSegment.end];
        }

        const firstResidueAtomCount = endIndices[0] - startIndices[0];
        const gridCellCount = 512 * firstResidueAtomCount;

        const { x, y, z } = model.atomicConformation;
        const position = { x, y, z, indices: SortedArray.ofSortedArray(startIndices) };
        const lookup = GridLookup3D(position, getBoundary(position), gridCellCount);
        const { offset, count, array } = lookup.buckets;

        const traits = (multiChain ? Unit.Trait.MultiChain : Unit.Trait.None) | (offset.length > 1 ? Unit.Trait.Partitioned : Unit.Trait.None);

        builder.beginChainGroup();
        for (let i = 0, _i = offset.length; i < _i; i++) {
            const start = offset[i];
            const set: number[] = [];
            for (let j = 0, _j = count[i]; j < _j; j++) {
                const k = array[start + j];
                for (let l = startIndices[k], _l = endIndices[k]; l < _l; l++) {
                    set[set.length] = l;
                }
            }
            builder.addUnit(Unit.Kind.Atomic, model, operator, SortedArray.ofSortedArray(new Int32Array(set)), traits);
        }
        builder.endChainGroup();
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
        return new Structure(units, { parent: s, coordinateSystem: newCS });
    }

    export class StructureBuilder {
        private units: Unit[] = [];
        private invariantId = idFactory()

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

        addWithOperator(unit: Unit, operator: SymmetryOperator): Unit {
            return this.add(unit.applyOperator(this.units.length, operator));
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
    export function areUnitAndIndicesEqual(a: Structure, b: Structure) {
        if (a === b) return true;

        if (a.elementCount !== b.elementCount) return false;
        const len = a.units.length;
        if (len !== b.units.length) return false;

        for (let i = 0; i < len; i++) {
            if (a.units[i].id !== b.units[i].id) return false;
        }

        for (let i = 0; i < len; i++) {
            if (!SortedArray.areEqual(a.units[i].elements, b.units[i].elements)) return false;
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

    const distVec = Vec3.zero();
    function unitElementMinDistance(unit: Unit, p: Vec3, eRadius: number) {
        const { elements, conformation: { position, r } } = unit, dV = distVec;
        let minD = Number.MAX_VALUE;
        for (let i = 0, _i = elements.length; i < _i; i++) {
            const e = elements[i];
            const d = Vec3.distance(p, position(e, dV)) - eRadius - r(e);
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

    const distPivot = Vec3.zero();
    export function distance(a: Structure, b: Structure) {
        if (a.elementCount === 0 || b.elementCount === 0) return 0;

        const { units } = a;
        let minD = Number.MAX_VALUE;

        for (let i = 0, _i = units.length; i < _i; i++) {
            const unit = units[i];
            const { elements, conformation: { position, r } } = unit;
            for (let i = 0, _i = elements.length; i < _i; i++) {
                const e = elements[i];
                const d = minDistanceToPoint(b, position(e, distPivot), r(e));
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
     * and unit pairs if within a max distance.
     */
    export function eachUnitPair(structure: Structure, callback: (unitA: Unit, unitB: Unit) => void, props: EachUnitPairProps) {
        const { maxRadius, validUnit, validUnitPair } = props;
        if (!structure.units.some(u => validUnit(u))) return;

        const lookup = structure.lookup3d;
        const imageCenter = Vec3.zero();

        for (const unit of structure.units) {
            if (!validUnit(unit)) continue;

            const bs = unit.lookup3d.boundary.sphere;
            Vec3.transformMat4(imageCenter, bs.center, unit.conformation.operator.matrix);
            const closeUnits = lookup.findUnitIndices(imageCenter[0], imageCenter[1], imageCenter[2], bs.radius + maxRadius);
            for (let i = 0; i < closeUnits.count; i++) {
                const other = structure.units[closeUnits.indices[i]];
                if (!validUnit(other) || unit.id >= other.id || !validUnitPair(unit, other)) continue;

                if (other.elements.length >= unit.elements.length) callback(unit, other);
                else callback(other, unit);
            }
        }
    }

    //

    export const DefaultSizeThresholds = {
        /** Must be lower to be small */
        smallResidueCount: 10,
        /** Must be lower to be medium */
        mediumResidueCount: 3000,
        /** Must be lower to be large (big ribosomes like 4UG0 should still be `large`) */
        largeResidueCount: 20000,
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
}

export default Structure;