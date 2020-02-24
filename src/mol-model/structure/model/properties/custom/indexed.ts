/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ResidueIndex, ChainIndex, ElementIndex, EntityIndex } from '../../indexing';
import { Unit, Structure, StructureElement } from '../../../structure';
import { Segmentation } from '../../../../../mol-data/int';
import { UUID } from '../../../../../mol-util';
import { CifWriter } from '../../../../../mol-io/writer/cif';
import { Model } from '../../model';

export interface IndexedCustomProperty<Idx extends IndexedCustomProperty.Index, T = any> {
    readonly id: UUID,
    readonly kind: Unit.Kind,
    readonly level: IndexedCustomProperty.Level,
    has(idx: Idx): boolean,
    get(idx: Idx): T | undefined,
    getElements(structure: Structure): IndexedCustomProperty.Elements<T>
}

export namespace IndexedCustomProperty {
    export type Index = ElementIndex | ResidueIndex | ChainIndex | EntityIndex
    export type Level = 'atom' | 'residue' | 'chain' | 'entity'

    export interface Elements<T> {
        elements: StructureElement.Location[],
        property(index: number): T
    }

    export function getCifDataSource<Idx extends Index, T>(structure: Structure, prop: IndexedCustomProperty<Idx, T> | undefined, cache: any): CifWriter.Category.Instance['source'][0] {
        if (!prop) return { rowCount: 0 };
        if (cache && cache[prop.id]) return cache[prop.id];
        const data = prop.getElements(structure);
        const ret = { data, rowCount: data.elements.length };
        if (cache) cache[prop.id] = ret;
        return ret;
    }

    export type Atom<T> = IndexedCustomProperty<ElementIndex, T>
    export function fromAtomMap<T>(map: Map<ElementIndex, T>): Atom<T> {
        return new ElementMappedCustomProperty(map);
    }

    export function fromAtomArray<T>(array: ArrayLike<T>): Atom<T> {
        // TODO: create "array based custom property" as optimization
        return new ElementMappedCustomProperty(arrayToMap(array));
    }

    export type Residue<T> = IndexedCustomProperty<ResidueIndex, T>
    const getResidueSegments = (model: Model) => model.atomicHierarchy.residueAtomSegments;
    export function fromResidueMap<T>(map: Map<ResidueIndex, T>): Residue<T> {
        return new SegmentedMappedIndexedCustomProperty('residue', map, getResidueSegments, Unit.Kind.Atomic);
    }

    export function fromResidueArray<T>(array: ArrayLike<T>): Residue<T> {
        // TODO: create "array based custom property" as optimization
        return new SegmentedMappedIndexedCustomProperty('residue', arrayToMap(array), getResidueSegments, Unit.Kind.Atomic);
    }

    export type Chain<T> = IndexedCustomProperty<ChainIndex, T>
    const getChainSegments = (model: Model) => model.atomicHierarchy.chainAtomSegments;
    export function fromChainMap<T>(map: Map<ChainIndex, T>): Chain<T> {
        return new SegmentedMappedIndexedCustomProperty('chain', map, getChainSegments, Unit.Kind.Atomic);
    }

    export function fromChainArray<T>(array: ArrayLike<T>): Chain<T> {
        // TODO: create "array based custom property" as optimization
        return new SegmentedMappedIndexedCustomProperty('chain', arrayToMap(array), getChainSegments, Unit.Kind.Atomic);
    }

    export type Entity<T> = IndexedCustomProperty<EntityIndex, T>
    export function fromEntityMap<T>(map: Map<EntityIndex, T>): Entity<T> {
        return new EntityMappedCustomProperty(map);
    }
}

function arrayToMap<Idx extends IndexedCustomProperty.Index, T>(array: ArrayLike<T>): Map<Idx, T> {
    const ret = new Map<Idx, T>();
    for (let i = 0 as Idx, _i = array.length; i < _i; i++) ret.set(i, array[i as number]);
    return ret;
}

class SegmentedMappedIndexedCustomProperty<Idx extends IndexedCustomProperty.Index, T = any> implements IndexedCustomProperty<Idx, T> {
    readonly id: UUID = UUID.create22();
    readonly kind: Unit.Kind;
    has(idx: Idx): boolean { return this.map.has(idx); }
    get(idx: Idx) { return this.map.get(idx); }

    private getStructureElements(structure: Structure) {
        const models = structure.models;
        if (models.length !== 1) throw new Error(`Only works on structures with a single model.`);

        const seenIndices = new Set<Idx>();
        const unitGroups = structure.unitSymmetryGroups;
        const loci: StructureElement.Location[] = [];

        const segments = this.segmentGetter(models[0]);

        for (const unitGroup of unitGroups) {
            const unit = unitGroup.units[0];
            if (unit.kind !== this.kind) {
                continue;
            }

            const chains = Segmentation.transientSegments(segments, unit.elements);
            while (chains.hasNext) {
                const seg = chains.move();
                if (!this.has(seg.index) || seenIndices.has(seg.index)) continue;
                seenIndices.add(seg.index);
                loci[loci.length] = StructureElement.Location.create(structure, unit, unit.elements[seg.start]);
            }
        }

        loci.sort((x, y) => x.element - y.element);
        return loci;
    }

    getElements(structure: Structure): IndexedCustomProperty.Elements<T> {
        const index = this.segmentGetter(structure.model).index;
        const elements = this.getStructureElements(structure);
        return { elements, property: i => this.get(index[elements[i].element])! };
    }

    constructor(public level: 'residue' | 'chain', private map: Map<Idx, T>, private segmentGetter: (model: Model) => Segmentation<ElementIndex, Idx>, kind: Unit.Kind) {
        this.kind = kind;
    }
}

class ElementMappedCustomProperty<T = any> implements IndexedCustomProperty<ElementIndex, T> {
    readonly id: UUID = UUID.create22();
    readonly kind: Unit.Kind;
    readonly level = 'atom';
    has(idx: ElementIndex): boolean { return this.map.has(idx); }
    get(idx: ElementIndex) { return this.map.get(idx); }

    private getStructureElements(structure: Structure) {
        const models = structure.models;
        if (models.length !== 1) throw new Error(`Only works on structures with a single model.`);

        const seenIndices = new Set<ElementIndex>();
        const unitGroups = structure.unitSymmetryGroups;
        const loci: StructureElement.Location[] = [];

        for (const unitGroup of unitGroups) {
            const unit = unitGroup.units[0];
            if (unit.kind !== this.kind) {
                continue;
            }

            const elements = unit.elements;
            for (let i = 0, _i = elements.length; i < _i; i++) {
                const e = elements[i];
                if (!this.has(e) || seenIndices.has(e)) continue;
                seenIndices.add(elements[i]);
                loci[loci.length] = StructureElement.Location.create(structure, unit, e);
            }
        }

        loci.sort((x, y) => x.element - y.element);
        return loci;
    }

    getElements(structure: Structure): IndexedCustomProperty.Elements<T> {
        const elements = this.getStructureElements(structure);
        return { elements, property: i => this.get(elements[i].element)! };
    }

    constructor(private map: Map<ElementIndex, T>) {
        this.kind = Unit.Kind.Atomic;
    }
}

class EntityMappedCustomProperty<T = any> implements IndexedCustomProperty<EntityIndex, T> {
    readonly id: UUID = UUID.create22();
    readonly kind: Unit.Kind;
    readonly level = 'entity';
    has(idx: EntityIndex): boolean { return this.map.has(idx); }
    get(idx: EntityIndex) { return this.map.get(idx); }

    private getStructureElements(structure: Structure) {
        const models = structure.models;
        if (models.length !== 1) throw new Error(`Only works on structures with a single model.`);

        const index = models[0].atomicHierarchy.index;
        const seenIndices = new Set<EntityIndex>();
        const unitGroups = structure.unitSymmetryGroups;
        const loci: StructureElement.Location[] = [];

        const segments = models[0].atomicHierarchy.chainAtomSegments;

        for (const unitGroup of unitGroups) {
            const unit = unitGroup.units[0];
            if (unit.kind !== this.kind) {
                continue;
            }

            const chains = Segmentation.transientSegments(segments, unit.elements);
            while (chains.hasNext) {
                const seg = chains.move();
                const eI = index.getEntityFromChain(seg.index);
                if (!this.has(eI) || seenIndices.has(eI)) continue;
                seenIndices.add(eI);
                loci[loci.length] = StructureElement.Location.create(structure, unit, unit.elements[seg.start]);
            }
        }

        loci.sort((x, y) => x.element - y.element);
        return loci;
    }

    getElements(structure: Structure): IndexedCustomProperty.Elements<T> {
        const elements = this.getStructureElements(structure);
        const chainIndex = structure.model.atomicHierarchy.chainAtomSegments.index;
        const index = structure.model.atomicHierarchy.index;
        return { elements, property: i => this.get(index.getEntityFromChain(chainIndex[elements[i].element]))! };
    }

    constructor(private map: Map<EntityIndex, T>) {
        this.kind = Unit.Kind.Atomic;
    }
}