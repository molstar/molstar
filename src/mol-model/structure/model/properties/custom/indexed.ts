/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ResidueIndex, ChainIndex, ElementIndex } from '../../indexing';
import { Unit, Structure, StructureElement } from '../../../structure';
import { Segmentation } from 'mol-data/int';
import { UUID } from 'mol-util';
import { CifWriter } from 'mol-io/writer/cif';
import { Model } from '../../model';

export class IndexedCustomProperty<Idx extends IndexedCustomProperty.Index, T = any> {
    readonly id: UUID = UUID.create();
    readonly kind: Unit.Kind;
    has(idx: Idx): boolean { return this.map.has(idx); }
    get(idx: Idx) { return this.map.get(idx); }

    private getStructureElements(structure: Structure) {
        const models = structure.models;
        if (models.length !== 1) throw new Error(`Only works on structures with a single model.`);

        const seenIndices = new Set<Idx>();
        const unitGroups = structure.unitSymmetryGroups;
        const loci: StructureElement[] = [];

        const segments = this.segmentGetter(models[0])

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
                loci[loci.length] = StructureElement.create(unit, unit.elements[seg.start]);
            }
        }

        loci.sort((x, y) => x.element - y.element);
        return loci;
    }

    getExportContext(structure: Structure): IndexedCustomProperty.ExportCtx<T> {
        const index = this.segmentGetter(structure.model).index;
        const elements = this.getStructureElements(structure);
        return { elements, property: i => this.get(index[elements[i].element])! };
    }

    constructor(private map: Map<Idx, T>, private segmentGetter: (model: Model) => Segmentation<ElementIndex, Idx>, kind: Unit.Kind) {
        this.kind = kind;
    }
}

export namespace IndexedCustomProperty {
    export type Index = ResidueIndex | ChainIndex

    export interface ExportCtx<T> {
        elements: StructureElement[],
        property(index: number): T
    }

    export function getCifDataSource<Idx extends Index, T>(structure: Structure, prop: IndexedCustomProperty<Idx, T> | undefined, cache: any): CifWriter.Category.Instance['source'][0] {
        if (!prop) return { rowCount: 0 };
        if (cache && cache[prop.id]) return cache[prop.id];
        const data = prop.getExportContext(structure);
        const ret = { data, rowCount: data.elements.length };
        if (cache) cache[prop.id] = ret;
        return ret;
    }

    export type Residue<T> = IndexedCustomProperty<ResidueIndex, T>
    const getResidueSegments = (model: Model) => model.atomicHierarchy.residueAtomSegments;
    export function fromResidueMap<T>(map: Map<ResidueIndex, T>, kind: Unit.Kind): Residue<T> {
        return new IndexedCustomProperty(map, getResidueSegments, kind);
    }

    export type Chain<T> = IndexedCustomProperty<ChainIndex, T>
    const getChainSegments = (model: Model) => model.atomicHierarchy.chainAtomSegments;
    export function fromChainMap<T>(map: Map<ChainIndex, T>, kind: Unit.Kind): Chain<T> {
        return new IndexedCustomProperty(map, getChainSegments, kind);
    }
}