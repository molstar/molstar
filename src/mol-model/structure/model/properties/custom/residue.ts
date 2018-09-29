/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ResidueIndex } from '../../indexing';
import { Unit, Structure, StructureElement } from '../../../structure';
import { Segmentation } from 'mol-data/int';
import { UUID } from 'mol-util';
import { CifWriter } from 'mol-io/writer/cif';

export interface ResidueCustomProperty<T = any> {
    readonly id: UUID,
    readonly kind: Unit.Kind,
    has(idx: ResidueIndex): boolean
    get(idx: ResidueIndex): T | undefined
}

export namespace ResidueCustomProperty {
    export interface ExportCtx<T> {
        elements: StructureElement[],
        property(index: number): T
    };

    function getExportCtx<T>(prop: ResidueCustomProperty<T>, structure: Structure): ExportCtx<T> {
        const residueIndex = structure.model.atomicHierarchy.residueAtomSegments.index;
        const elements = getStructureElements(structure, prop);
        return { elements, property: i => prop.get(residueIndex[elements[i].element])! };
    }

    export function getCifDataSource<T>(structure: Structure, prop: ResidueCustomProperty<T> | undefined, cache: any): CifWriter.Category.Instance['source'][0] {
        if (!prop) return { rowCount: 0 };
        if (cache && cache[prop.id]) return cache[prop.id];
        const data = getExportCtx(prop, structure);
        const ret = { data, rowCount: data.elements.length };
        if (cache) cache[prop.id] = ret;
        return ret;
    }

    class FromMap<T> implements ResidueCustomProperty<T> {
        readonly id = UUID.create();

        has(idx: ResidueIndex): boolean {
            return this.map.has(idx);
        }

        get(idx: ResidueIndex) {
            return this.map.get(idx);
        }

        constructor(private map: Map<ResidueIndex, T>, public kind: Unit.Kind) {
        }
    }

    export function fromMap<T>(map: Map<ResidueIndex, T>, kind: Unit.Kind) {
        return new FromMap(map, kind);
    }

    /**
     * Gets all StructureElements that correspond to 1st atoms of residues that have an property assigned.
     * Only works correctly for structures with a single model.
     */
    export function getStructureElements(structure: Structure, property: ResidueCustomProperty) {
        const models = structure.models;
        if (models.length !== 1) throw new Error(`Only works on structures with a single model.`);

        const seenResidues = new Set<ResidueIndex>();
        const unitGroups = structure.unitSymmetryGroups;
        const loci: StructureElement[] = [];

        for (const unitGroup of unitGroups) {
            const unit = unitGroup.units[0];
            if (unit.kind !== property.kind) {
                continue;
            }

            const residues = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
            while (residues.hasNext) {
                const seg = residues.move();
                if (!property.has(seg.index) || seenResidues.has(seg.index)) continue;

                seenResidues.add(seg.index);
                loci[loci.length] = StructureElement.create(unit, unit.elements[seg.start]);
            }
        }

        loci.sort((x, y) => x.element - y.element);
        return loci;
    }
}