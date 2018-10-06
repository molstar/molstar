// /**
//  * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// import { ChainIndex } from '../../indexing';
// import { Unit, Structure, StructureElement } from '../../../structure';
// import { Segmentation } from 'mol-data/int';
// import { UUID } from 'mol-util';
// import { CifWriter } from 'mol-io/writer/cif';

// export interface ChainCustomProperty<T = any> {
//     readonly id: UUID,
//     readonly kind: Unit.Kind,
//     has(idx: ChainIndex): boolean
//     get(idx: ChainIndex): T | undefined
// }

// export namespace ChainCustomProperty {
//     export interface ExportCtx<T> {
//         elements: StructureElement[],
//         property(index: number): T
//     };

//     function getExportCtx<T>(prop: ChainCustomProperty<T>, structure: Structure): ExportCtx<T> {
//         const chainIndex = structure.model.atomicHierarchy.chainAtomSegments.index;
//         const elements = getStructureElements(structure, prop);
//         return { elements, property: i => prop.get(chainIndex[elements[i].element])! };
//     }

//     export function getCifDataSource<T>(structure: Structure, prop: ChainCustomProperty<T> | undefined, cache: any): CifWriter.Category.Instance['source'][0] {
//         if (!prop) return { rowCount: 0 };
//         if (cache && cache[prop.id]) return cache[prop.id];
//         const data = getExportCtx(prop, structure);
//         const ret = { data, rowCount: data.elements.length };
//         if (cache) cache[prop.id] = ret;
//         return ret;
//     }

//     class FromMap<T> implements ChainCustomProperty<T> {
//         readonly id = UUID.create();

//         has(idx: ChainIndex): boolean {
//             return this.map.has(idx);
//         }

//         get(idx: ChainIndex) {
//             return this.map.get(idx);
//         }

//         constructor(private map: Map<ChainIndex, T>, public kind: Unit.Kind) {
//         }
//     }

//     export function fromMap<T>(map: Map<ChainIndex, T>, kind: Unit.Kind) {
//         return new FromMap(map, kind);
//     }

//     /**
//      * Gets all StructureElements that correspond to 1st atoms of residues that have an property assigned.
//      * Only works correctly for structures with a single model.
//      */
//     export function getStructureElements(structure: Structure, property: ChainCustomProperty) {
//         const models = structure.models;
//         if (models.length !== 1) throw new Error(`Only works on structures with a single model.`);

//         const seenChains = new Set<ChainIndex>();
//         const unitGroups = structure.unitSymmetryGroups;
//         const loci: StructureElement[] = [];

//         for (const unitGroup of unitGroups) {
//             const unit = unitGroup.units[0];
//             if (unit.kind !== property.kind) {
//                 continue;
//             }

//             const chains = Segmentation.transientSegments(unit.model.atomicHierarchy.chainAtomSegments, unit.elements);
//             while (chains.hasNext) {
//                 const seg = chains.move();
//                 if (!property.has(seg.index) || seenChains.has(seg.index)) continue;

//                 seenChains.add(seg.index);
//                 loci[loci.length] = StructureElement.create(unit, unit.elements[seg.start]);
//             }
//         }

//         loci.sort((x, y) => x.element - y.element);
//         return loci;
//     }
// }