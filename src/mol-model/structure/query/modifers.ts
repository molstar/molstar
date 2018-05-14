// /**
//  * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// import Query from './query'
// import Selection from './selection'
// import P from './properties'
// import { Element, Unit } from '../structure'
// import { OrderedSet, Segmentation } from 'mol-data/int'
// import { LinearGroupingBuilder } from './utils/builders';

// export function wholeResidues(query: Query, isFlat: boolean): Query.Provider {
//     return async (structure, ctx) => {
//         const selection = query(structure).runAsChild(ctx);
//         const { units } = structure;
//         const l = Element.Location();
//         const builder = structure.subsetBuilder(true);

//         for (const unit of units) {
//             l.unit = unit;
//             const elements = unit.elements;

//             builder.beginUnit(unit.id);
//             for (let j = 0, _j = elements.length; j < _j; j++) {
//                 l.element = elements[j];
//                 if (atomTest(l)) builder.addElement(l.element);
//             }
//             builder.commitUnit();

//             if (ctx.shouldUpdate) await ctx.update({ message: 'Atom Groups', current: 0, max: units.length });
//         }

//         return Selection.Singletons(structure, builder.getStructure());
//     };
// }

// export interface IncludeSurroundingsParams {
//     selection: Selection,
//     radius: number,
//     atomRadius?: number,
//     wholeResidues?: boolean
// }