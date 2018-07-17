/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// import { StructureQuery } from '../query'
// import { StructureSelection } from '../selection'
// import { Unit, StructureProperties as P, Structure } from '../../structure'
// import { Segmentation, SortedArray } from 'mol-data/int'
// import { LinearGroupingBuilder } from '../utils/builders';
// import { QueryPredicate, QueryFn, QueryContextView } from '../context';

// export function pick(query: StructureQuery, pred: QueryPredicate): StructureQuery {
//     return async ctx => {
//         const sel = await query(ctx);

//         if (StructureSelection.isSingleton(sel)) {
//             const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
//             for (const unit of ctx.inputStructure.units) {
//                 const { elements } = unit;
//                 for (let i = 0, _i = elements.length; i < _i; i++) {
//                     // TODO: optimize this somehow???
//                     const s = Structure.create([unit.getChild(SortedArray.ofSingleton(elements[i]))]);
//                     ctx.lockCurrentStructure(s);
//                     if (pred(ctx)) ret.add(s);
//                     ctx.unlockCurrentStructure();
//                 }
//             }
//             return ret.getSelection();
//         } else {
//             const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
//             for (const s of sel.structures) {
//                 ctx.lockCurrentStructure(s);
//                 if (pred(ctx)) ret.add(s);
//                 ctx.unlockCurrentStructure();
//             }
//             return ret.getSelection();
//         }
//     };
// }
