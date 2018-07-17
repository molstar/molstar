/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray } from 'mol-data/int';
import { Structure } from '../../structure';
import { QueryPredicate } from '../context';
import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';

export function pick(query: StructureQuery, pred: QueryPredicate): StructureQuery {
    return ctx => {
        const sel = query(ctx);

        if (StructureSelection.isSingleton(sel)) {
            const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
            for (const unit of ctx.inputStructure.units) {
                const { elements } = unit;
                for (let i = 0, _i = elements.length; i < _i; i++) {
                    // TODO: optimize this somehow???
                    const s = Structure.create([unit.getChild(SortedArray.ofSingleton(elements[i]))]);
                    ctx.lockCurrentStructure(s);
                    if (pred(ctx)) ret.add(s);
                    ctx.unlockCurrentStructure();
                }
            }
            return ret.getSelection();
        } else {
            const ret = StructureSelection.LinearBuilder(ctx.inputStructure);
            for (const s of sel.structures) {
                ctx.lockCurrentStructure(s);
                if (pred(ctx)) ret.add(s);
                ctx.unlockCurrentStructure();
            }
            return ret.getSelection();
        }
    };
}
