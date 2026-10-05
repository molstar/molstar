/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Koya Sakuma
 * Adapted from MolQL implemtation of atom-set.ts
 *
 * Copyright (c) 2017 MolQL contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StructureQuery } from '../query.js';
import { StructureSelection } from '../selection.js';
import { getCurrentStructureProperties } from './filters.js';
import type { QueryContext, QueryFn } from '../context.js';


export function atomCount(ctx: QueryContext) {
    return ctx.currentStructure.elementCount;
}


export function countQuery(query: StructureQuery) {
    return (ctx: QueryContext) => {
        const sel = query(ctx);
        return StructureSelection.structureCount(sel);
    };
}

export function propertySet(prop: QueryFn<any>) {
    return (ctx: QueryContext) => {
        const set = new Set();
        return getCurrentStructureProperties(ctx, prop, set);
    };
}

