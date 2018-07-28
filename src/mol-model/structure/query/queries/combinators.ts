/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';
import { none } from './generators';

export function merge(queries: ArrayLike<StructureQuery>): StructureQuery {
    if (queries.length === 0) {
        return none;
    } else if (queries.length === 1) {
        return queries[0];
    }
    return ctx => {
        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        for (let i = 0; i < queries.length; i++) {
            StructureSelection.forEach(queries[i](ctx), (s, j) => {
                ret.add(s);
                if (i % 100) ctx.throwIfTimedOut();
            });
        }
        return ret.getSelection();
    }
}

// TODO: intersect, distanceCluster