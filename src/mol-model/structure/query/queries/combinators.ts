/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureQuery } from '../query';
import { StructureSelection } from '../selection';

export function merge(queries: ArrayLike<StructureQuery>): StructureQuery {
    return ctx => {
        const ret = StructureSelection.UniqueBuilder(ctx.inputStructure);
        for (let i = 0; i < queries.length; i++) {
            StructureSelection.forEach(queries[i](ctx), s => {
                ret.add(s);
            });
        }
        return ret.getSelection();
    }
}

// TODO: intersect, distanceCluster