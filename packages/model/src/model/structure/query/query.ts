/**
 * Copyright (c) 2017-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Structure } from '../structure.js';
import { StructureSelection } from './selection.js';
import { QueryContext, type QueryFn, type QueryContextOptions } from './context.js';

interface StructureQuery extends QueryFn<StructureSelection> { }
namespace StructureQuery {
    export function run(query: StructureQuery, structure: Structure, options?: QueryContextOptions) {
        return query(new QueryContext(structure, options));
    }

    export function loci(query: StructureQuery, structure: Structure, options?: QueryContextOptions) {
        const sel = query(new QueryContext(structure, options));
        return StructureSelection.toLociWithSourceUnits(sel);
    }
}

export { StructureQuery };