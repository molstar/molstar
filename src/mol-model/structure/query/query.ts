/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RuntimeContext } from 'mol-task'
import { Structure } from '../structure'
import { StructureSelection } from './selection'
import { QueryContext } from './context';

interface StructureQuery { (ctx: QueryContext): Promise<StructureSelection> }
namespace StructureQuery {
    export function run(query: StructureQuery, structure: Structure, ctx?: RuntimeContext) {
        return query(new QueryContext(structure, ctx || RuntimeContext.Synchronous))
    }
}

export { StructureQuery }