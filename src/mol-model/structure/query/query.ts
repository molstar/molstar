/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from '../structure'
import { StructureSelection } from './selection'
import { QueryContext, QueryFn, QueryContextOptions } from './context';
import Expression from '../../../mol-script/language/expression';
import { compile } from '../../../mol-script/runtime/query/compiler';
import { MolScriptBuilder } from '../../../mol-script/language/builder';

interface StructureQuery extends QueryFn<StructureSelection> { }
namespace StructureQuery {
    export function run(query: StructureQuery, structure: Structure, options?: QueryContextOptions) {
        return query(new QueryContext(structure, options));
    }

    export function runExpr(expr: Expression | ((builder: typeof MolScriptBuilder) => Expression), structure: Structure, options?: QueryContextOptions) {
        const e = typeof expr === 'function' ? expr(MolScriptBuilder) : expr;
        const query = compile<StructureSelection>(e);
        return query(new QueryContext(structure, options));
    }
}

export { StructureQuery }