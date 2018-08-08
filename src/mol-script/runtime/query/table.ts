/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MolScriptSymbolTable as MolScript } from '../../language/symbol-table';
import { DefaultQueryRuntimeTable, QuerySymbolRuntime } from './compiler';

import C = QuerySymbolRuntime.Const
import D = QuerySymbolRuntime.Dynamic
import { Queries, StructureProperties } from 'mol-model/structure';
import { ElementSymbol } from 'mol-model/structure/model/types';

const symbols = [
    C(MolScript.core.math.add, (ctx, xs) => {
        let ret = 0;
        if (typeof xs.length === 'number') {
            for (let i = 0, _i = xs.length; i < _i; i++) ret += xs[i](ctx);
        } else {
            for (const k of Object.keys(xs)) ret += xs[k](ctx);
        }
        return ret;
    }),

    // ============= RELATIONAL ================
    C(MolScript.core.rel.eq, (ctx, v) => v[0](ctx) === v[1](ctx)),
    C(MolScript.core.rel.neq, (ctx, v) => v[0](ctx) !== v[1](ctx)),
    C(MolScript.core.rel.lt, (ctx, v) => v[0](ctx) < v[1](ctx)),
    C(MolScript.core.rel.lte, (ctx, v) => v[0](ctx) <= v[1](ctx)),
    C(MolScript.core.rel.gr, (ctx, v) => v[0](ctx) > v[1](ctx)),
    C(MolScript.core.rel.gre, (ctx, v) => v[0](ctx) >= v[1](ctx)),

    // ============= TYPES ================
    C(MolScript.structureQuery.type.elementSymbol, (ctx, v) => ElementSymbol(v[0](ctx))),

    // ============= GENERATORS ================
    D(MolScript.structureQuery.generator.atomGroups, (ctx, xs) => Queries.generators.atoms({
        entityTest: xs['entity-test'],
        chainTest: xs['chain-test'],
        residueTest: xs['residue-test'],
        atomTest: xs['atom-test'],
        groupBy: xs['group-by']
    })(ctx)),

    // ============= ATOM PROPERTIES ================
    D(MolScript.structureQuery.atomProperty.macromolecular.label_comp_id, (ctx, _) => StructureProperties.residue.label_comp_id(ctx.element)),
    D(MolScript.structureQuery.atomProperty.core.elementSymbol, (ctx, _) => StructureProperties.atom.type_symbol(ctx.element))

    // Symbol(MolQL.core.rel.neq, staticAttr)((env, v) => v[0](env) !== v[1](env)),
    // Symbol(MolQL.core.rel.lt, staticAttr)((env, v) => v[0](env) < v[1](env)),
    // Symbol(MolQL.core.rel.lte, staticAttr)((env, v) => v[0](env) <= v[1](env)),
    // Symbol(MolQL.core.rel.gr, staticAttr)((env, v) => v[0](env) > v[1](env)),
    // Symbol(MolQL.core.rel.gre, staticAttr)((env, v) => v[0](env) >= v[1](env)),
    // Symbol(MolQL.core.rel.inRange, staticAttr)((env, v) => {
    //     const x = v[0](env);
    //     return x >= v[1](env) && x <= v[2](env)
    // }),
];

(function () {
    for (const s of symbols) {
        DefaultQueryRuntimeTable.addSymbol(s);
    }
})();