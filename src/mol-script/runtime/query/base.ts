/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Expression from '../../language/expression';
import { QueryContext, QueryFn, Structure } from '../../../mol-model/structure';
import { MSymbol } from '../../language/symbol';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';

export class QueryRuntimeTable {
    private map = new Map<string, QuerySymbolRuntime>();

    removeSymbol(runtime: QuerySymbolRuntime) {
        this.map.delete(runtime.symbol.id);
    }

    addSymbol(runtime: QuerySymbolRuntime) {
        if (this.map.has(runtime.symbol.id)) {
            console.warn(`Symbol '${runtime.symbol.id}' already added. Call removeSymbol/removeCustomProps re-adding the symbol.`);
        }
        this.map.set(runtime.symbol.id, runtime);
    }

    addCustomProp(desc: CustomPropertyDescriptor<any>) {
        if (!desc.symbols) return;

        for (const k of Object.keys(desc.symbols)) {
            this.addSymbol((desc.symbols as any)[k]);
        }
    }

    removeCustomProp(desc: CustomPropertyDescriptor<any>) {
        if (!desc.symbols) return;

        for (const k of Object.keys(desc.symbols)) {
            this.removeSymbol((desc.symbols as any)[k]);
        }
    }

    getRuntime(id: string) {
        return this.map.get(id);
    }
}

export const DefaultQueryRuntimeTable = new QueryRuntimeTable();

export class QueryCompilerCtx {
    constQueryContext: QueryContext = new QueryContext(Structure.Empty);

    constructor(public table: QueryRuntimeTable) {

    }
}

export type ConstQuerySymbolFn<S extends MSymbol = MSymbol> = (ctx: QueryContext, args: QueryRuntimeArguments<S>) => any
export type QuerySymbolFn<S extends MSymbol = MSymbol> = (ctx: QueryContext, args: QueryRuntimeArguments<S>) => any


export type QueryCompiledSymbolRuntime = { kind: 'const', value: any } | { kind: 'dynamic', runtime: QuerySymbolFn }

export type CompiledQueryFn<T = any> = { isConst: boolean, fn: QueryFn }

export namespace QueryCompiledSymbol {
    export function Const(value: any): QueryCompiledSymbolRuntime  {
        return { kind: 'const', value };
    }

    export function Dynamic(runtime: QuerySymbolFn): QueryCompiledSymbolRuntime {
        return { kind: 'dynamic', runtime };
    }
}

export namespace CompiledQueryFn {
    export function Const(value: any): CompiledQueryFn  {
        return { isConst: true, fn: function CompiledQueryFn_Const(ctx) { return value; } };
    }

    export function Dynamic(fn: QueryFn): CompiledQueryFn {
        return { isConst: false, fn };
    }
}

export interface QuerySymbolRuntime {
    symbol: MSymbol,
    compile(ctx: QueryCompilerCtx, args?: Expression.Arguments): CompiledQueryFn
}

export type QueryRuntimeArguments<S extends MSymbol> =
    { length?: number } & { [P in keyof S['args']['@type']]: QueryFn<S['args']['@type'][P]> }

export namespace QueryRuntimeArguments {
    export function forEachEval<S extends MSymbol, Ctx>(xs: QueryRuntimeArguments<S>, queryCtx: QueryContext, f: (arg: any, i: number, ctx: Ctx) => void, ctx: Ctx): Ctx {
        if (typeof xs.length === 'number') {
            for (let i = 0, _i = xs.length; i < _i; i++) f((xs as any)[i](queryCtx), i, ctx);
        } else {
            let i = 0;
            for (const k of Object.keys(xs)) f((xs as any)[k](queryCtx), i++, ctx);
        }
        return ctx;
    }
}

export namespace QuerySymbolRuntime {
    export function Const<S extends MSymbol<any>>(symbol: S, fn: ConstQuerySymbolFn<S>): QuerySymbolRuntime {
        return new SymbolRuntimeImpl(symbol, fn, true);
    }

    export function Dynamic<S extends MSymbol<any>>(symbol: S, fn: QuerySymbolFn<S>): QuerySymbolRuntime {
        return new SymbolRuntimeImpl(symbol, fn, false);
    }
}

class SymbolRuntimeImpl<S extends MSymbol> implements QuerySymbolRuntime {
    compile(ctx: QueryCompilerCtx, inputArgs?: Expression.Arguments): CompiledQueryFn {
        let args: any, constArgs = false;
        if (!inputArgs) {
            args = void 0;
            constArgs = true;
        } else if (Expression.isArgumentsArray(inputArgs)) {
            args = [];
            constArgs = false;
            for (const arg of inputArgs) {
                const compiled = _compile(ctx, arg);
                constArgs = constArgs && compiled.isConst;
                args.push(compiled.fn);
            }
        } else {
            args = Object.create(null);
            constArgs = false;
            for (const key of Object.keys(inputArgs)) {
                const compiled = _compile(ctx, inputArgs[key]);
                constArgs = constArgs && compiled.isConst;
                args[key] = compiled.fn;
            }
        }

        if (this.isConst) {
            if (this.isConst && constArgs) {
                return CompiledQueryFn.Const(this.fn(ctx.constQueryContext, args));
            }

            return CompiledQueryFn.Dynamic(createDynamicFn(this.fn, args));
        }

        return CompiledQueryFn.Dynamic(createDynamicFn(this.fn, args));
    }

    constructor(public symbol: S, private fn: QuerySymbolFn<S>, private isConst: boolean) {

    }
}

function createDynamicFn<S extends MSymbol>(fn: QuerySymbolFn<S>, args: any): QueryFn {
    return function DynamicFn(ctx) { return fn(ctx, args); };
}

function _compile(ctx: QueryCompilerCtx, expression: Expression): CompiledQueryFn {
    if (Expression.isLiteral(expression)) {
        return CompiledQueryFn.Const(expression);
    }

    if (Expression.isSymbol(expression)) {
        const runtime = ctx.table.getRuntime(expression.name);
        if (!runtime) return CompiledQueryFn.Const(expression.name);

        return runtime.compile(ctx);
    }

    if (!Expression.isSymbol(expression.head)) {
        throw new Error('Can only apply symbols.');
    }

    const compiler = ctx.table.getRuntime(expression.head.name);
    if (!compiler) {
        throw new Error(`Symbol '${expression.head.name}' is not implemented.`);
    }

    return compiler.compile(ctx, expression.args);
}

export function compile<T = any>(expression: Expression): QueryFn<T> {
    const ctx = new QueryCompilerCtx(DefaultQueryRuntimeTable);
    return _compile(ctx, expression).fn;
}