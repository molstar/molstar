/*
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

type Expression =
    | Expression.Literal
    | Expression.Symbol
    | Expression.Apply

namespace Expression {
    export type Literal = string | number | boolean
    export type Symbol = { kind: 'symbol', name: string }
    export type Arguments = Expression[] | { [name: string]: Expression }
    export interface Apply { readonly head: Expression, readonly args?: Arguments }

    export function Symbol(name: string): Symbol { return { kind: 'symbol', name }; }
    export function Apply(head: Expression, args?: Arguments): Apply { return args ? { head, args } : { head }; }

    export function isArgumentsArray(e: Arguments): e is Expression[] { return e instanceof Array; }
    export function isArgumentsMap(e: Arguments): e is { [name: string]: Expression } { return !(e instanceof Array); }
    export function isLiteral(e: Expression): e is Expression.Literal { return !isApply(e); }
    export function isApply(e: Expression): e is Expression.Apply { return !!e && !!(e as Expression.Apply).head && typeof e === 'object'; }
    export function isSymbol(e: Expression): e is Expression.Symbol { return !!e && (e as any).kind === 'symbol' }
}

export default Expression