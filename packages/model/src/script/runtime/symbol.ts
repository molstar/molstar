// /**
//  * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// import Environment from './environment'
// import RuntimeExpression from './expression'
// import Expression from '../language/expression';

// type SymbolRuntime = SymbolRuntime.Dynamic | SymbolRuntime.Static

// namespace SymbolRuntime {
//     export interface Static {
//         kind: 'static',
//         readonly runtime: (ctx: any, args: Arguments) => any,
//         readonly attributes: Attributes
//     }

//     export interface Dynamic {
//         kind: 'dynamic',
//         readonly compile: (env: Environment, args: Expression.Arguments) => RuntimeExpression
//     }

//     export interface Attributes { isStatic: boolean }

//     export type Table = Map<string, SymbolRuntime>

//     export type Arguments = ArrayLike<RuntimeExpression> | { [name: string]: RuntimeExpression | undefined }
// }

// export { SymbolRuntime }