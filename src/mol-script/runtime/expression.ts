// /**
//  * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// import Environment from './environment'

// type RuntimeExpression<T = any> = (env: Environment) => T

// export interface ExpressionInfo {
//     isConst?: boolean
// }

// namespace RuntimeExpression {
//     export function constant<T>(c: T): RuntimeExpression<T> {
//         return env => c;
//     }

//     export function func<T>(f: (env: Environment) => T): RuntimeExpression<T> {
//         return env => f(env);
//     }
// }

// export default RuntimeExpression