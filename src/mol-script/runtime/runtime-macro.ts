// /**
//  * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// import Expression from '../language/expression';

// interface Macro {
//     readonly argNames: ReadonlyArray<string>,
//     readonly argIndex: { [name: string]: number },
//     readonly expression: Expression
// }

// namespace Macro {
//     export type Table = Map<string, Macro>

//     function subst(table: Table, expr: Expression, argIndex: { [name: string]: number }, args: ArrayLike<Expression>): Expression {
//         if (Expression.isLiteral(expr)) return expr;
//         if (Expression.isSymbol(expr)) {
//             const idx = argIndex[expr.name];
//             if (typeof idx !== 'undefined') return args[idx];

//             if (table.has(expr.name)) {
//                 const macro = table.get(expr.name)!;
//                 if (macro.argNames.length === 0) return macro.expression;
//             }

//             return expr;
//         }

//         const head = subst(table, expr.head, argIndex, args);
//         const headChanged = head !== expr.head;
//         if (!expr.args) {
//             return headChanged ? Expression.Apply(head) : expr;
//         }

//         let argsChanged = false;

//         if (Expression.isArgumentsArray(expr.args)) {
//             let newArgs: Expression[] = [];
//             for (let i = 0, _i = expr.args.length; i < _i; i++) {
//                 const oldArg = expr.args[i];
//                 const newArg = subst(table, oldArg, argIndex, args);
//                 if (oldArg !== newArg) argsChanged = true;
//                 newArgs[newArgs.length] = newArg;
//             }
//             if (!argsChanged) newArgs = expr.args;

//             if (Expression.isSymbol(head) && table.has(head.name)) {
//                 const macro = table.get(head.name)!;
//                 if (macro.argNames.length === newArgs.length) {
//                     return subst(table, macro.expression, macro.argIndex, newArgs);
//                 }
//             }

//             if (!headChanged && !argsChanged) return expr;
//             return Expression.Apply(head, newArgs);
//         } else {

//             let newArgs: any = {}
//             for (const key of Object.keys(expr.args)) {
//                 const oldArg = expr.args[key];
//                 const newArg = subst(table, oldArg, argIndex, args);
//                 if (oldArg !== newArg) argsChanged = true;
//                 newArgs[key] = newArg;
//             }
//             if (!headChanged && !argsChanged) return expr;
//             if (!argsChanged) newArgs = expr.args;

//             return Expression.Apply(head, newArgs);
//         }
//     }

//     export function substitute(table: Table, macro: Macro, args: ArrayLike<Expression>) {
//         if (args.length === 0) return macro.expression;
//         return subst(table, macro.expression, macro.argIndex, args);
//     }
// }

// export { Macro }