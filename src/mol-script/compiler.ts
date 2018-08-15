// /**
//  * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// import Expression from './expression'
// import Environment from './runtime/environment'
// import RuntimeExpression from './runtime/expression'
// import SymbolRuntime, { RuntimeArguments } from './runtime/symbol'

// export type CompiledExpression<T> = () => T

// type Compiler = <T>(expression: Expression) => CompiledExpression<T>
// function Compiler<C>(env: Environment): Compiler {
//     return expression => compile(env, expression).runtime;
// }

// type CompileResult = { isConst: boolean, runtime: RuntimeExpression }

// namespace CompileResult {
//     export function Const(value: any): CompileResult { return { isConst: true, runtime: RuntimeExpression.constant(value) } }
//     export function Dynamic(runtime: RuntimeExpression): CompileResult { return { isConst: false, runtime } }
// }

// function wrap<T>(envProvider: (ctx?: C) => Environment, runtime: RuntimeExpression<T>) {
//     return (ctx: C) => runtime(envProvider(ctx));
// }

// function noRuntimeFor(symbol: string) {
//     throw new Error(`Could not find runtime for symbol '${symbol}'.`);
// }

// function applySymbolStatic(runtime: SymbolRuntime, args: RuntimeArguments) {
//     return CompileResult.Dynamic(env => runtime(env, args))
// }

// function applySymbolDynamic(head: RuntimeExpression, args: RuntimeArguments) {
//     return CompileResult.Dynamic(env => {
//         const value = head(env);
//         const symbol = env.runtimeTable[value];
//         if (!symbol) noRuntimeFor(value);
//         return symbol.runtime(env, args);
//     })
// }

// function apply(env: Environment, head: CompileResult, args: RuntimeArguments, constArgs: boolean): CompileResult {
//     if (head.isConst) {
//         const value = head.runtime(env);
//         const symbol = env.runtimeTable[value];
//         if (!symbol) throw new Error(`Could not find runtime for symbol '${value}'.`);
//         if (symbol.attributes.isStatic && constArgs) return CompileResult.Const(symbol.runtime(env, args));
//         return applySymbolStatic(symbol.runtime, args);
//     }
//     return applySymbolDynamic(head.runtime, args);
// }

// function compile(env: Environment, expression: Expression): CompileResult {
//     if (Expression.isLiteral(expression)) {
//         return CompileResult.Const(expression);
//     }

//     if (Expression.isSymbol(expression)) {
//         // TOTO: this needs to look up in the symbol table.
//         return 0 as any;
//     }

//     const head = compile(env, expression.head);
//     if (!expression.args) {
//         return apply(env, head, [], true);
//     } else if (Expression.isArgumentsArray(expression.args)) {
//         const args = [];
//         let constArgs = false;
//         for (const arg of expression.args) {
//             const compiled = compile(env, arg);
//             constArgs = constArgs && compiled.isConst;
//             args.push(compiled.runtime);
//         }
//         return apply(env, head, args as any, constArgs);
//     } else {
//         const args = Object.create(null);
//         let constArgs = false;
//         for (const key of Object.keys(expression.args)) {
//             const compiled = compile(env, expression.args[key]);
//             constArgs = constArgs && compiled.isConst;
//             args[key] = compiled.runtime;
//         }
//         return apply(env, head, args, constArgs);
//     }
// }

// export default Compiler