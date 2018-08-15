// /**
//  * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
//  *
//  * @author David Sehnal <david.sehnal@gmail.com>
//  */

// import { MSymbol } from '../language/symbol'
// import { SymbolRuntime } from './symbol'
// import { Macro } from './macro';
// import Expression from '../language/expression';

// class Environment {
//     readonly runtimeTable: SymbolRuntime.Table;
//     readonly macroTable: Macro.Table = new Map<string, Macro>();

//     addMacro(name: string, expression: Expression, argNames: ReadonlyArray<string>): Macro {
//         const argIndex: Macro['argIndex'] = {};
//         for (let i = 0; i < argNames.length; i++) argIndex[argNames[i]] = i;
//         const macro: Macro = { expression, argIndex, argNames };
//         this.macroTable.set(name, macro);
//         return macro;
//     }

//     removeMacro(name: string) {
//         this.macroTable.delete(name);
//     }

//     addSymbolRuntime(symbol: MSymbol, runtime: SymbolRuntime) {
//         this.runtimeTable.set(symbol.id, runtime);
//     }

//     removeSymbolRuntime(symbol: MSymbol) {
//         this.runtimeTable.delete(symbol.id);
//     }
// }

// export default Environment