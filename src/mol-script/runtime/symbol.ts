/*
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Environment from './environment'
import RuntimeExpression from './expression'

export type RuntimeArguments = ArrayLike<RuntimeExpression> | { [name: string]: RuntimeExpression | undefined }

type SymbolRuntime = (env: Environment, args: RuntimeArguments) => any

namespace SymbolRuntime {
    export interface Info {
        readonly runtime: SymbolRuntime,
        readonly attributes: Attributes
    }

    export interface Attributes { isStatic: boolean }
}

function SymbolRuntime(symbol: Symbol, attributes: Partial<SymbolRuntime.Attributes> = {}) {
    const { isStatic = false } = attributes;
    return (runtime: SymbolRuntime): SymbolRuntime.Info => {
        return ({ runtime, attributes: { isStatic } });
    };
}

export type SymbolRuntimeTable = { readonly [id: string]: SymbolRuntime.Info }

export default SymbolRuntime