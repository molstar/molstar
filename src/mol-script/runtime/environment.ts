/*
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SymbolRuntimeTable } from './symbol'

interface Environment<T = any> {
    readonly symbolTable: SymbolRuntimeTable,
    readonly context: T
}

export default Environment