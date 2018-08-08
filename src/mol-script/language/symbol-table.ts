/**
 * Copyright (c) 2017 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import core from './symbol-table/core'
import structureQuery from './symbol-table/structure-query'
import { normalizeTable, symbolList } from './helpers'
import { MSymbol } from './symbol'

const MolScriptSymbolTable = { core, structureQuery };

normalizeTable(MolScriptSymbolTable);

export const SymbolList = symbolList(MolScriptSymbolTable);

export const SymbolMap = (function() {
    const map: { [id: string]: MSymbol | undefined } = Object.create(null);
    for (const s of SymbolList) map[s.id] = s;
    return map;
})();

export { MolScriptSymbolTable }