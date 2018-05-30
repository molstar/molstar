/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Expression from './expression'
import Symbol from './symbol'
import SymbolTable from './symbol-table'

namespace Builder {
    export const core = SymbolTable.core;
    export const struct = SymbolTable.structureQuery;

    export function atomName(s: string) { return struct.type.atomName([s]); }
    export function es(s: string) { return struct.type.elementSymbol([s]); }
    export function list(...xs: Expression[]) { return core.type.list(xs); }
    export function set(...xs: Expression[]) { return core.type.set(xs); }
    export function fn(x: Expression) { return core.ctrl.fn([x]); }
    export function evaluate(x: Expression) { return core.ctrl.eval([x]); }

    const _acp = struct.atomProperty.core, _ammp = struct.atomProperty.macromolecular, _atp = struct.atomProperty.topology;

    // atom core property
    export function acp(p: keyof typeof _acp) { return (_acp[p] as Symbol<any>)() };

    // atom topology property
    export function atp(p: keyof typeof _atp) { return (_atp[p] as Symbol<any>)() };

    // atom macromolecular property
    export function ammp(p: keyof typeof _ammp) { return (_ammp[p] as Symbol<any>)() };

    // atom property sets
    const _aps = struct.atomSet.propertySet
    export function acpSet(p: keyof typeof _acp) { return _aps([ acp(p) ]) };
    export function atpSet(p: keyof typeof _atp) { return _aps([ atp(p) ]) };
    export function ammpSet(p: keyof typeof _ammp) { return _aps([ ammp(p) ]) };
}

export default Builder