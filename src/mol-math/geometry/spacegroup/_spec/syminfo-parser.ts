/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../../linear-algebra';
import { coordinateExpressionToOperator, wrap01 } from '../common';
import { SyminfoEntry, syminfoEntries } from './syminfo.lib';

export type { SyminfoEntry };

/**
 * The full operator set for a `SyminfoEntry`: the cross product of its
 * `symops` (point-group coset representatives) and `cenops` (lattice
 * centering translations), translations wrapped into `[0, 1)`.
 */
export function syminfoOperators(entry: SyminfoEntry): Mat4[] {
    const symopMats = entry.symops.map(coordinateExpressionToOperator);
    const cenopMats = entry.cenops.map(coordinateExpressionToOperator);

    const result: Mat4[] = [];
    for (const s of symopMats) {
        for (const c of cenopMats) {
            const op = Mat4.clone(s);
            op[12] = wrap01(op[12] + c[12]);
            op[13] = wrap01(op[13] + c[13]);
            op[14] = wrap01(op[14] + c[14]);
            result.push(op);
        }
    }
    return result;
}

/** All 540 `SyminfoEntry` records, grouped by ITA spacegroup `number`. */
export const SyminfoEntriesByNumber: ReadonlyMap<number, SyminfoEntry[]> = (function () {
    const map = new Map<number, SyminfoEntry[]>();
    for (const entry of syminfoEntries) {
        const list = map.get(entry.number);
        if (list) list.push(entry);
        else map.set(entry.number, [entry]);
    }
    return map;
}());

/** Flattened lookup for the entries with a nonzero `ccp4` code. */
export const SyminfoEntriesByCcp4Number: ReadonlyMap<number, SyminfoEntry> = (function () {
    const map = new Map<number, SyminfoEntry>();
    for (const entry of syminfoEntries) {
        if (entry.ccp4 !== 0) map.set(entry.ccp4, entry);
    }
    return map;
}());
