/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, ReadonlyMat4 } from '../../../linear-algebra/3d/mat4';
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

/**
 * Canonical key for a symmetry operator, rounding the rotation part to the
 * nearest integer and the translation part to the nearest 1/24 (see
 * `common.ts`'s `GroupDenominator`), so operators that differ only by
 * floating-point drift compare equal.
 */
export function opKey(op: ReadonlyMat4): string {
    const parts: number[] = [];
    for (const i of [0, 1, 2, 4, 5, 6, 8, 9, 10]) parts.push(Math.round(op[i]));
    for (const i of [12, 13, 14]) {
        let v = op[i] % 1;
        if (v < 0) v += 1;
        parts.push(Math.round(v * 24) % 24);
    }
    return parts.join(',');
}

/** Whether two operator lists contain exactly the same operators (as multisets), via `opKey`. */
export function sameOperatorSet(a: ReadonlyArray<ReadonlyMat4>, b: ReadonlyArray<ReadonlyMat4>) {
    if (a.length !== b.length) return false;
    const ka = a.map(opKey).sort();
    const kb = b.map(opKey).sort();
    return ka.every((k, i) => k === kb[i]);
}

/**
 * Like `sameOperatorSet`, but ignores duplicate multiplicity. Needed for
 * basisop transforms with a non-unimodular determinant (e.g. the R/H
 * rhombohedral<->hexagonal settings, determinant 3): conjugating the
 * redundant hex-centering copies of an operator legitimately collapses onto
 * the same matrix once wrapped into the smaller rhombohedral cell.
 */
export function sameOperatorKeySet(a: ReadonlyArray<ReadonlyMat4>, b: ReadonlyArray<ReadonlyMat4>) {
    const ka = new Set(a.map(opKey));
    const kb = new Set(b.map(opKey));
    if (ka.size !== kb.size) return false;
    for (const k of ka) if (!kb.has(k)) return false;
    return true;
}
