/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../../linear-algebra/3d/mat4';
import { coordinateExpressionToOperator, transformOperators } from '../common';
import { operatorsFromHall } from '../notation';
import { getHallSymbol, SpacegroupData } from '../tables';
import { sameOperatorKeySet, sameOperatorSet, SyminfoEntriesByCcp4Number, SyminfoEntriesByNumber, syminfoOperators } from './utils';

/**
 * Reconstructs operators directly from CCP4 `syminfo.lib` (see `./syminfo.lib`), kept as an independent regression oracle now that
 * Spacegroup.create itself derives operators from Hall symbols (so comparing
 * against Spacegroup.create would be circular).
 */
function syminfoOperatorsFor(spacegroupNumber: number): Mat4[] {
    const entry = SyminfoEntriesByCcp4Number.get(spacegroupNumber);
    if (!entry) throw new Error(`no syminfo.lib entry with ccp4 number ${spacegroupNumber}`);
    return syminfoOperators(entry);
}

/**
 * (molstar spacegroup number, Hall symbol) — a diverse hand-picked sample used
 * to validate the interpreter before the full 230-entry table is introduced.
 */
const Samples: [number, string][] = [
    [1, 'P 1'],
    [2, '-P 1'],
    [4, 'P 2yb'], // P 1 21 1 (monoclinic screw)
    [5, 'C 2y'], // C 1 2 1
    [14, '-P 2ybc'], // P 21/c
    [19, 'P 2ac 2ab'], // P 21 21 21
    [29, 'P 2c -2ac'], // P c a 21 (non-symmorphic ortho)
    [43, 'F 2 -2d'], // F d d 2 (d-glide)
    [48, 'P 2 2 -1n'], // P n n n (origin choice 1)
    [70, 'F 2 2 -1d'], // F d d d (origin choice 1)
    [92, 'P 4abw 2nw'], // P 41 21 2 (tetragonal screw)
    [85, 'P 4ab -1ab'], // P 4/n (origin choice 1)
    [143, 'P 3'],
    [146, 'R 3'], // molstar 'H 3' hexagonal R-centred setting
    [151, 'P 31 2 (0 0 4)'], // P 31 1 2 (origin shift)
    [178, 'P 61 2 (0 0 5)'], // P 61 2 2 (origin shift)
    [198, 'P 2ac 2ab 3'], // P 21 3 (cubic)
    [227, 'F 4d 2 3 -1d'], // F d -3 m (origin choice 1)
    [230, '-I 4bd 2c 3'], // I a -3 d
];

describe('spacegroup Hall interpreter', () => {
    it('operatorsFromHall matches syminfo.lib operators for a diverse sample', () => {
        const failures: string[] = [];
        for (const [num, hall] of Samples) {
            const expected = syminfoOperatorsFor(num);
            let actual;
            try {
                actual = operatorsFromHall(hall);
            } catch (e) {
                failures.push(`spg ${num} '${hall}': threw ${e}`);
                continue;
            }
            if (!sameOperatorSet(actual, expected)) {
                failures.push(`spg ${num} '${hall}': ${actual.length} ops vs ${expected.length} expected`);
            }
        }
        expect(failures).toEqual([]);
    });
});

describe('spacegroup Hall symbol table', () => {
    it('getHallSymbol + operatorsFromHall reproduces syminfo.lib operators for every numbered spacegroup in SpacegroupData', () => {
        const failures: string[] = [];
        for (const entry of SpacegroupData) {
            if (entry.ccp4Number === 0) continue;
            const number = entry.ccp4Number;
            const hall = getHallSymbol(number);
            if (!hall) {
                failures.push(`num ${number}: no Hall symbol`);
                continue;
            }
            const expected = syminfoOperatorsFor(number);
            let actual;
            try {
                actual = operatorsFromHall(hall);
            } catch (e) {
                failures.push(`num ${number}, '${hall}': threw ${e}`);
                continue;
            }
            if (!sameOperatorSet(actual, expected)) {
                failures.push(`num ${number}, '${hall}': ${actual.length} ops vs ${expected.length} expected`);
            }
        }
        expect(failures).toEqual([]);
    });
});

describe('spacegroup Hall symbol table (alternate settings)', () => {
    it('transforming the reference setting by basisop reproduces syminfo.lib operators for every setting, including non-CCP4 (ccp4=0) ones', () => {
        const failures: string[] = [];
        for (const [number, entries] of SyminfoEntriesByNumber) {
            // Prefer the identity-basisop ('x,y,z') sibling as the anchor; for
            // origin-choice spacegroups (e.g. 48, 201) that sibling has no
            // CCP4 number (ccp4=0, not in SpacegroupData), so fall back to
            // any sibling that does have a resolvable Hall symbol.
            const identity = entries.find(e => e.basisop === 'x,y,z');
            const anchor = (identity && getHallSymbol(identity.ccp4))
                ? identity
                : entries.find(e => e.ccp4 !== 0 && getHallSymbol(e.ccp4));
            if (!anchor) {
                failures.push(`num ${number}: no sibling setting with a resolvable Hall symbol`);
                continue;
            }

            let canonicalOps: Mat4[];
            try {
                const anchorOps = operatorsFromHall(getHallSymbol(anchor.ccp4)!) as Mat4[];
                const anchorP = coordinateExpressionToOperator(anchor.basisop);
                // Bring anchorOps back into the common (basisop-identity) frame:
                // `transformOperators(ops, Q)` computes Q⁻¹·op·Q, so passing the
                // anchor's own basisop as `Q` undoes it.
                canonicalOps = transformOperators(anchorOps, anchorP);
            } catch (e) {
                failures.push(`num ${number}: failed to build reference operators from anchor ccp4 ${anchor.ccp4}: ${e}`);
                continue;
            }

            for (const entry of entries) {
                try {
                    const P = coordinateExpressionToOperator(entry.basisop);
                    const Pinv = Mat4.invert(Mat4(), P);
                    if (!Pinv) throw new Error('non-invertible basisop');
                    // Re-apply forward via the inverse (see note above): this
                    // conjugates canonicalOps by `P`, landing in entry's frame.
                    const actual = transformOperators(canonicalOps, Pinv);
                    const expected = syminfoOperators(entry);
                    if (!sameOperatorKeySet(actual, expected)) {
                        failures.push(`num ${number} setting '${entry.basisop}' (ccp4 ${entry.ccp4}, hall '${entry.hall}'): ${actual.length} ops vs ${expected.length} expected`);
                    }
                } catch (e) {
                    failures.push(`num ${number} setting '${entry.basisop}' (ccp4 ${entry.ccp4}, hall '${entry.hall}'): threw ${e}`);
                }
            }
        }
        expect(failures).toEqual([]);
    });
});
