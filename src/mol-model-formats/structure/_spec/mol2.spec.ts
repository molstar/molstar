/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseMol2 } from '../../../mol-io/reader/mol2/parser';
import { transformOperators } from '../../../mol-math/geometry/spacegroup/common';
import { Spacegroup, SpacegroupCell } from '../../../mol-math/geometry/spacegroup/construction';
import { Mat4 } from '../../../mol-math/linear-algebra/3d/mat4';
import { Vec3 } from '../../../mol-math/linear-algebra/3d/vec3';
import { getOperatorsForSetting, SettingsByNumber, trajectoryFromMol2 } from '../mol2';
import { ModelSymmetry } from '../property/symmetry';

// real CSD-derived CRYSIN record (spaceGroup 29, setting 5), taken from
// src/mol-io/reader/_spec/mol2.spec.ts 'crysin' test fixture
const Mol2StringCrysin = `@<TRIPOS>MOLECULE
1144204
    12    11     2     0     0
SMALL
USER_CHARGES
****
Generated from the CSD

@<TRIPOS>ATOM
     1 Cl1      0.0925   3.6184   1.9845   Cl        1 RES1  -1.0000
     2 C1      -4.7391   0.3350   0.4215   C.ar      2 RES2   0.0000
     3 C2      -3.4121   0.2604   0.9351   C.ar      2 RES2   0.0000
     4 C3      -2.9169   1.2555   1.7726   C.ar      2 RES2   0.0000
     5 C4      -3.7118   2.3440   2.1099   C.ar      2 RES2   0.0000
     6 C5      -5.0314   2.4052   1.6209   C.ar      2 RES2   0.0000
     7 C6      -5.5372   1.4057   0.7962   C.ar      2 RES2   0.0000
     8 C7      -6.9925   1.4547   0.3334   C.3       2 RES2   0.0000
     9 C8      -7.8537   0.5554   1.1859   C.3       2 RES2   0.0000
    10 N1      -9.3089   0.7134   0.8192   N.3       2 RES2   1.0000
    11 O1      -2.6613  -0.8147   0.5707   O.3       2 RES2   0.0000
    12 O2      -1.6204   1.0919   2.2584   O.3       2 RES2   0.0000
@<TRIPOS>BOND
     1     2     3   ar
     2     3     4   ar
     3     4     5   ar
     4     5     6   ar
     5     6     7   ar
     6     7     2   ar
     7     8     7    1
     8     9     8    1
     9    10     9    1
    10    11     3    1
    11    12     4    1
@<TRIPOS>SUBSTRUCTURE
     1 RES1        1 GROUP             0 ****  ****    0
     2 RES2        2 GROUP             0 ****  ****    0
@<TRIPOS>CRYSIN
   10.5150   11.1300    7.9380   90.0000   90.0000   90.0000    29     5
`;


const Angles90 = Vec3.create(Math.PI / 2, Math.PI / 2, Math.PI / 2);

function baseOperators(spacegroupNumber: number) {
    const cell = SpacegroupCell.create(spacegroupNumber, Vec3.create(1, 1, 1), Angles90);
    return Spacegroup.create(cell).operators;
}

function operatorsByName(name: string) {
    const cell = SpacegroupCell.create(name as any, Vec3.create(1, 1, 1), Angles90);
    return Spacegroup.create(cell).operators;
}

/** Compares two operator lists as sets (order independent), matching translations up to 1/12. */
function operatorsEqualAsSets(a: ReadonlyArray<ArrayLike<number>>, b: ReadonlyArray<ArrayLike<number>>) {
    if (a.length !== b.length) return false;
    const used = new Array(b.length).fill(false);
    for (const opA of a) {
        const idx = b.findIndex((opB, i) => !used[i] && Array.from(opA).every((v: any, k: number) => Math.abs(v - opB[k]) < 1e-4));
        if (idx < 0) return false;
        used[idx] = true;
    }
    return true;
}

// The six proper axis-permutation change-of-basis matrices relating the standard
// monoclinic/orthorhombic settings (ITfC Vol. A Table 1.5.4.4).
const AxisPermutations: Mat4[] = [
    Mat4.identity(),
    Mat4.ofRows([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]),
    Mat4.ofRows([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
    Mat4.ofRows([[0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1]]),
    Mat4.ofRows([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
    Mat4.ofRows([[1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]]),
];

/** True if `settingOps` is an axis-permutation of the canonical operators (i.e. a genuine alternate setting). */
function isAxisPermutationOfCanonical(settingOps: ReadonlyArray<Mat4>, canonicalOps: ReadonlyArray<Mat4>) {
    return AxisPermutations.some(P => operatorsEqualAsSets(transformOperators(canonicalOps, P), settingOps));
}

describe('mol2 crysin setting', () => {
    it('resolves spacegroup 29, setting 5 to the "P b c 21" alternate description', async () => {
        const parsed = await parseMol2(Mol2StringCrysin, '').run();
        if (parsed.isError) throw new Error(parsed.message);

        const trajectory = await trajectoryFromMol2(parsed.result).run();
        const model = trajectory.representative;

        const symmetry = ModelSymmetry.Provider.get(model);
        expect(symmetry).toBeDefined();
        expect(symmetry!.spacegroup.name).toBe('P b c 21');
        expect(symmetry!.spacegroup.num).toBe(29);
        expect(symmetry!.spacegroup.operators.length).toBe(4);
    });

    it('resolves the real-world CSD case: spacegroup 29, setting 5 -> "P b c 21"', () => {
        const result = getOperatorsForSetting(29, 5, baseOperators(29));
        expect(result).toBeDefined();
        expect(result!.name).toBe('P b c 21');
        expect(result!.basisop).toBe('y,x,-z');
        expect(result!.operators.length).toBe(4);
    });

    it('matches existing hardcoded alternate settings (calibration / regression)', () => {
        // spg 17 'P 2 2 21' -> settings 2 & 3 (existing GroupData indices 256/257)
        {
            const base = baseOperators(17);
            const s2 = getOperatorsForSetting(17, 2, base)!;
            expect(s2.name).toBe('P 21 2 2');
            expect(operatorsEqualAsSets(s2.operators, operatorsByName('P 21 2 2'))).toBe(true);

            const s3 = getOperatorsForSetting(17, 3, base)!;
            expect(s3.name).toBe('P 2 21 2');
            expect(operatorsEqualAsSets(s3.operators, operatorsByName('P 2 21 2'))).toBe(true);
        }
        // spg 18 'P 21 21 2' -> settings 2 & 3 (existing GroupData indices 260/259)
        {
            const base = baseOperators(18);
            const s2 = getOperatorsForSetting(18, 2, base)!;
            expect(s2.name).toBe('P 2 21 21');
            expect(operatorsEqualAsSets(s2.operators, operatorsByName('P 2 21 21'))).toBe(true);

            const s3 = getOperatorsForSetting(18, 3, base)!;
            expect(s3.name).toBe('P 21 2 21');
            expect(operatorsEqualAsSets(s3.operators, operatorsByName('P 21 2 21'))).toBe(true);
        }
        // spg 5 'C 1 2 1' -> setting 3 "B 1 1 2" (existing GroupData index 239)
        {
            const base = baseOperators(5);
            const s3 = getOperatorsForSetting(5, 3, base)!;
            expect(s3.name).toBe('B 1 1 2');
            expect(operatorsEqualAsSets(s3.operators, operatorsByName('B 1 1 2'))).toBe(true);
        }
        // spg 3 'P 1 2 1' -> setting 3 "P 1 1 2" (existing GroupData index 237)
        {
            const base = baseOperators(3);
            const s3 = getOperatorsForSetting(3, 3, base)!;
            expect(s3.name).toBe('P 1 1 2');
            expect(operatorsEqualAsSets(s3.operators, operatorsByName('P 1 1 2'))).toBe(true);
        }
    });

    it('every non-blank cell in SettingsByNumber for spg 3-74 resolves to a correct axis-permuted operator set', () => {
        const failures: string[] = [];
        for (let spg = 3; spg <= 74; spg++) {
            const row = SettingsByNumber[spg];
            if (!row) continue;
            const base = baseOperators(spg);
            for (let setting = 2; setting <= 6; setting++) {
                const name = row[setting - 2]?.[0];
                if (!name) continue;

                const result = getOperatorsForSetting(spg, setting, base);
                if (!result) {
                    failures.push(`spg ${spg} setting ${setting} (${name}): no matching transform found`);
                } else if (result.operators.length !== base.length) {
                    failures.push(`spg ${spg} setting ${setting} (${name}): operator count mismatch (${result.operators.length} vs ${base.length})`);
                } else if (!isAxisPermutationOfCanonical(result.operators, base)) {
                    failures.push(`spg ${spg} setting ${setting} (${name}): operators are not an axis-permutation of the canonical group`);
                }
            }
        }
        expect(failures).toEqual([]);
    });

    it('tetragonal setting 2 doubles operator count with C-centering', () => {
        const base = baseOperators(75); // 'P 4'
        const result = getOperatorsForSetting(75, 2, base)!;
        expect(result.name).toBe('C 4');
        expect(result.basisop).toBe('x+y,-x+y,z');
        expect(result.operators.length).toBe(base.length * 2);
    });
});
