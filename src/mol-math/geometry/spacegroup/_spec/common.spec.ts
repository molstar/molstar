/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '../../../linear-algebra';
import {
    AxisPermutations, closeGroup, coordinateExpressionToOperator, seitz,
    TetragonalCenteringBasisop, transformOperators, transformTetragonalCentering, wrap01
} from '../common';
import { sameOperatorSet, SyminfoEntriesByCcp4Number, SyminfoEntriesByNumber, syminfoOperators } from './utils';

const IdentityRows = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

describe('wrap01', () => {
    it('leaves in-range twelfths unchanged', () => {
        expect(wrap01(0)).toBe(0);
        expect(wrap01(1 / 2)).toBeCloseTo(1 / 2, 10);
        expect(wrap01(1 / 12)).toBeCloseTo(1 / 12, 10);
        expect(wrap01(5 / 6)).toBeCloseTo(5 / 6, 10);
    });

    it('wraps negative values up into [0, 1)', () => {
        expect(wrap01(-0.25)).toBeCloseTo(0.75, 10);
    });

    it('wraps values >= 1 down into [0, 1)', () => {
        expect(wrap01(1.5)).toBeCloseTo(0.5, 10);
    });

    it('snaps floating-point drift onto the nearest 1/12th', () => {
        expect(wrap01(1 / 3 + 1e-9)).toBeCloseTo(1 / 3, 9);
    });

    it('snaps a value that rounds up to 1 back to 0', () => {
        expect(wrap01(0.999999999)).toBe(0);
    });
});

describe('seitz', () => {
    it('builds the identity operator from identity rotation + zero translation', () => {
        const m = seitz(IdentityRows, [0, 0, 0]);
        expect(Mat4.areEqual(m, Mat4.identity(), 1e-6)).toBe(true);
    });

    it('combines rotation and translation consistently with Vec3.transformMat4', () => {
        // swap x/y, keep z, then translate by (0.5, 0, 0.25)
        const rot = [[0, 1, 0], [1, 0, 0], [0, 0, 1]];
        const m = seitz(rot, [0.5, 0, 0.25]);
        const out = Vec3.transformMat4(Vec3(), Vec3.create(1, 2, 3), m);
        expect(Vec3.equals(out, Vec3.create(2.5, 1, 3.25))).toBe(true);
    });
});

describe('coordinateExpressionToOperator', () => {
    const cases: [string, number[][], number[]][] = [
        ['x,y,z', IdentityRows, [0, 0, 0]],
        ['-x,y,-z', [[-1, 0, 0], [0, 1, 0], [0, 0, -1]], [0, 0, 0]],
        ['-x+1/2,y,-z+1/2', [[-1, 0, 0], [0, 1, 0], [0, 0, -1]], [0.5, 0, 0.5]],
        [' -x + 1/2 , y , -z + 1/2 ', [[-1, 0, 0], [0, 1, 0], [0, 0, -1]], [0.5, 0, 0.5]],
        ['z,x,y', [[0, 0, 1], [1, 0, 0], [0, 1, 0]], [0, 0, 0]],
        [TetragonalCenteringBasisop, [[1, 1, 0], [-1, 1, 0], [0, 0, 1]], [0, 0, 0]],
        ['y,x,-z', [[0, 1, 0], [1, 0, 0], [0, 0, -1]], [0, 0, 0]],
        ['x,z,-y', [[1, 0, 0], [0, 0, 1], [0, -1, 0]], [0, 0, 0]],
    ];

    for (const [expr, rot, tran] of cases) {
        it(`parses '${expr}'`, () => {
            const actual = coordinateExpressionToOperator(expr);
            const expected = seitz(rot, tran);
            expect(Mat4.areEqual(actual, expected, 1e-6)).toBe(true);
        });
    }

    it("matches each AxisPermutations entry's precomputed matrix", () => {
        for (const { basisop, matrix } of AxisPermutations) {
            expect(Mat4.areEqual(coordinateExpressionToOperator(basisop), matrix, 1e-6)).toBe(true);
        }
    });

    it('throws for a triplet without exactly 3 comma-separated parts', () => {
        expect(() => coordinateExpressionToOperator('x,y')).toThrow(/expected 3 comma-separated/);
    });

    it('throws for an empty coordinate expression', () => {
        expect(() => coordinateExpressionToOperator(',,')).toThrow(/empty coordinate expression/);
    });
});

describe('transformOperators', () => {
    // the two symops of spacegroup 3 (P2, 'x,y,z' basisop / canonical setting)
    const canonicalOps = [
        coordinateExpressionToOperator('x,y,z'),
        coordinateExpressionToOperator('-x,y,-z'),
    ];

    it('leaves operators unchanged when P is the identity', () => {
        const out = transformOperators(canonicalOps, Mat4.identity());
        for (let i = 0; i < canonicalOps.length; i++) {
            expect(Mat4.areEqual(out[i], canonicalOps[i], 1e-6)).toBe(true);
        }
    });

    it('is undone by conjugating again with the inverse of P (round-trip)', () => {
        const P = AxisPermutations[4].matrix; // 'z,x,y'
        const Pinv = Mat4.invert(Mat4(), P)!;
        const transformed = transformOperators(canonicalOps, P);
        const roundTrip = transformOperators(transformed, Pinv);
        for (let i = 0; i < canonicalOps.length; i++) {
            expect(Mat4.areEqual(roundTrip[i], canonicalOps[i], 1e-6)).toBe(true);
        }
    });

    it('reproduces syminfo.lib operators for a real alternate axis setting (spacegroup 3, z,x,y)', () => {
        const entries = SyminfoEntriesByNumber.get(3)!;
        const anchor = entries.find(e => e.basisop === 'x,y,z')!;
        const target = entries.find(e => e.basisop === 'z,x,y')!;
        expect(anchor).toBeDefined();
        expect(target).toBeDefined();

        const anchorOps = syminfoOperators(anchor);
        const P = coordinateExpressionToOperator(target.basisop);
        const Pinv = Mat4.invert(Mat4(), P)!;
        const actual = transformOperators(anchorOps, Pinv);
        const expected = syminfoOperators(target);
        expect(sameOperatorSet(actual, expected)).toBe(true);
    });

    it('throws for a non-invertible P', () => {
        const singular = seitz([[0, 0, 0], [1, 0, 0], [0, 1, 0]], [0, 0, 0]);
        expect(() => transformOperators([Mat4.identity()], singular)).toThrow(/non-invertible setting transform/);
    });
});

describe('closeGroup', () => {
    it('returns just the identity for an empty generator set', () => {
        const result = closeGroup([]);
        expect(result.length).toBe(1);
        expect(Mat4.areEqual(result[0], Mat4.identity(), 1e-6)).toBe(true);
    });

    it('closes a single order-2 generator into a 2-element group', () => {
        const gen = seitz([[-1, 0, 0], [0, -1, 0], [0, 0, 1]], [0, 0, 0]); // 2-fold about z
        const result = closeGroup([gen]);
        expect(sameOperatorSet(result, [Mat4.identity(), gen])).toBe(true);
    });

    it('reproduces syminfo.lib operators for a real primitive spacegroup (P4, ccp4 75) from a single generator', () => {
        const entry = SyminfoEntriesByCcp4Number.get(75)!;
        expect(entry).toBeDefined();
        const gen = coordinateExpressionToOperator('-y,x,z'); // 4-fold about z
        const result = closeGroup([gen]);
        const expected = syminfoOperators(entry); // cenops = ['x,y,z'] only (P lattice)
        expect(sameOperatorSet(result, expected)).toBe(true);
    });

    it('throws when the group would exceed maxOrder', () => {
        const gen = coordinateExpressionToOperator('-y,x,z');
        expect(() => closeGroup([gen], 0)).toThrow(/too many elements/);
    });
});

describe('transformTetragonalCentering', () => {
    it('doubles a single identity operator into [identity, +(0.5,0.5,0)]', () => {
        const result = transformTetragonalCentering([Mat4.identity()]);
        expect(result.length).toBe(2);
        expect(Mat4.areEqual(result[0], Mat4.identity(), 1e-6)).toBe(true);
        const expectedShifted = seitz(IdentityRows, [0.5, 0.5, 0]);
        expect(Mat4.areEqual(result[1], expectedShifted, 1e-6)).toBe(true);
    });

    it('doubles a primitive tetragonal (C4) point group and shifts the copy by (0.5, 0.5, 0)', () => {
        const baseOperators = ['x,y,z', '-y,x,z', '-x,-y,z', 'y,-x,z'].map(coordinateExpressionToOperator);
        const result = transformTetragonalCentering(baseOperators);
        expect(result.length).toBe(2 * baseOperators.length);

        const P = coordinateExpressionToOperator(TetragonalCenteringBasisop);
        const expectedFirstHalf = transformOperators(baseOperators, P);
        for (let i = 0; i < baseOperators.length; i++) {
            expect(Mat4.areEqual(result[i], expectedFirstHalf[i], 1e-6)).toBe(true);
        }

        for (let i = 0; i < baseOperators.length; i++) {
            const a = result[i];
            const b = result[baseOperators.length + i];
            // same rotation part
            for (const k of [0, 1, 2, 4, 5, 6, 8, 9, 10]) {
                expect(b[k]).toBeCloseTo(a[k], 6);
            }
            // translation shifted by exactly (0.5, 0.5, 0) mod 1
            expect(b[12]).toBeCloseTo(wrap01(a[12] + 0.5), 6);
            expect(b[13]).toBeCloseTo(wrap01(a[13] + 0.5), 6);
            expect(b[14]).toBeCloseTo(a[14], 6);
        }
    });
});
