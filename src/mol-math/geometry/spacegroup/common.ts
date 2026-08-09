/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../linear-algebra';

/**
 * Wraps a fractional coordinate into the `[0, 1)` interval, snapping to the
 * nearest 1/12th to avoid floating-point drift accumulating (crystallographic
 * translations are always simple fractions like 1/2, 1/3, 1/4, 1/6).
 */
export function wrap01(x: number): number {
    let v = x % 1;
    if (v < 0) v += 1;
    const snapped = Math.round(v * 12) / 12;
    return Math.abs(snapped - 1) < 1e-6 ? 0 : snapped;
}

/**
 * Applies a change-of-basis (conjugation `P⁻¹ · W · P`) to a set of
 * symmetry operators, re-expressing them in the new basis described by `P`
 * (whose columns are the new basis vectors expressed in the old basis).
 * Works for any invertible `P`, including non-unimodular ones (e.g. the
 * tetragonal centering transform, determinant 2).
 */
export function transformOperators(operators: ReadonlyArray<Mat4>, P: Mat4): Mat4[] {
    const Pinv = Mat4();
    if (!Mat4.tryInvert(Pinv, P)) throw new Error('non-invertible setting transform');

    const out: Mat4[] = [];
    for (const op of operators) {
        const tmp = Mat4.mul(Mat4(), op, P);
        const res = Mat4.mul(Mat4(), Pinv, tmp);
        res[12] = wrap01(res[12]);
        res[13] = wrap01(res[13]);
        res[14] = wrap01(res[14]);
        out.push(res);
    }
    return out;
}

/**
 * Builds a fractional-space Seitz operator (`Mat4`) from a 3x3 rotation given
 * as rows together with a translation vector.
 */
export function seitz(rot: ReadonlyArray<ReadonlyArray<number>>, tran: ReadonlyArray<number>): Mat4 {
    return Mat4.ofRows([
        [rot[0][0], rot[0][1], rot[0][2], tran[0]],
        [rot[1][0], rot[1][1], rot[1][2], tran[1]],
        [rot[2][0], rot[2][1], rot[2][2], tran[2]],
        [0, 0, 0, 1],
    ]);
}

function parseFraction(s: string): number {
    const slashIndex = s.indexOf('/');
    if (slashIndex < 0) return parseFloat(s);
    return parseFloat(s.slice(0, slashIndex)) / parseFloat(s.slice(slashIndex + 1));
}

interface AffineRow {
    readonly coeffs: readonly [number, number, number];
    readonly translation: number;
}

/** Parses a single signed affine coordinate expression, e.g. `'-x+y-1/2'` or `'x-y'`. */
function parseAffineExpression(expr: string): AffineRow {
    const cleaned = expr.replace(/\s+/g, '');
    const signed = cleaned[0] === '+' || cleaned[0] === '-' ? cleaned : `+${cleaned}`;
    const terms = signed.match(/[+-][^+-]+/g);
    if (!terms) throw new Error(`empty coordinate expression '${expr}'`);

    const coeffs: [number, number, number] = [0, 0, 0];
    let translation = 0;
    for (const term of terms) {
        const sign = term[0] === '-' ? -1 : 1;
        const body = term.slice(1);
        const axisIndex = 'xyz'.indexOf(body.toLowerCase());
        if (axisIndex >= 0) coeffs[axisIndex] += sign;
        else translation += sign * parseFraction(body);
    }
    return { coeffs, translation };
}

/**
 * Parses a CCP4/ITA `basisop`-style affine coordinate triplet (e.g.
 * `'x,y,z'`, `'-x+1/2,y,-z+1/2'`, `'z,x,y'`) into a fractional-space `Mat4`.
 * Coefficients of `x`/`y`/`z` are always ±1 in this notation.
 */
export function coordinateExpressionToOperator(triplet: string): Mat4 {
    const parts = triplet.split(',');
    if (parts.length !== 3) throw new Error(`expected 3 comma-separated coordinates in '${triplet}'`);
    const rows = parts.map(parseAffineExpression);
    return seitz(rows.map(r => r.coeffs), rows.map(r => r.translation));
}

/**
 * Common denominator for translation parts during group generation. 24 keeps
 * every crystallographic fraction (1/2, 1/3, 1/4, 1/6, 1/8) exact.
 */
const GroupDenominator = 24;

function snapFraction(x: number): number {
    let v = x % 1;
    if (v < 0) v += 1;
    return (Math.round(v * GroupDenominator) % GroupDenominator) / GroupDenominator;
}

function wrapTranslation(op: Mat4): Mat4 {
    const m = Mat4.clone(op);
    m[12] = snapFraction(m[12]);
    m[13] = snapFraction(m[13]);
    m[14] = snapFraction(m[14]);
    return m;
}

function operatorKey(op: Mat4): string {
    const parts: number[] = [];
    for (const i of [0, 1, 2, 4, 5, 6, 8, 9, 10]) parts.push(Math.round(op[i]));
    for (const i of [12, 13, 14]) {
        let v = op[i] % 1;
        if (v < 0) v += 1;
        parts.push(Math.round(v * GroupDenominator) % GroupDenominator);
    }
    return parts.join(',');
}

/**
 * Closes a set of generator operators into the full group by repeated
 * multiplication (orbit of the identity), wrapping translations into `[0, 1)`
 * at 1/24 resolution. Throws if the group exceeds `maxOrder` (which indicates
 * bad or inconsistent generators).
 */
export function closeGroup(generators: ReadonlyArray<Mat4>, maxOrder = 1024): Mat4[] {
    const identity = Mat4.identity();
    const result: Mat4[] = [identity];
    const seen = new Set<string>([operatorKey(identity)]);
    const gens = generators.map(wrapTranslation);
    for (let i = 0; i < result.length; i++) {
        for (const g of gens) {
            const prod = wrapTranslation(Mat4.mul(Mat4(), result[i], g));
            const key = operatorKey(prod);
            if (!seen.has(key)) {
                seen.add(key);
                result.push(prod);
                if (result.length > maxOrder) throw new Error('closeGroup: too many elements - bad generators');
            }
        }
    }
    return result;
}

/**
 * A standard ITA axis-permutation change-of-basis: `basisop` is the
 * xyz-style affine coordinate expression (same notation as
 * `SpacegroupEntry.basisop` in `./tables`, e.g. `'y,x,-z'`) and `matrix` is
 * the corresponding change-of-basis `Mat4` (derived from `basisop` via
 * `coordinateExpressionToOperator`).
 */
export interface AxisPermutation {
    readonly basisop: string,
    readonly matrix: Mat4
}

/**
 * The six standard axis settings of the monoclinic and orthorhombic space
 * groups (ITfC Vol. A Table 1.5.4.4), expressed as xyz-style `basisop`
 * coordinate expressions rather than the classic basis-vector notation
 * (`'abc'`, `'ba-c'`, ...) - e.g. `'y,x,-z'` is the `'ba-c'` setting. A
 * non-default axis setting's symmetry operators are obtained by conjugating
 * the canonical (setting `'x,y,z'`) operators with the matching `matrix`
 * (see `createWithSetting` + `transformOperators`) - this is standard ITA
 * data, independent of any particular file format's setting-numbering
 * convention.
 */
export const AxisPermutations: readonly AxisPermutation[] = ([
    'x,y,z',
    'y,x,-z',
    'y,z,x',
    'z,y,-x',
    'z,x,y',
    'x,z,-y',
] as const).map(basisop => ({ basisop, matrix: coordinateExpressionToOperator(basisop) }));

/**
 * Standard ITA change-of-basis `basisop` for the doubled-cell "C"/"F"
 * description of a primitive tetragonal lattice (det = 2): `'x+y,-x+y,z'`
 * (a general linear combination, not a pure axis permutation - so it isn't
 * one of the six `AxisPermutations` entries).
 */
export const TetragonalCenteringBasisop = 'x+y,-x+y,z';

/** a'=a-b, b'=a+b, c'=c: the doubled-cell "C" description of a primitive tetragonal lattice (det = 2). */
const TetragonalCentering = coordinateExpressionToOperator(TetragonalCenteringBasisop);

/** Doubled-cell "C"/"F"-centered description of a primitive tetragonal group (setting 2 only). */
export function transformTetragonalCentering(baseOperators: ReadonlyArray<Mat4>): Mat4[] {
    const transformed = transformOperators(baseOperators, TetragonalCentering);
    const out = transformed.slice();
    for (const op of transformed) {
        const copy = Mat4.clone(op);
        copy[12] = wrap01(copy[12] + 0.5);
        copy[13] = wrap01(copy[13] + 0.5);
        out.push(copy);
    }
    return out;
}
