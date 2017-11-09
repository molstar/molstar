/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3, Mat4 } from '../linear-algebra/3d'

interface SymmetryOperator {
    readonly name: string,
    readonly hkl: Vec3,

    readonly matrix: Mat4,
    // cache the inverse of the transform
    readonly inverse: Mat4,
    // optimize the identity case
    readonly isIdentity: boolean
}

namespace SymmetryOperator {
    export const DefaultName = '1_555'
    export const Default: SymmetryOperator = create(DefaultName, Mat4.identity());

    const RotationEpsilon = 0.0001;

    export function create(name: string, matrix: Mat4, hkl?: number[]): SymmetryOperator {
        const _hkl = hkl ? Vec3.create(hkl[0], hkl[1], hkl[2]) : Vec3.zero();
        if (Mat4.isIdentity(matrix)) return { name, matrix, inverse: Mat4.identity(), isIdentity: true, hkl: _hkl };
        if (!Mat4.isRotationAndTranslation(matrix, RotationEpsilon)) throw new Error(`Symmetry operator (${name}) must be a composition of rotation and translation.`);
        return { name, matrix, inverse: Mat4.invert(Mat4.zero(), matrix), isIdentity: false, hkl: _hkl };
    }

    export interface CoordinateMapper { (index: number, slot: Vec3): Vec3 }
    export interface ArrayMapping {
        readonly invariantPosition: CoordinateMapper,
        readonly position: CoordinateMapper,
        x(index: number): number,
        y(index: number): number,
        z(index: number): number
    }

    export interface Coordinates { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> }

    export function createMapping(operator: SymmetryOperator, coords: Coordinates) {
        const invariantPosition = SymmetryOperator.createCoordinateMapper(SymmetryOperator.Default, coords);
        const position = operator.isIdentity ? invariantPosition : SymmetryOperator.createCoordinateMapper(operator, coords);
        const { x, y, z } = createProjections(operator, coords);
        return { invariantPosition, position, x, y, z };
    }

    export function createCoordinateMapper(t: SymmetryOperator, coords: Coordinates): CoordinateMapper {
        if (t.isIdentity) return identityPosition(coords);
        return generalPosition(t, coords);
    }
}

export default SymmetryOperator

interface Projections { x(index: number): number, y(index: number): number, z(index: number): number }

function createProjections(t: SymmetryOperator, coords: SymmetryOperator.Coordinates): Projections {
    if (t.isIdentity) return { x: projectCoord(coords.x), y: projectCoord(coords.y), z: projectCoord(coords.z) };
    return { x: projectX(t, coords), y: projectY(t, coords), z: projectZ(t, coords) };
}

function projectCoord(xs: ArrayLike<number>) {
    return (i: number) => xs[i];
}

function isW1(m: Mat4) {
    return m[3] === 0 && m[7] === 0 && m[11] === 0 && m[15] === 1;
}

function projectX({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs}: SymmetryOperator.Coordinates) {
    const xx = m[0], yy = m[4], zz = m[8], tx = m[12];

    if (isW1(m)) {
        // this should always be the case.
        return (i: number) => xx * xs[i] + yy * ys[i] + zz * zs[i] + tx;
    }

    return (i: number) => {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + tx) / w;
    }
}

function projectY({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs}: SymmetryOperator.Coordinates) {
    const xx = m[1], yy = m[5], zz = m[9], ty = m[13];

    if (isW1(m)) {
        // this should always be the case.
        return (i: number) => xx * xs[i] + yy * ys[i] + zz * zs[i] + ty;
    }

    return (i: number) => {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + ty) / w;
    }
}

function projectZ({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs}: SymmetryOperator.Coordinates) {
    const xx = m[2], yy = m[6], zz = m[10], tz = m[14];

    if (isW1(m)) {
        // this should always be the case.
        return (i: number) => xx * xs[i] + yy * ys[i] + zz * zs[i] + tz;
    }

    return (i: number) => {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + tz) / w;
    }
}

function identityPosition({ x, y, z }: SymmetryOperator.Coordinates): SymmetryOperator.CoordinateMapper {
    return (i, s) => {
        s[0] = x[i];
        s[1] = y[i];
        s[2] = z[i];
        return s;
    }
}

function generalPosition({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs }: SymmetryOperator.Coordinates) {
    if (isW1(m)) {
        // this should always be the case.
        return (i: number, r: Vec3): Vec3 => {
            const x = xs[i], y = ys[i], z = zs[i];
            r[0] = m[0] * x + m[4] * y + m[8] * z + m[12];
            r[1] = m[1] * x + m[5] * y + m[9] * z + m[13];
            r[2] = m[2] * x + m[6] * y + m[10] * z + m[14];
            return r;
        }
    }
    return (i: number, r: Vec3): Vec3 => {
        r[0] = xs[i];
        r[1] = ys[i];
        r[2] = zs[i];
        Vec3.transformMat4(r, r, m);
        return r;
    }
}