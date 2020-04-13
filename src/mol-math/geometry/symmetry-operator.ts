/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3, Mat4, Mat3, Quat } from '../linear-algebra/3d';
import { lerp as scalar_lerp } from '../../mol-math/interpolate';
import { defaults } from '../../mol-util';

interface SymmetryOperator {
    readonly name: string,

    readonly assembly?: {
        /** pointer to `pdbx_struct_assembly.id` or empty string */
        readonly id: string,
        /** pointers to `pdbx_struct_oper_list.id` or empty list */
        readonly operList: string[],
        /** (arbitrary) unique id of the operator to be used in suffix */
        readonly operId: number
    },

    /** pointer to `struct_ncs_oper.id` or empty string */
    readonly ncsId: string,

    readonly hkl: Vec3,
    /** spacegroup symmetry operator index, -1 if not applicable */
    readonly spgrOp: number,

    readonly matrix: Mat4,
    // cache the inverse of the transform
    readonly inverse: Mat4,
    // optimize the identity case
    readonly isIdentity: boolean,

    /**
     * Suffix based on operator type.
     * - Assembly: _assembly.operId
     * - Crytal: -op_ijk
     * - ncs: _ncsId
     */
    readonly suffix: string
}

namespace SymmetryOperator {
    export const DefaultName = '1_555';
    export const Default: SymmetryOperator = create(DefaultName, Mat4.identity());

    export const RotationTranslationEpsilon = 0.005;

    export type CreateInfo = { assembly?: SymmetryOperator['assembly'], ncsId?: string, hkl?: Vec3, spgrOp?: number }
    export function create(name: string, matrix: Mat4, info?: CreateInfo): SymmetryOperator {
        let { assembly, ncsId, hkl, spgrOp } = info || { };
        const _hkl = hkl ? Vec3.clone(hkl) : Vec3.zero();
        spgrOp = defaults(spgrOp, -1);
        ncsId = ncsId || '';
        const suffix = getSuffix(info);
        if (Mat4.isIdentity(matrix)) return { name, assembly, matrix, inverse: Mat4.identity(), isIdentity: true, hkl: _hkl, spgrOp, ncsId, suffix };
        if (!Mat4.isRotationAndTranslation(matrix, RotationTranslationEpsilon)) throw new Error(`Symmetry operator (${name}) must be a composition of rotation and translation.`);
        return { name, assembly, matrix, inverse: Mat4.invert(Mat4.zero(), matrix), isIdentity: false, hkl: _hkl, spgrOp, ncsId, suffix };
    }

    function getSuffix(info?: CreateInfo) {
        if (!info) return '';

        if (info.assembly) {
            return `_${info.assembly.operId}`;
        }

        if (typeof info.spgrOp !== 'undefined' && typeof info.hkl !== 'undefined' && info.spgrOp !== -1) {
            const [i, j, k] = info.hkl;
            return `-${info.spgrOp + 1}_${5 + i}${5 + j}${5 + k}`;
        }

        if (info.ncsId) {
            return `_${info.ncsId}`;
        }

        return '';
    }

    export function checkIfRotationAndTranslation(rot: Mat3, offset: Vec3) {
        const matrix = Mat4.identity();
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                Mat4.setValue(matrix, i, j, Mat3.getValue(rot, i, j));
            }
        }
        Mat4.setTranslation(matrix, offset);
        return Mat4.isRotationAndTranslation(matrix, RotationTranslationEpsilon);
    }

    export function ofRotationAndOffset(name: string, rot: Mat3, offset: Vec3, ncsId?: string) {
        const t = Mat4.identity();
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                Mat4.setValue(t, i, j, Mat3.getValue(rot, i, j));
            }
        }
        Mat4.setTranslation(t, offset);
        return create(name, t, { ncsId });
    }

    const _q1 = Quat.identity(), _q2 = Quat.zero(), _q3 = Quat.zero(), _axis = Vec3.zero();
    export function lerpFromIdentity(out: Mat4, op: SymmetryOperator, t: number): Mat4 {
        const m = op.inverse;
        if (op.isIdentity) return Mat4.copy(out, m);

        const _t = 1 - t;
        // interpolate rotation
        Mat4.getRotation(_q2, m);
        Quat.slerp(_q2, _q1, _q2, _t);
        const angle = Quat.getAxisAngle(_axis, _q2);
        Mat4.fromRotation(out, angle, _axis);

        // interpolate translation
        Mat4.setValue(out, 0, 3, _t * Mat4.getValue(m, 0, 3));
        Mat4.setValue(out, 1, 3, _t * Mat4.getValue(m, 1, 3));
        Mat4.setValue(out, 2, 3, _t * Mat4.getValue(m, 2, 3));

        return out;
    }

    export function slerp(out: Mat4, src: Mat4, tar: Mat4, t: number): Mat4 {
        if (Math.abs(t) <= 0.00001) return Mat4.copy(out, src);
        if (Math.abs(t - 1) <= 0.00001) return Mat4.copy(out, tar);

        // interpolate rotation
        Mat4.getRotation(_q2, src);
        Mat4.getRotation(_q3, tar);
        Quat.slerp(_q3, _q2, _q3, t);
        const angle = Quat.getAxisAngle(_axis, _q3);
        Mat4.fromRotation(out, angle, _axis);

        // interpolate translation
        Mat4.setValue(out, 0, 3, scalar_lerp(Mat4.getValue(src, 0, 3), Mat4.getValue(tar, 0, 3), t));
        Mat4.setValue(out, 1, 3, scalar_lerp(Mat4.getValue(src, 1, 3), Mat4.getValue(tar, 1, 3), t));
        Mat4.setValue(out, 2, 3, scalar_lerp(Mat4.getValue(src, 2, 3), Mat4.getValue(tar, 2, 3), t));

        return out;
    }

    /**
     * Apply the 1st and then 2nd operator. ( = second.matrix * first.matrix).
     * Keep `name`, `assembly`, `ncsId`, `hkl` and `spgrOpId` properties from second.
     */
    export function compose(first: SymmetryOperator, second: SymmetryOperator) {
        const matrix = Mat4.mul(Mat4.zero(), second.matrix, first.matrix);
        return create(second.name, matrix, second);
    }

    export interface CoordinateMapper<T extends number> { (index: T, slot: Vec3): Vec3 }
    export interface ArrayMapping<T extends number> {
        readonly operator: SymmetryOperator,
        readonly invariantPosition: CoordinateMapper<T>,
        readonly position: CoordinateMapper<T>,
        x(index: T): number,
        y(index: T): number,
        z(index: T): number,
        r(index: T): number
    }

    export interface Coordinates { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> }

    export function createMapping<T extends number>(operator: SymmetryOperator, coords: Coordinates, radius: ((index: T) => number) | undefined): ArrayMapping<T> {
        const invariantPosition = SymmetryOperator.createCoordinateMapper(SymmetryOperator.Default, coords);
        const position = operator.isIdentity ? invariantPosition : SymmetryOperator.createCoordinateMapper(operator, coords);
        const { x, y, z } = createProjections(operator, coords);
        return { operator, invariantPosition, position, x, y, z, r: radius ? radius : _zeroRadius };
    }

    export function createCoordinateMapper<T extends number>(t: SymmetryOperator, coords: Coordinates): CoordinateMapper<T> {
        if (t.isIdentity) return identityPosition(coords);
        return generalPosition(t, coords);
    }
}

export { SymmetryOperator };

function _zeroRadius(i: number) { return 0; }

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

function projectX({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs }: SymmetryOperator.Coordinates) {
    const xx = m[0], yy = m[4], zz = m[8], tx = m[12];

    if (isW1(m)) {
        // this should always be the case.
        return (i: number) => xx * xs[i] + yy * ys[i] + zz * zs[i] + tx;
    }

    return (i: number) => {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + tx) / w;
    };
}

function projectY({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs }: SymmetryOperator.Coordinates) {
    const xx = m[1], yy = m[5], zz = m[9], ty = m[13];

    if (isW1(m)) {
        // this should always be the case.
        return (i: number) => xx * xs[i] + yy * ys[i] + zz * zs[i] + ty;
    }

    return (i: number) => {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + ty) / w;
    };
}

function projectZ({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs }: SymmetryOperator.Coordinates) {
    const xx = m[2], yy = m[6], zz = m[10], tz = m[14];

    if (isW1(m)) {
        // this should always be the case.
        return (i: number) => xx * xs[i] + yy * ys[i] + zz * zs[i] + tz;
    }

    return (i: number) => {
        const x = xs[i], y = ys[i], z = zs[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        return (xx * x + yy * y + zz * z + tz) / w;
    };
}

function identityPosition<T extends number>({ x, y, z }: SymmetryOperator.Coordinates): SymmetryOperator.CoordinateMapper<T> {
    return (i, s) => {
        s[0] = x[i];
        s[1] = y[i];
        s[2] = z[i];
        return s;
    };
}

function generalPosition<T extends number>({ matrix: m }: SymmetryOperator, { x: xs, y: ys, z: zs }: SymmetryOperator.Coordinates) {
    if (isW1(m)) {
        // this should always be the case.
        return (i: T, r: Vec3): Vec3 => {
            const x = xs[i], y = ys[i], z = zs[i];
            r[0] = m[0] * x + m[4] * y + m[8] * z + m[12];
            r[1] = m[1] * x + m[5] * y + m[9] * z + m[13];
            r[2] = m[2] * x + m[6] * y + m[10] * z + m[14];
            return r;
        };
    }
    return (i: T, r: Vec3): Vec3 => {
        r[0] = xs[i];
        r[1] = ys[i];
        r[2] = zs[i];
        Vec3.transformMat4(r, r, m);
        return r;
    };
}