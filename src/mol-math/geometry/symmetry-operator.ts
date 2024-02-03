/**
 * Copyright (c) 2017-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { lerp as scalar_lerp } from '../../mol-math/interpolate';
import { Mat3 } from '../linear-algebra/3d/mat3';
import { Mat4 } from '../linear-algebra/3d/mat4';
import { Quat } from '../linear-algebra/3d/quat';
import { Vec3 } from '../linear-algebra/3d/vec3';

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

    /** pointer to `struct_ncs_oper.id` or -1 */
    readonly ncsId: number,

    readonly hkl: Vec3,
    /** spacegroup symmetry operator index, -1 if not applicable */
    readonly spgrOp: number,
    /** unique (external) key, -1 if not available */
    readonly key: number,

    readonly matrix: Mat4,
    /** cache the inverse of the transform */
    readonly inverse: Mat4,
    /** optimize the identity case */
    readonly isIdentity: boolean,

    /**
     * Suffix based on operator type.
     * - Assembly: _assembly.operId
     * - Crystal: -op_ijk
     * - ncs: _ncsId
     */
    readonly suffix: string
}

namespace SymmetryOperator {
    export const DefaultName = '1_555';
    export const Default: SymmetryOperator = create(DefaultName, Mat4.identity());

    export const RotationTranslationEpsilon = 0.005;

    export type CreateInfo = { assembly?: SymmetryOperator['assembly'], ncsId?: number, hkl?: Vec3, spgrOp?: number, key?: number }
    export function create(name: string, matrix: Mat4, info?: CreateInfo | SymmetryOperator): SymmetryOperator {
        let { assembly, ncsId, hkl, spgrOp, key } = info || { };
        const _hkl = hkl ? Vec3.clone(hkl) : Vec3();
        spgrOp = spgrOp ?? -1;
        key = key ?? -1;
        ncsId = ncsId || -1;
        const isIdentity = Mat4.isIdentity(matrix);
        const suffix = getSuffix(info, isIdentity);
        if (isIdentity) return { name, assembly, matrix, inverse: Mat4.identity(), isIdentity: true, hkl: _hkl, spgrOp, ncsId, suffix, key };
        if (!Mat4.isRotationAndTranslation(matrix, RotationTranslationEpsilon)) {
            console.warn(`Symmetry operator (${name}) should be a composition of rotation and translation.`);
        }
        return { name, assembly, matrix, inverse: Mat4.invert(Mat4(), matrix), isIdentity: false, hkl: _hkl, spgrOp, key, ncsId, suffix };
    }

    function isSymmetryOperator(x: any): x is SymmetryOperator {
        return !!x && !!x.matrix && !!x.inverse && typeof x.name === 'string';
    }

    function getSuffix(info?: CreateInfo | SymmetryOperator, isIdentity?: boolean) {
        if (!info) return '';

        if (info.assembly) {
            if (isSymmetryOperator(info)) return info.suffix;
            return isIdentity ? '' : `_${info.assembly.operId}`;
        }

        if (typeof info.spgrOp !== 'undefined' && typeof info.hkl !== 'undefined' && info.spgrOp !== -1) {
            const [i, j, k] = info.hkl;
            return `-${info.spgrOp + 1}_${5 + i}${5 + j}${5 + k}`;
        }

        if (info.ncsId !== -1) {
            return `_${info.ncsId}`;
        }

        return '';
    }

    const _m = Mat4();
    export function checkIfRotationAndTranslation(rot: Mat3, offset: Vec3) {
        Mat4.setIdentity(_m);
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                Mat4.setValue(_m, i, j, Mat3.getValue(rot, i, j));
            }
        }
        Mat4.setTranslation(_m, offset);
        return Mat4.isRotationAndTranslation(_m, RotationTranslationEpsilon);
    }

    export function ofRotationAndOffset(name: string, rot: Mat3, offset: Vec3, ncsId?: number) {
        const t = Mat4.identity();
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                Mat4.setValue(t, i, j, Mat3.getValue(rot, i, j));
            }
        }
        Mat4.setTranslation(t, offset);
        return create(name, t, { ncsId });
    }

    const _q1 = Quat.identity(), _q2 = Quat(), _q3 = Quat(), _axis = Vec3();
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
        const matrix = Mat4.mul(Mat4(), second.matrix, first.matrix);
        return create(second.name, matrix, second);
    }

    export interface CoordinateMapper<T extends number> { (index: T, slot: Vec3): Vec3 }
    export interface ArrayMapping<T extends number> {
        readonly coordinates: Coordinates,
        readonly operator: SymmetryOperator,
        invariantPosition(this: ArrayMapping<T>, i: T, s: Vec3): Vec3,
        position(this: ArrayMapping<T>, i: T, s: Vec3): Vec3,
        x(this: ArrayMapping<T>, index: T): number,
        y(this: ArrayMapping<T>, index: T): number,
        z(this: ArrayMapping<T>, index: T): number,
        readonly r: (this: ArrayMapping<T>, index: T) => number
    }

    class _ArrayMapping<T extends number> implements ArrayMapping<T> {
        private readonly _x: ArrayLike<number>;
        private readonly _y: ArrayLike<number>;
        private readonly _z: ArrayLike<number>;
        private readonly _m: Mat4;

        constructor(readonly operator: SymmetryOperator, readonly coordinates: Coordinates, readonly r: ((index: T) => number) = _zeroRadius) {
            this._x = coordinates.x;
            this._y = coordinates.y;
            this._z = coordinates.z;
            this._m = operator.matrix;
        }

        invariantPosition(i: T, s: Vec3): Vec3 {
            s[0] = this._x[i];
            s[1] = this._y[i];
            s[2] = this._z[i];
            return s;
        }

        position(i: T, s: Vec3): Vec3 {
            s[0] = this._x[i];
            s[1] = this._y[i];
            s[2] = this._z[i];
            Vec3.transformMat4(s, s, this._m);
            return s;
        }

        x(i: T): number {
            const m = this._m;
            const xx = m[0], yy = m[4], zz = m[8], tx = m[12];

            const x = this._x[i], y = this._y[i], z = this._z[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
            return (xx * x + yy * y + zz * z + tx) / w;
        }

        y(i: T): number {
            const m = this._m;
            const xx = m[1], yy = m[5], zz = m[9], ty = m[13];

            const x = this._x[i], y = this._y[i], z = this._z[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
            return (xx * x + yy * y + zz * z + ty) / w;
        }

        z(i: T): number {
            const m = this._m;
            const xx = m[2], yy = m[6], zz = m[10], tz = m[14];

            const x = this._x[i], y = this._y[i], z = this._z[i], w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
            return (xx * x + yy * y + zz * z + tz) / w;
        }
    }

    class _ArrayMappingW1<T extends number> implements ArrayMapping<T> {
        private readonly _x: ArrayLike<number>;
        private readonly _y: ArrayLike<number>;
        private readonly _z: ArrayLike<number>;
        private readonly _m: Mat4;

        constructor(readonly operator: SymmetryOperator, readonly coordinates: Coordinates, readonly r: ((index: T) => number) = _zeroRadius) {
            this._x = coordinates.x;
            this._y = coordinates.y;
            this._z = coordinates.z;
            this._m = operator.matrix;
        }

        invariantPosition(i: T, s: Vec3): Vec3 {
            s[0] = this._x[i];
            s[1] = this._y[i];
            s[2] = this._z[i];
            return s;
        }

        position(i: T, s: Vec3): Vec3 {
            s[0] = this.x(i);
            s[1] = this.y(i);
            s[2] = this.z(i);
            return s;
        }

        x(i: T): number {
            const m = this._m;
            return m[0] * this._x[i] + m[4] * this._y[i] + m[8] * this._z[i] + m[12];
        }

        y(i: T): number {
            const m = this._m;
            return m[1] * this._x[i] + m[5] * this._y[i] + m[9] * this._z[i] + m[13];
        }

        z(i: T): number {
            const m = this._m;
            return m[2] * this._x[i] + m[6] * this._y[i] + m[10] * this._z[i] + m[14];
        }
    }

    class _ArrayMappingIdentity<T extends number> implements ArrayMapping<T> {
        private readonly _x: ArrayLike<number>;
        private readonly _y: ArrayLike<number>;
        private readonly _z: ArrayLike<number>;

        constructor(readonly operator: SymmetryOperator, readonly coordinates: Coordinates, readonly r: ((index: T) => number) = _zeroRadius) {
            this._x = coordinates.x;
            this._y = coordinates.y;
            this._z = coordinates.z;
        }

        invariantPosition(i: T, s: Vec3): Vec3 {
            s[0] = this._x[i];
            s[1] = this._y[i];
            s[2] = this._z[i];
            return s;
        }

        position(i: T, s: Vec3): Vec3 {
            s[0] = this._x[i];
            s[1] = this._y[i];
            s[2] = this._z[i];
            return s;
        }

        x(i: T): number {
            return this._x[i];
        }

        y(i: T): number {
            return this._y[i];
        }

        z(i: T): number {
            return this._z[i];
        }
    }

    export interface Coordinates { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> }

    export function createMapping<T extends number>(operator: SymmetryOperator, coords: Coordinates, radius: ((index: T) => number) = _zeroRadius): ArrayMapping<T> {
        return Mat4.isIdentity(operator.matrix)
            ? new _ArrayMappingIdentity(operator, coords, radius)
            : isW1(operator.matrix)
                ? new _ArrayMappingW1(operator, coords, radius)
                : new _ArrayMapping(operator, coords, radius);
    }
}

export { SymmetryOperator };

function _zeroRadius(_i: number) { return 0; }

function isW1(m: Mat4) {
    return m[3] === 0 && m[7] === 0 && m[11] === 0 && m[15] === 1;
}
