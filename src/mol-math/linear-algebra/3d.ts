/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/*
 * This code has been modified from https://github.com/toji/gl-matrix/,
 * copyright (c) 2015, Brandon Jones, Colin MacKenzie IV.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 */

export interface Mat4 extends Array<number> { [d: number]: number, '@type': 'mat4', length: 16 }
export interface Mat3 extends Array<number> { [d: number]: number, '@type': 'mat3', length: 9 }
export interface Vec3 extends Array<number> { [d: number]: number, '@type': 'vec3', length: 3 }
export interface Vec4 extends Array<number> { [d: number]: number, '@type': 'vec4', length: 4 }
export interface Quat extends Array<number> { [d: number]: number, '@type': 'quat', length: 4 }

const enum EPSILON { Value = 0.000001 }

export function Mat4() {
    return Mat4.zero();
}

export function Quat() {
    return Quat.zero();
}

/**
 * Stores a 4x4 matrix in a column major (j * 4 + i indexing) format.
 */
export namespace Mat4 {
    export function zero(): Mat4 {
        // force double backing array by 0.1.
        const ret = [0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        ret[0] = 0.0;
        return ret as any;
    }

    export function identity(): Mat4 {
        const out = zero();
        out[0] = 1;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = 1;
        out[6] = 0;
        out[7] = 0;
        out[8] = 0;
        out[9] = 0;
        out[10] = 1;
        out[11] = 0;
        out[12] = 0;
        out[13] = 0;
        out[14] = 0;
        out[15] = 1;
        return out;
    }

    export function setIdentity(mat: Mat4): Mat4 {
        mat[0] = 1;
        mat[1] = 0;
        mat[2] = 0;
        mat[3] = 0;
        mat[4] = 0;
        mat[5] = 1;
        mat[6] = 0;
        mat[7] = 0;
        mat[8] = 0;
        mat[9] = 0;
        mat[10] = 1;
        mat[11] = 0;
        mat[12] = 0;
        mat[13] = 0;
        mat[14] = 0;
        mat[15] = 1;
        return mat;
    }

    export function ofRows(rows: number[][]): Mat4 {
        const out = zero();
        for (let i = 0; i < 4; i++) {
            const r = rows[i];
            for (let j = 0; j < 4; j++) {
                out[4 * j + i] = r[j];
            }
        }
        return out;
    }

    const _id = identity();
    export function isIdentity(m: Mat4, eps?: number) {
        return areEqual(m, _id, typeof eps === 'undefined' ? EPSILON.Value : eps);
    }

    export function areEqual(a: Mat4, b: Mat4, eps: number) {
        for (let i = 0; i < 16; i++) {
            if (Math.abs(a[i] - b[i]) > eps) return false;
        }
        return true;
    }

    export function setValue(a: Mat4, i: number, j: number, value: number) {
        a[4 * j + i] = value;
    }

    export function toArray(a: Mat4, out: Helpers.NumberArray, offset: number) {
        out[offset + 0] = a[0];
        out[offset + 1] = a[1];
        out[offset + 2] = a[2];
        out[offset + 3] = a[3];
        out[offset + 4] = a[4];
        out[offset + 5] = a[5];
        out[offset + 6] = a[6];
        out[offset + 7] = a[7];
        out[offset + 8] = a[8];
        out[offset + 9] = a[9];
        out[offset + 10] = a[10];
        out[offset + 11] = a[11];
        out[offset + 12] = a[12];
        out[offset + 13] = a[13];
        out[offset + 14] = a[14];
        out[offset + 15] = a[15];
    }

    export function fromArray(a: Mat4, array: Helpers.NumberArray, offset: number) {
        a[0] = array[offset + 0]
        a[1] = array[offset + 1]
        a[2] = array[offset + 2]
        a[3] = array[offset + 3]
        a[4] = array[offset + 4]
        a[5] = array[offset + 5]
        a[6] = array[offset + 6]
        a[7] = array[offset + 7]
        a[8] = array[offset + 8]
        a[9] = array[offset + 9]
        a[10] = array[offset + 10]
        a[11] = array[offset + 11]
        a[12] = array[offset + 12]
        a[13] = array[offset + 13]
        a[14] = array[offset + 14]
        a[15] = array[offset + 15]
    }

    export function copy(out: Mat4, a: Mat4) {
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        out[3] = a[3];
        out[4] = a[4];
        out[5] = a[5];
        out[6] = a[6];
        out[7] = a[7];
        out[8] = a[8];
        out[9] = a[9];
        out[10] = a[10];
        out[11] = a[11];
        out[12] = a[12];
        out[13] = a[13];
        out[14] = a[14];
        out[15] = a[15];
        return out;
    }

    export function clone(a: Mat4) {
        return Mat4.copy(Mat4.zero(), a);
    }

    export function transpose(out: Mat4, a: Mat4) {
        // If we are transposing ourselves we can skip a few steps but have to cache some values
        if (out === a) {
            const a01 = a[1], a02 = a[2], a03 = a[3];
            const a12 = a[6], a13 = a[7];
            const a23 = a[11];
            out[1] = a[4];
            out[2] = a[8];
            out[3] = a[12];
            out[4] = a01;
            out[6] = a[9];
            out[7] = a[13];
            out[8] = a02;
            out[9] = a12;
            out[11] = a[14];
            out[12] = a03;
            out[13] = a13;
            out[14] = a23;
        } else {
            out[0] = a[0];
            out[1] = a[4];
            out[2] = a[8];
            out[3] = a[12];
            out[4] = a[1];
            out[5] = a[5];
            out[6] = a[9];
            out[7] = a[13];
            out[8] = a[2];
            out[9] = a[6];
            out[10] = a[10];
            out[11] = a[14];
            out[12] = a[3];
            out[13] = a[7];
            out[14] = a[11];
            out[15] = a[15];
        }
        return out;
    }

    export function invert(out: Mat4, a: Mat4) {
        const a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
            a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7],
            a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11],
            a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15],

            b00 = a00 * a11 - a01 * a10,
            b01 = a00 * a12 - a02 * a10,
            b02 = a00 * a13 - a03 * a10,
            b03 = a01 * a12 - a02 * a11,
            b04 = a01 * a13 - a03 * a11,
            b05 = a02 * a13 - a03 * a12,
            b06 = a20 * a31 - a21 * a30,
            b07 = a20 * a32 - a22 * a30,
            b08 = a20 * a33 - a23 * a30,
            b09 = a21 * a32 - a22 * a31,
            b10 = a21 * a33 - a23 * a31,
            b11 = a22 * a33 - a23 * a32;

        // Calculate the determinant
        let det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

        if (!det) {
            console.warn('non-invertible matrix.', a);
            return out;
        }
        det = 1.0 / det;

        out[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
        out[1] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
        out[2] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
        out[3] = (a22 * b04 - a21 * b05 - a23 * b03) * det;
        out[4] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
        out[5] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
        out[6] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
        out[7] = (a20 * b05 - a22 * b02 + a23 * b01) * det;
        out[8] = (a10 * b10 - a11 * b08 + a13 * b06) * det;
        out[9] = (a01 * b08 - a00 * b10 - a03 * b06) * det;
        out[10] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
        out[11] = (a21 * b02 - a20 * b04 - a23 * b00) * det;
        out[12] = (a11 * b07 - a10 * b09 - a12 * b06) * det;
        out[13] = (a00 * b09 - a01 * b07 + a02 * b06) * det;
        out[14] = (a31 * b01 - a30 * b03 - a32 * b00) * det;
        out[15] = (a20 * b03 - a21 * b01 + a22 * b00) * det;

        return out;
    }

    export function mul(out: Mat4, a: Mat4, b: Mat4) {
        const a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
            a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7],
            a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11],
            a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15];

        // Cache only the current line of the second matrix
        let b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        out[0] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
        out[1] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
        out[2] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
        out[3] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

        b0 = b[4]; b1 = b[5]; b2 = b[6]; b3 = b[7];
        out[4] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
        out[5] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
        out[6] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
        out[7] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

        b0 = b[8]; b1 = b[9]; b2 = b[10]; b3 = b[11];
        out[8] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
        out[9] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
        out[10] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
        out[11] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

        b0 = b[12]; b1 = b[13]; b2 = b[14]; b3 = b[15];
        out[12] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
        out[13] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
        out[14] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
        out[15] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;
        return out;
    }

    export function mul3(out: Mat4, a: Mat4, b: Mat4, c: Mat4) {
        return mul(out, mul(out, a, b), c);
    }

    export function translate(out: Mat4, a: Mat4, v: Vec3) {
        const x = v[0], y = v[1], z = v[2];
        let a00: number, a01: number, a02: number, a03: number,
            a10: number, a11: number, a12: number, a13: number,
            a20: number, a21: number, a22: number, a23: number;

        if (a === out) {
            out[12] = a[0] * x + a[4] * y + a[8] * z + a[12];
            out[13] = a[1] * x + a[5] * y + a[9] * z + a[13];
            out[14] = a[2] * x + a[6] * y + a[10] * z + a[14];
            out[15] = a[3] * x + a[7] * y + a[11] * z + a[15];
        } else {
            a00 = a[0]; a01 = a[1]; a02 = a[2]; a03 = a[3];
            a10 = a[4]; a11 = a[5]; a12 = a[6]; a13 = a[7];
            a20 = a[8]; a21 = a[9]; a22 = a[10]; a23 = a[11];

            out[0] = a00; out[1] = a01; out[2] = a02; out[3] = a03;
            out[4] = a10; out[5] = a11; out[6] = a12; out[7] = a13;
            out[8] = a20; out[9] = a21; out[10] = a22; out[11] = a23;

            out[12] = a00 * x + a10 * y + a20 * z + a[12];
            out[13] = a01 * x + a11 * y + a21 * z + a[13];
            out[14] = a02 * x + a12 * y + a22 * z + a[14];
            out[15] = a03 * x + a13 * y + a23 * z + a[15];
        }

        return out;
    }

    export function fromTranslation(out: Mat4, v: Vec3) {
        out[0] = 1;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = 1;
        out[6] = 0;
        out[7] = 0;
        out[8] = 0;
        out[9] = 0;
        out[10] = 1;
        out[11] = 0;
        out[12] = v[0];
        out[13] = v[1];
        out[14] = v[2];
        out[15] = 1;
        return out;
    }

    export function setTranslation(out: Mat4, v: Vec3) {
        out[12] = v[0];
        out[13] = v[1];
        out[14] = v[2];
        return out;
    }

    export function rotate(out: Mat4, a: Mat4, rad: number, axis: Mat4) {
        let x = axis[0], y = axis[1], z = axis[2],
            len = Math.sqrt(x * x + y * y + z * z),
            s, c, t,
            a00, a01, a02, a03,
            a10, a11, a12, a13,
            a20, a21, a22, a23,
            b00, b01, b02,
            b10, b11, b12,
            b20, b21, b22;

        if (Math.abs(len) < EPSILON.Value) {
            return Mat4.identity();
        }

        len = 1 / len;
        x *= len;
        y *= len;
        z *= len;

        s = Math.sin(rad);
        c = Math.cos(rad);
        t = 1 - c;

        a00 = a[0]; a01 = a[1]; a02 = a[2]; a03 = a[3];
        a10 = a[4]; a11 = a[5]; a12 = a[6]; a13 = a[7];
        a20 = a[8]; a21 = a[9]; a22 = a[10]; a23 = a[11];

        // Construct the elements of the rotation matrix
        b00 = x * x * t + c; b01 = y * x * t + z * s; b02 = z * x * t - y * s;
        b10 = x * y * t - z * s; b11 = y * y * t + c; b12 = z * y * t + x * s;
        b20 = x * z * t + y * s; b21 = y * z * t - x * s; b22 = z * z * t + c;

        // Perform rotation-specific matrix multiplication
        out[0] = a00 * b00 + a10 * b01 + a20 * b02;
        out[1] = a01 * b00 + a11 * b01 + a21 * b02;
        out[2] = a02 * b00 + a12 * b01 + a22 * b02;
        out[3] = a03 * b00 + a13 * b01 + a23 * b02;
        out[4] = a00 * b10 + a10 * b11 + a20 * b12;
        out[5] = a01 * b10 + a11 * b11 + a21 * b12;
        out[6] = a02 * b10 + a12 * b11 + a22 * b12;
        out[7] = a03 * b10 + a13 * b11 + a23 * b12;
        out[8] = a00 * b20 + a10 * b21 + a20 * b22;
        out[9] = a01 * b20 + a11 * b21 + a21 * b22;
        out[10] = a02 * b20 + a12 * b21 + a22 * b22;
        out[11] = a03 * b20 + a13 * b21 + a23 * b22;

        if (a !== out) { // If the source and destination differ, copy the unchanged last row
            out[12] = a[12];
            out[13] = a[13];
            out[14] = a[14];
            out[15] = a[15];
        }
        return out;
    }

    export function fromRotation(out: Mat4, rad: number, axis: Vec3) {
        let x = axis[0], y = axis[1], z = axis[2],
            len = Math.sqrt(x * x + y * y + z * z),
            s, c, t;

        if (Math.abs(len) < EPSILON.Value) { return setIdentity(out); }

        len = 1 / len;
        x *= len;
        y *= len;
        z *= len;

        s = Math.sin(rad);
        c = Math.cos(rad);
        t = 1 - c;

        // Perform rotation-specific matrix multiplication
        out[0] = x * x * t + c;
        out[1] = y * x * t + z * s;
        out[2] = z * x * t - y * s;
        out[3] = 0;
        out[4] = x * y * t - z * s;
        out[5] = y * y * t + c;
        out[6] = z * y * t + x * s;
        out[7] = 0;
        out[8] = x * z * t + y * s;
        out[9] = y * z * t - x * s;
        out[10] = z * z * t + c;
        out[11] = 0;
        out[12] = 0;
        out[13] = 0;
        out[14] = 0;
        out[15] = 1;
        return out;
    }

    export function scale(out: Mat4, a: Mat4, v: Vec3) {
        const x = v[0], y = v[1], z = v[2];

        out[0] = a[0] * x;
        out[1] = a[1] * x;
        out[2] = a[2] * x;
        out[3] = a[3] * x;
        out[4] = a[4] * y;
        out[5] = a[5] * y;
        out[6] = a[6] * y;
        out[7] = a[7] * y;
        out[8] = a[8] * z;
        out[9] = a[9] * z;
        out[10] = a[10] * z;
        out[11] = a[11] * z;
        out[12] = a[12];
        out[13] = a[13];
        out[14] = a[14];
        out[15] = a[15];
        return out;
    }

    export function fromScaling(out: Mat4, v: Vec3) {
        out[0] = v[0];
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = v[1];
        out[6] = 0;
        out[7] = 0;
        out[8] = 0;
        out[9] = 0;
        out[10] = v[2];
        out[11] = 0;
        out[12] = 0;
        out[13] = 0;
        out[14] = 0;
        out[15] = 1;
        return out;
    }

    export function makeTable(m: Mat4) {
        let ret = '';
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                ret += m[4 * j + i].toString();
                if (j < 3) ret += ' ';
            }
            if (i < 3) ret += '\n';
        }
        return ret;
    }

    export function determinant(a: Mat4) {
        const a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
            a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7],
            a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11],
            a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15],

            b00 = a00 * a11 - a01 * a10,
            b01 = a00 * a12 - a02 * a10,
            b02 = a00 * a13 - a03 * a10,
            b03 = a01 * a12 - a02 * a11,
            b04 = a01 * a13 - a03 * a11,
            b05 = a02 * a13 - a03 * a12,
            b06 = a20 * a31 - a21 * a30,
            b07 = a20 * a32 - a22 * a30,
            b08 = a20 * a33 - a23 * a30,
            b09 = a21 * a32 - a22 * a31,
            b10 = a21 * a33 - a23 * a31,
            b11 = a22 * a33 - a23 * a32;

        // Calculate the determinant
        return b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
    }

    /**
     * Check if the matrix has the form
     * [ Rotation    Translation ]
     * [ 0           1           ]
     */
    export function isRotationAndTranslation(a: Mat4, eps?: number) {
        return _isRotationAndTranslation(a, typeof eps !== 'undefined' ? eps : EPSILON.Value)
    }

    function _isRotationAndTranslation(a: Mat4, eps: number) {
        const a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
            a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7],
            a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11],
            /* a30 = a[12], a31 = a[13], a32 = a[14],*/ a33 = a[15];

        if (a33 !== 1 || a03 !== 0 || a13 !== 0 || a23 !== 0) {
            return false;
        }
        const det3x3 = a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20) + a02 * (a10 * a21 - a11 * a20);
        if (det3x3 < 1 - eps || det3x3 > 1 + eps) {
            return false;
        }
        return true;
    }

    export function fromQuat(out: Mat4, q: Quat) {
        const x = q[0], y = q[1], z = q[2], w = q[3];
        const x2 = x + x;
        const y2 = y + y;
        const z2 = z + z;

        const xx = x * x2;
        const yx = y * x2;
        const yy = y * y2;
        const zx = z * x2;
        const zy = z * y2;
        const zz = z * z2;
        const wx = w * x2;
        const wy = w * y2;
        const wz = w * z2;

        out[0] = 1 - yy - zz;
        out[1] = yx + wz;
        out[2] = zx - wy;
        out[3] = 0;

        out[4] = yx - wz;
        out[5] = 1 - xx - zz;
        out[6] = zy + wx;
        out[7] = 0;

        out[8] = zx + wy;
        out[9] = zy - wx;
        out[10] = 1 - xx - yy;
        out[11] = 0;

        out[12] = 0;
        out[13] = 0;
        out[14] = 0;
        out[15] = 1;

        return out;
    }

    /**
     * Generates a frustum matrix with the given bounds
     */
    export function frustum(out: Mat4, left: number, right: number, bottom: number, top: number, near: number, far: number) {
        let rl = 1 / (right - left);
        let tb = 1 / (top - bottom);
        let nf = 1 / (near - far);
        out[0] = (near * 2) * rl;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = (near * 2) * tb;
        out[6] = 0;
        out[7] = 0;
        out[8] = (right + left) * rl;
        out[9] = (top + bottom) * tb;
        out[10] = (far + near) * nf;
        out[11] = -1;
        out[12] = 0;
        out[13] = 0;
        out[14] = (far * near * 2) * nf;
        out[15] = 0;
        return out;
    }

    /**
     * Generates a perspective projection matrix with the given bounds
     */
    export function perspective(out: Mat4, fovy: number, aspect: number, near: number, far: number) {
        let f = 1.0 / Math.tan(fovy / 2);
        let nf = 1 / (near - far);
        out[0] = f / aspect;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = f;
        out[6] = 0;
        out[7] = 0;
        out[8] = 0;
        out[9] = 0;
        out[10] = (far + near) * nf;
        out[11] = -1;
        out[12] = 0;
        out[13] = 0;
        out[14] = (2 * far * near) * nf;
        out[15] = 0;
        return out;
    }

    /**
     * Generates a orthogonal projection matrix with the given bounds
     */
    export function ortho(out: Mat4, left: number, right: number, bottom: number, top: number, near: number, far: number) {
        let lr = 1 / (left - right);
        let bt = 1 / (bottom - top);
        let nf = 1 / (near - far);
        out[0] = -2 * lr;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = -2 * bt;
        out[6] = 0;
        out[7] = 0;
        out[8] = 0;
        out[9] = 0;
        out[10] = 2 * nf;
        out[11] = 0;
        out[12] = (left + right) * lr;
        out[13] = (top + bottom) * bt;
        out[14] = (far + near) * nf;
        out[15] = 1;
        return out;
    }

    /**
     * Generates a look-at matrix with the given eye position, focal point, and up axis
     */
    export function lookAt(out: Mat4, eye: Vec3, center: Vec3, up: Vec3) {
        let x0, x1, x2, y0, y1, y2, z0, z1, z2, len;
        let eyex = eye[0];
        let eyey = eye[1];
        let eyez = eye[2];
        let upx = up[0];
        let upy = up[1];
        let upz = up[2];
        let centerx = center[0];
        let centery = center[1];
        let centerz = center[2];

        if (Math.abs(eyex - centerx) < EPSILON.Value &&
            Math.abs(eyey - centery) < EPSILON.Value &&
            Math.abs(eyez - centerz) < EPSILON.Value
        ) {
            return setIdentity(out);
        }

        z0 = eyex - centerx;
        z1 = eyey - centery;
        z2 = eyez - centerz;

        len = 1 / Math.sqrt(z0 * z0 + z1 * z1 + z2 * z2);
        z0 *= len;
        z1 *= len;
        z2 *= len;

        x0 = upy * z2 - upz * z1;
        x1 = upz * z0 - upx * z2;
        x2 = upx * z1 - upy * z0;
        len = Math.sqrt(x0 * x0 + x1 * x1 + x2 * x2);
        if (!len) {
            x0 = 0;
            x1 = 0;
            x2 = 0;
        } else {
            len = 1 / len;
            x0 *= len;
            x1 *= len;
            x2 *= len;
        }

        y0 = z1 * x2 - z2 * x1;
        y1 = z2 * x0 - z0 * x2;
        y2 = z0 * x1 - z1 * x0;

        len = Math.sqrt(y0 * y0 + y1 * y1 + y2 * y2);
        if (!len) {
            y0 = 0;
            y1 = 0;
            y2 = 0;
        } else {
            len = 1 / len;
            y0 *= len;
            y1 *= len;
            y2 *= len;
        }

        out[0] = x0;
        out[1] = y0;
        out[2] = z0;
        out[3] = 0;
        out[4] = x1;
        out[5] = y1;
        out[6] = z1;
        out[7] = 0;
        out[8] = x2;
        out[9] = y2;
        out[10] = z2;
        out[11] = 0;
        out[12] = -(x0 * eyex + x1 * eyey + x2 * eyez);
        out[13] = -(y0 * eyex + y1 * eyey + y2 * eyez);
        out[14] = -(z0 * eyex + z1 * eyey + z2 * eyez);
        out[15] = 1;

        return out;
    }
}

export namespace Mat3 {
    export function zero(): Mat3 {
        // force double backing array by 0.1.
        const ret = [0.1, 0, 0, 0, 0, 0, 0, 0, 0];
        ret[0] = 0.0;
        return ret as any;
    }
}

export namespace Vec3 {
    export function zero(): Vec3 {
        const out = [0.1, 0.0, 0.0];
        out[0] = 0;
        return out as any;
    }

    export function clone(a: Vec3): Vec3 {
        const out = zero();
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        return out;
    }

    export function fromObj(v: { x: number, y: number, z: number }): Vec3 {
        return create(v.x, v.y, v.z);
    }

    export function toObj(v: Vec3) {
        return { x: v[0], y: v[1], z: v[2] };
    }

    export function fromArray(v: Vec3, array: Helpers.NumberArray, offset: number) {
        v[0] = array[offset + 0]
        v[1] = array[offset + 1]
        v[2] = array[offset + 2]
    }

    export function toArray(v: Vec3, out: Helpers.NumberArray, offset: number) {
        out[offset + 0] = v[0]
        out[offset + 1] = v[1]
        out[offset + 2] = v[2]
    }

    export function create(x: number, y: number, z: number): Vec3 {
        const out = zero();
        out[0] = x;
        out[1] = y;
        out[2] = z;
        return out;
    }

    export function set(out: Vec3, x: number, y: number, z: number): Vec3 {
        out[0] = x;
        out[1] = y;
        out[2] = z;
        return out;
    }

    export function copy(out: Vec3, a: Vec3) {
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        return out;
    }

    export function add(out: Vec3, a: Vec3, b: Vec3) {
        out[0] = a[0] + b[0];
        out[1] = a[1] + b[1];
        out[2] = a[2] + b[2];
        return out;
    }

    export function sub(out: Vec3, a: Vec3, b: Vec3) {
        out[0] = a[0] - b[0];
        out[1] = a[1] - b[1];
        out[2] = a[2] - b[2];
        return out;
    }

    export function scale(out: Vec3, a: Vec3, b: number) {
        out[0] = a[0] * b;
        out[1] = a[1] * b;
        out[2] = a[2] * b;
        return out;
    }

    export function scaleAndAdd(out: Vec3, a: Vec3, b: Vec3, scale: number) {
        out[0] = a[0] + (b[0] * scale);
        out[1] = a[1] + (b[1] * scale);
        out[2] = a[2] + (b[2] * scale);
        return out;
    }

    export function distance(a: Vec3, b: Vec3) {
        const x = b[0] - a[0],
            y = b[1] - a[1],
            z = b[2] - a[2];
        return Math.sqrt(x * x + y * y + z * z);
    }

    export function squaredDistance(a: Vec3, b: Vec3) {
        const x = b[0] - a[0],
            y = b[1] - a[1],
            z = b[2] - a[2];
        return x * x + y * y + z * z;
    }

    export function magnitude(a: Vec3) {
        const x = a[0],
            y = a[1],
            z = a[2];
        return Math.sqrt(x * x + y * y + z * z);
    }

    export function squaredMagnitude(a: Vec3) {
        const x = a[0],
            y = a[1],
            z = a[2];
        return x * x + y * y + z * z;
    }

    export function normalize(out: Vec3, a: Vec3) {
        const x = a[0],
            y = a[1],
            z = a[2];
        let len = x * x + y * y + z * z;
        if (len > 0) {
            len = 1 / Math.sqrt(len);
            out[0] = a[0] * len;
            out[1] = a[1] * len;
            out[2] = a[2] * len;
        }
        return out;
    }

    export function dot(a: Vec3, b: Vec3) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    export function cross(out: Vec3, a: Vec3, b: Vec3) {
        const ax = a[0], ay = a[1], az = a[2],
            bx = b[0], by = b[1], bz = b[2];

        out[0] = ay * bz - az * by;
        out[1] = az * bx - ax * bz;
        out[2] = ax * by - ay * bx;
        return out;
    }

    export function lerp(out: Vec3, a: Vec3, b: Vec3, t: number) {
        const ax = a[0],
            ay = a[1],
            az = a[2];
        out[0] = ax + t * (b[0] - ax);
        out[1] = ay + t * (b[1] - ay);
        out[2] = az + t * (b[2] - az);
        return out;
    }

    export function transformMat4(out: Vec3, a: Vec3, m: Mat4) {
        const x = a[0], y = a[1], z = a[2],
            w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        out[0] = (m[0] * x + m[4] * y + m[8] * z + m[12]) / w;
        out[1] = (m[1] * x + m[5] * y + m[9] * z + m[13]) / w;
        out[2] = (m[2] * x + m[6] * y + m[10] * z + m[14]) / w;
        return out;
    }

    const angleTempA = zero(), angleTempB = zero();
    export function angle(a: Vec3, b: Vec3) {
        copy(angleTempA, a);
        copy(angleTempB, b);

        normalize(angleTempA, angleTempA);
        normalize(angleTempB, angleTempB);

        const cosine = dot(angleTempA, angleTempB);

        if (cosine > 1.0) {
            return 0;
        }
        else if (cosine < -1.0) {
            return Math.PI;
        } else {
            return Math.acos(cosine);
        }
    }

    const rotTemp = zero();
    export function makeRotation(mat: Mat4, a: Vec3, b: Vec3): Mat4 {
        const by = angle(a, b);
        if (Math.abs(by) < 0.0001) return Mat4.setIdentity(mat);
        const axis = cross(rotTemp, a, b);
        return Mat4.fromRotation(mat, by, axis);
    }
}

export namespace Vec4 {
    export function zero(): Vec4 {
        // force double backing array by 0.1.
        const ret = [0.1, 0, 0, 0];
        ret[0] = 0.0;
        return ret as any;
    }

    export function clone(a: Vec4) {
        const out = zero();
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        out[3] = a[3];
        return out;
    }

    export function create(x: number, y: number, z: number, w: number) {
        const out = zero();
        out[0] = x;
        out[1] = y;
        out[2] = z;
        out[3] = w;
        return out;
    }

    export function copy(out: Vec4, a: Vec4) {
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        out[3] = a[3];
        return out;
    }

    export function set(out: Vec4, x: number, y: number, z: number, w: number) {
        out[0] = x;
        out[1] = y;
        out[2] = z;
        out[3] = w;
        return out;
    }

    export function add(out: Quat, a: Quat, b: Quat) {
        out[0] = a[0] + b[0];
        out[1] = a[1] + b[1];
        out[2] = a[2] + b[2];
        out[3] = a[3] + b[3];
        return out;
    }

    export function distance(a: Vec4, b: Vec4) {
        const x = b[0] - a[0],
            y = b[1] - a[1],
            z = b[2] - a[2],
            w = b[3] - a[3];
        return Math.sqrt(x * x + y * y + z * z + w * w);
    }

    export function squaredDistance(a: Vec4, b: Vec4) {
        const x = b[0] - a[0],
            y = b[1] - a[1],
            z = b[2] - a[2],
            w = b[3] - a[3];
        return x * x + y * y + z * z + w * w;
    }

    export function norm(a: Vec4) {
        const x = a[0],
            y = a[1],
            z = a[2],
            w = a[3];
        return Math.sqrt(x * x + y * y + z * z + w * w);
    }

    export function squaredNorm(a: Vec4) {
        const x = a[0],
            y = a[1],
            z = a[2],
            w = a[3];
        return x * x + y * y + z * z + w * w;
    }

    export function transform(out: Vec4, a: Vec4, m: Mat4) {
        const x = a[0], y = a[1], z = a[2], w = a[3];
        out[0] = m[0] * x + m[4] * y + m[8] * z + m[12] * w;
        out[1] = m[1] * x + m[5] * y + m[9] * z + m[13] * w;
        out[2] = m[2] * x + m[6] * y + m[10] * z + m[14] * w;
        out[3] = m[3] * x + m[7] * y + m[11] * z + m[15] * w;
        return out;
    }
}

export namespace Quat {
    export function zero(): Quat {
        // force double backing array by 0.1.
        const ret = [0.1, 0, 0, 0];
        ret[0] = 0.0;
        return ret as any;
    }

    export function identity(): Quat {
        const out = zero();
        out[3] = 1;
        return out;
    }

    export function create(x: number, y: number, z: number, w: number) {
        const out = identity();
        out[0] = x;
        out[1] = y;
        out[2] = z;
        out[3] = w;
        return out;
    }

    export function setAxisAngle(out: Quat, axis: Vec3, rad: number) {
        rad = rad * 0.5;
        let s = Math.sin(rad);
        out[0] = s * axis[0];
        out[1] = s * axis[1];
        out[2] = s * axis[2];
        out[3] = Math.cos(rad);
        return out;
    }

    /**
     * Gets the rotation axis and angle for a given
     *  quaternion. If a quaternion is created with
     *  setAxisAngle, this method will return the same
     *  values as providied in the original parameter list
     *  OR functionally equivalent values.
     * Example: The quaternion formed by axis [0, 0, 1] and
     *  angle -90 is the same as the quaternion formed by
     *  [0, 0, 1] and 270. This method favors the latter.
     */
    export function getAxisAngle(out_axis: Vec3, q: Quat) {
        let rad = Math.acos(q[3]) * 2.0;
        let s = Math.sin(rad / 2.0);
        if (s !== 0.0) {
            out_axis[0] = q[0] / s;
            out_axis[1] = q[1] / s;
            out_axis[2] = q[2] / s;
        } else {
            // If s is zero, return any axis (no rotation - axis does not matter)
            out_axis[0] = 1;
            out_axis[1] = 0;
            out_axis[2] = 0;
        }
        return rad;
    }

    export function multiply(out: Quat, a: Quat, b: Quat) {
        let ax = a[0], ay = a[1], az = a[2], aw = a[3];
        let bx = b[0], by = b[1], bz = b[2], bw = b[3];

        out[0] = ax * bw + aw * bx + ay * bz - az * by;
        out[1] = ay * bw + aw * by + az * bx - ax * bz;
        out[2] = az * bw + aw * bz + ax * by - ay * bx;
        out[3] = aw * bw - ax * bx - ay * by - az * bz;
        return out;
    }

    export function rotateX(out: Quat, a: Quat, rad: number) {
        rad *= 0.5;

        let ax = a[0], ay = a[1], az = a[2], aw = a[3];
        let bx = Math.sin(rad), bw = Math.cos(rad);

        out[0] = ax * bw + aw * bx;
        out[1] = ay * bw + az * bx;
        out[2] = az * bw - ay * bx;
        out[3] = aw * bw - ax * bx;
        return out;
    }

    export function rotateY(out: Quat, a: Quat, rad: number) {
        rad *= 0.5;

        let ax = a[0], ay = a[1], az = a[2], aw = a[3];
        let by = Math.sin(rad), bw = Math.cos(rad);

        out[0] = ax * bw - az * by;
        out[1] = ay * bw + aw * by;
        out[2] = az * bw + ax * by;
        out[3] = aw * bw - ay * by;
        return out;
    }

    export function rotateZ(out: Quat, a: Quat, rad: number) {
        rad *= 0.5;

        let ax = a[0], ay = a[1], az = a[2], aw = a[3];
        let bz = Math.sin(rad), bw = Math.cos(rad);

        out[0] = ax * bw + ay * bz;
        out[1] = ay * bw - ax * bz;
        out[2] = az * bw + aw * bz;
        out[3] = aw * bw - az * bz;
        return out;
    }

    /**
     * Calculates the W component of a quat from the X, Y, and Z components.
     * Assumes that quaternion is 1 unit in length.
     * Any existing W component will be ignored.
     */
    export function calculateW(out: Quat, a: Quat) {
        let x = a[0], y = a[1], z = a[2];

        out[0] = x;
        out[1] = y;
        out[2] = z;
        out[3] = Math.sqrt(Math.abs(1.0 - x * x - y * y - z * z));
        return out;
    }

    /**
     * Performs a spherical linear interpolation between two quat
     */
    export function slerp(out: Quat, a: Quat, b: Quat, t: number) {
        // benchmarks:
        //    http://jsperf.com/quaternion-slerp-implementations
        let ax = a[0], ay = a[1], az = a[2], aw = a[3];
        let bx = b[0], by = b[1], bz = b[2], bw = b[3];

        let omega, cosom, sinom, scale0, scale1;

        // calc cosine
        cosom = ax * bx + ay * by + az * bz + aw * bw;
        // adjust signs (if necessary)
        if ( cosom < 0.0 ) {
        cosom = -cosom;
        bx = - bx;
        by = - by;
        bz = - bz;
        bw = - bw;
        }
        // calculate coefficients
        if ( (1.0 - cosom) > 0.000001 ) {
            // standard case (slerp)
            omega  = Math.acos(cosom);
            sinom  = Math.sin(omega);
            scale0 = Math.sin((1.0 - t) * omega) / sinom;
            scale1 = Math.sin(t * omega) / sinom;
        } else {
            // "from" and "to" quaternions are very close
            //  ... so we can do a linear interpolation
            scale0 = 1.0 - t;
            scale1 = t;
        }
        // calculate final values
        out[0] = scale0 * ax + scale1 * bx;
        out[1] = scale0 * ay + scale1 * by;
        out[2] = scale0 * az + scale1 * bz;
        out[3] = scale0 * aw + scale1 * bw;

        return out;
    }

    export function invert(out: Quat, a: Quat) {
        let a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        let dot = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
        let invDot = dot ? 1.0/dot : 0;

        // TODO: Would be faster to return [0,0,0,0] immediately if dot == 0

        out[0] = -a0 * invDot;
        out[1] = -a1 * invDot;
        out[2] = -a2 * invDot;
        out[3] = a3 * invDot;
        return out;
    }

    /**
     * Calculates the conjugate of a quat
     * If the quaternion is normalized, this function is faster than quat.inverse and produces the same result.
     */
    export function conjugate(out: Quat, a: Quat) {
        out[0] = -a[0];
        out[1] = -a[1];
        out[2] = -a[2];
        out[3] = a[3];
        return out;
    }

    /**
     * Creates a quaternion from the given 3x3 rotation matrix.
     *
     * NOTE: The resultant quaternion is not normalized, so you should be sure
     * to renormalize the quaternion yourself where necessary.
     */
    export function fromMat3(out: Quat, m: Mat3) {
        // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
        // article "Quaternion Calculus and Fast Animation".
        const fTrace = m[0] + m[4] + m[8];
        let fRoot;

        if ( fTrace > 0.0 ) {
            // |w| > 1/2, may as well choose w > 1/2
            fRoot = Math.sqrt(fTrace + 1.0);  // 2w
            out[3] = 0.5 * fRoot;
            fRoot = 0.5/fRoot;  // 1/(4w)
            out[0] = (m[5]-m[7])*fRoot;
            out[1] = (m[6]-m[2])*fRoot;
            out[2] = (m[1]-m[3])*fRoot;
            } else {
            // |w| <= 1/2
            let i = 0;
            if ( m[4] > m[0] ) i = 1;
            if ( m[8] > m[i*3+i] ) i = 2;
            let j = (i+1)%3;
            let k = (i+2)%3;

            fRoot = Math.sqrt(m[i*3+i]-m[j*3+j]-m[k*3+k] + 1.0);
            out[i] = 0.5 * fRoot;
            fRoot = 0.5 / fRoot;
            out[3] = (m[j*3+k] - m[k*3+j]) * fRoot;
            out[j] = (m[j*3+i] + m[i*3+j]) * fRoot;
            out[k] = (m[k*3+i] + m[i*3+k]) * fRoot;
        }

        return out;
    }

    export function clone(a: Quat) {
        const out = zero();
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        out[3] = a[3];
        return out;
    }

    export function copy(out: Quat, a: Quat) {
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        out[3] = a[3];
        return out;
    }

    export function set(out: Quat, x: number, y: number, z: number, w: number) {
        out[0] = x;
        out[1] = y;
        out[2] = z;
        out[3] = w;
        return out;
    }

    export function add(out: Quat, a: Quat, b: Quat) {
        out[0] = a[0] + b[0];
        out[1] = a[1] + b[1];
        out[2] = a[2] + b[2];
        out[3] = a[3] + b[3];
        return out;
    }

    export function normalize(out: Quat, a: Quat) {
        let x = a[0];
        let y = a[1];
        let z = a[2];
        let w = a[3];
        let len = x*x + y*y + z*z + w*w;
        if (len > 0) {
            len = 1 / Math.sqrt(len);
            out[0] = x * len;
            out[1] = y * len;
            out[2] = z * len;
            out[3] = w * len;
        }
        return out;
    }

    /**
     * Sets a quaternion to represent the shortest rotation from one
     * vector to another.
     *
     * Both vectors are assumed to be unit length.
     */
    const rotTmpVec3 = Vec3.zero();
    const rotTmpVec3UnitX = Vec3.create(1, 0, 0);
    const rotTmpVec3UnitY = Vec3.create(0, 1, 0);
    export function rotationTo(out: Quat, a: Vec3, b: Vec3) {
        let dot = Vec3.dot(a, b);
        if (dot < -0.999999) {
            Vec3.cross(rotTmpVec3, rotTmpVec3UnitX, a);
            if (Vec3.magnitude(rotTmpVec3) < 0.000001)
            Vec3.cross(rotTmpVec3, rotTmpVec3UnitY, a);
            Vec3.normalize(rotTmpVec3, rotTmpVec3);
            setAxisAngle(out, rotTmpVec3, Math.PI);
            return out;
        } else if (dot > 0.999999) {
            out[0] = 0;
            out[1] = 0;
            out[2] = 0;
            out[3] = 1;
            return out;
        } else {
            Vec3.cross(rotTmpVec3, a, b);
            out[0] = rotTmpVec3[0];
            out[1] = rotTmpVec3[1];
            out[2] = rotTmpVec3[2];
            out[3] = 1 + dot;
            return normalize(out, out);
        }
    }

    /**
     * Performs a spherical linear interpolation with two control points
     */
    let sqlerpTemp1 = Quat.zero();
    let sqlerpTemp2 = Quat.zero();
    export function sqlerp(out: Quat, a: Quat, b: Quat, c: Quat, d: Quat, t: number) {
        slerp(sqlerpTemp1, a, d, t);
        slerp(sqlerpTemp2, b, c, t);
        slerp(out, sqlerpTemp1, sqlerpTemp2, 2 * t * (1 - t));
        return out;
    }

    /**
     * Sets the specified quaternion with values corresponding to the given
     * axes. Each axis is a vec3 and is expected to be unit length and
     * perpendicular to all other specified axes.
     */
    const axesTmpMat = Mat3.zero();
    export function setAxes(out: Quat, view: Vec3, right: Vec3, up: Vec3) {
        axesTmpMat[0] = right[0];
        axesTmpMat[3] = right[1];
        axesTmpMat[6] = right[2];

        axesTmpMat[1] = up[0];
        axesTmpMat[4] = up[1];
        axesTmpMat[7] = up[2];

        axesTmpMat[2] = -view[0];
        axesTmpMat[5] = -view[1];
        axesTmpMat[8] = -view[2];

        return normalize(out, Quat.fromMat3(out, axesTmpMat));
    }
}