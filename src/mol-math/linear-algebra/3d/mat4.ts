/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

import { EPSILON, equalEps } from './common';
import Vec3 from './vec3';
import Quat from './quat';
import { degToRad } from '../../misc';
import { NumberArray } from '../../../mol-util/type-helpers';
import Mat3 from './mat3';

interface Mat4 extends Array<number> { [d: number]: number, '@type': 'mat4', length: 16 }
interface ReadonlyMat4 extends Array<number> { readonly [d: number]: number, '@type': 'mat4', length: 16 }

function Mat4() {
    return Mat4.zero();
}

/**
 * Stores a 4x4 matrix in a column major (j * 4 + i indexing) format.
 */
namespace Mat4 {
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

    export function setZero(mat: Mat4): Mat4 {
        for (let i = 0; i < 16; i++) mat[i] = 0;
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
        return areEqual(m, _id, typeof eps === 'undefined' ? EPSILON : eps);
    }

    export function hasNaN(m: Mat4) {
        for (let i = 0; i < 16; i++) if (isNaN(m[i])) return true;
        return false;
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

    export function getValue(a: Mat4, i: number, j: number) {
        return a[4 * j + i];
    }

    export function toArray(a: Mat4, out: NumberArray, offset: number) {
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
        return out;
    }

    export function fromArray(a: Mat4, array: NumberArray, offset: number) {
        a[0] = array[offset + 0];
        a[1] = array[offset + 1];
        a[2] = array[offset + 2];
        a[3] = array[offset + 3];
        a[4] = array[offset + 4];
        a[5] = array[offset + 5];
        a[6] = array[offset + 6];
        a[7] = array[offset + 7];
        a[8] = array[offset + 8];
        a[9] = array[offset + 9];
        a[10] = array[offset + 10];
        a[11] = array[offset + 11];
        a[12] = array[offset + 12];
        a[13] = array[offset + 13];
        a[14] = array[offset + 14];
        a[15] = array[offset + 15];
        return a;
    }

    export function fromBasis(a: Mat4, x: Vec3, y: Vec3, z: Vec3) {
        Mat4.setZero(a);
        Mat4.setValue(a, 0, 0, x[0]);
        Mat4.setValue(a, 1, 0, x[1]);
        Mat4.setValue(a, 2, 0, x[2]);
        Mat4.setValue(a, 0, 1, y[0]);
        Mat4.setValue(a, 1, 1, y[1]);
        Mat4.setValue(a, 2, 1, y[2]);
        Mat4.setValue(a, 0, 2, z[0]);
        Mat4.setValue(a, 1, 2, z[1]);
        Mat4.setValue(a, 2, 2, z[2]);
        Mat4.setValue(a, 3, 3, 1);
        return a;
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

    /**
     * Returns the translation vector component of a transformation matrix.
     */
    export function getTranslation(out: Vec3, mat: Mat4) {
        out[0] = mat[12];
        out[1] = mat[13];
        out[2] = mat[14];
        return out;
    }

    /**
     * Returns the scaling factor component of a transformation matrix.
     */
    export function getScaling(out: Vec3, mat: Mat4) {
        let m11 = mat[0];
        let m12 = mat[1];
        let m13 = mat[2];
        let m21 = mat[4];
        let m22 = mat[5];
        let m23 = mat[6];
        let m31 = mat[8];
        let m32 = mat[9];
        let m33 = mat[10];
        out[0] = Math.sqrt(m11 * m11 + m12 * m12 + m13 * m13);
        out[1] = Math.sqrt(m21 * m21 + m22 * m22 + m23 * m23);
        out[2] = Math.sqrt(m31 * m31 + m32 * m32 + m33 * m33);
        return out;
    }

    /**
     * Returns a quaternion representing the rotational component of a transformation matrix.
     */
    export function getRotation(out: Quat, mat: Mat4) {
        // Algorithm taken from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
        let trace = mat[0] + mat[5] + mat[10];
        let S = 0;

        if (trace > 0) {
            S = Math.sqrt(trace + 1.0) * 2;
            out[3] = 0.25 * S;
            out[0] = (mat[6] - mat[9]) / S;
            out[1] = (mat[8] - mat[2]) / S;
            out[2] = (mat[1] - mat[4]) / S;
        } else if ((mat[0] > mat[5]) && (mat[0] > mat[10])) {
            S = Math.sqrt(1.0 + mat[0] - mat[5] - mat[10]) * 2;
            out[3] = (mat[6] - mat[9]) / S;
            out[0] = 0.25 * S;
            out[1] = (mat[1] + mat[4]) / S;
            out[2] = (mat[8] + mat[2]) / S;
        } else if (mat[5] > mat[10]) {
            S = Math.sqrt(1.0 + mat[5] - mat[0] - mat[10]) * 2;
            out[3] = (mat[8] - mat[2]) / S;
            out[0] = (mat[1] + mat[4]) / S;
            out[1] = 0.25 * S;
            out[2] = (mat[6] + mat[9]) / S;
        } else {
            S = Math.sqrt(1.0 + mat[10] - mat[0] - mat[5]) * 2;
            out[3] = (mat[1] - mat[4]) / S;
            out[0] = (mat[8] + mat[2]) / S;
            out[1] = (mat[6] + mat[9]) / S;
            out[2] = 0.25 * S;
        }

        return out;
    }

    export function extractRotation(out: Mat4, mat: Mat4) {
        const scaleX = 1 / Math.sqrt(mat[0] * mat[0] + mat[1] * mat[1] + mat[2] * mat[2]);
        const scaleY = 1 / Math.sqrt(mat[4] * mat[4] + mat[5] * mat[5] + mat[6] * mat[6]);
        const scaleZ = 1 / Math.sqrt(mat[8] * mat[8] + mat[9] * mat[9] + mat[10] * mat[10]);

        out[0] = mat[0] * scaleX;
        out[1] = mat[1] * scaleX;
        out[2] = mat[2] * scaleX;
        out[3] = 0;
        out[4] = mat[4] * scaleY;
        out[5] = mat[5] * scaleY;
        out[6] = mat[6] * scaleY;
        out[7] = 0;
        out[8] = mat[8] * scaleZ;
        out[9] = mat[9] * scaleZ;
        out[10] = mat[10] * scaleZ;
        out[11] = 0;
        out[12] = 0;
        out[13] = 0;
        out[14] = 0;
        out[15] = 1;

        return out;
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

    /**
     * Like `mul` but with offsets into arrays
     */
    export function mulOffset(out: NumberArray, a: NumberArray, b: NumberArray, oOut: number, oA: number, oB: number) {
        const a00 = a[0 + oA], a01 = a[1 + oA], a02 = a[2 + oA], a03 = a[3 + oA],
            a10 = a[4 + oA], a11 = a[5 + oA], a12 = a[6 + oA], a13 = a[7 + oA],
            a20 = a[8 + oA], a21 = a[9 + oA], a22 = a[10 + oA], a23 = a[11 + oA],
            a30 = a[12 + oA], a31 = a[13 + oA], a32 = a[14 + oA], a33 = a[15 + oA];

        // Cache only the current line of the second matrix
        let b0 = b[0 + oB], b1 = b[1 + oB], b2 = b[2 + oB], b3 = b[3 + oB];
        out[0 + oOut] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
        out[1 + oOut] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
        out[2 + oOut] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
        out[3 + oOut] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

        b0 = b[4 + oB]; b1 = b[5 + oB]; b2 = b[6 + oB]; b3 = b[7 + oB];
        out[4 + oOut] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
        out[5 + oOut] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
        out[6 + oOut] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
        out[7 + oOut] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

        b0 = b[8 + oB]; b1 = b[9 + oB]; b2 = b[10 + oB]; b3 = b[11 + oB];
        out[8 + oOut] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
        out[9 + oOut] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
        out[10 + oOut] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
        out[11 + oOut] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;

        b0 = b[12 + oB]; b1 = b[13 + oB]; b2 = b[14 + oB]; b3 = b[15 + oB];
        out[12 + oOut] = b0 * a00 + b1 * a10 + b2 * a20 + b3 * a30;
        out[13 + oOut] = b0 * a01 + b1 * a11 + b2 * a21 + b3 * a31;
        out[14 + oOut] = b0 * a02 + b1 * a12 + b2 * a22 + b3 * a32;
        out[15 + oOut] = b0 * a03 + b1 * a13 + b2 * a23 + b3 * a33;
        return out;
    }

    export function mul3(out: Mat4, a: Mat4, b: Mat4, c: Mat4) {
        return mul(out, mul(out, a, b), c);
    }

    /** Translate a Mat4 by the given Vec3 */
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

    /**
     * Sets the specified quaternion with values corresponding to the given
     * axes. Each axis is a vec3 and is expected to be unit length and
     * perpendicular to all other specified axes.
     */
    export function setAxes(out: Mat4, view: Vec3, right: Vec3, up: Vec3) {
        out[0] = right[0];
        out[4] = right[1];
        out[8] = right[2];
        out[1] = up[0];
        out[5] = up[1];
        out[9] = up[2];
        out[2] = view[0];
        out[6] = view[1];
        out[10] = view[2];
        return out;
    }

    export function rotate(out: Mat4, a: Mat4, rad: number, axis: Vec3) {
        let x = axis[0], y = axis[1], z = axis[2],
            len = Math.sqrt(x * x + y * y + z * z),
            s, c, t,
            a00, a01, a02, a03,
            a10, a11, a12, a13,
            a20, a21, a22, a23,
            b00, b01, b02,
            b10, b11, b12,
            b20, b21, b22;

        if (Math.abs(len) < EPSILON) {
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

        if (Math.abs(len) < EPSILON) { return setIdentity(out); }

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

    export function scaleUniformly(out: Mat4, a: Mat4, scale: number) {
        out[0] = a[0] * scale;
        out[1] = a[1] * scale;
        out[2] = a[2] * scale;
        out[3] = a[3] * scale;
        out[4] = a[4] * scale;
        out[5] = a[5] * scale;
        out[6] = a[6] * scale;
        out[7] = a[7] * scale;
        out[8] = a[8] * scale;
        out[9] = a[9] * scale;
        out[10] = a[10] * scale;
        out[11] = a[11] * scale;
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

    export function fromUniformScaling(out: Mat4, scale: number) {
        out[0] = scale;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = scale;
        out[6] = 0;
        out[7] = 0;
        out[8] = 0;
        out[9] = 0;
        out[10] = scale;
        out[11] = 0;
        out[12] = 0;
        out[13] = 0;
        out[14] = 0;
        out[15] = 1;
        return out;
    }

    /**
     * Copies the mat3 into upper-left 3x3 values.
     */
    export function fromMat3(out: Mat4, a: Mat3) {
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        out[4] = a[3];
        out[5] = a[4];
        out[6] = a[5];
        out[8] = a[6];
        out[9] = a[7];
        out[10] = a[8];
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
     *
     * Allows for improper rotations
     */
    export function isRotationAndTranslation(a: Mat4, eps?: number) {
        return _isRotationAndTranslation(a, typeof eps !== 'undefined' ? eps : EPSILON);
    }

    function _isRotationAndTranslation(a: Mat4, eps: number) {
        const a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
            a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7],
            a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11],
            a33 = a[15];

        if (!equalEps(a33, 1, eps) || !equalEps(a03, 0, eps) || !equalEps(a13, 0, eps) || !equalEps(a23, 0, eps)) {
            return false;
        }

        // use `abs` to allow for improper rotations
        const det3x3 = Math.abs(a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20) + a02 * (a10 * a21 - a11 * a20));
        if (!equalEps(det3x3, 1, eps)) {
            return false;
        }
        return true;
    }

    /**
     * Check if the matrix has only translation and uniform scaling
     * [ S  0  0  X ]
     * [ 0  S  0  Y ]
     * [ 0  0  S  Z ]
     * [ 0  0  0  1 ]
     */
    export function isTranslationAndUniformScaling(a: Mat4, eps?: number) {
        return _isTranslationAndUniformScaling(a, typeof eps !== 'undefined' ? eps : EPSILON);
    }

    function _isTranslationAndUniformScaling(a: Mat4, eps: number) {
        const a00 = a[0];
        return (
            // 0 base scaling
            equalEps(a[1], 0, eps) &&
            equalEps(a[2], 0, eps) &&
            equalEps(a[3], 0, eps) &&
            equalEps(a[4], 0, eps) &&
            equalEps(a[5], a00, eps) &&
            equalEps(a[6], 0, eps) &&
            equalEps(a[7], 0, eps) &&
            equalEps(a[8], 0, eps) &&
            equalEps(a[9], 0, eps) &&
            equalEps(a[10], a00, eps) &&
            equalEps(a[11], 0, eps) &&
            // 12, 13, 14 translation can be anything
            equalEps(a[15], 1, eps)
        );
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
     * Generates a perspective projection (frustum) matrix with the given bounds
     */
    export function perspective(out: Mat4, left: number, right: number, top: number, bottom: number, near: number, far: number) {
        const x = 2 * near / (right - left);
        const y = 2 * near / (top - bottom);

        const a = (right + left) / (right - left);
        const b = (top + bottom) / (top - bottom);
        const c = -(far + near) / (far - near);
        const d = -2 * far * near / (far - near);

        out[0] = x;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = y;
        out[6] = 0;
        out[7] = 0;
        out[8] = a;
        out[9] = b;
        out[10] = c;
        out[11] = -1;
        out[12] = 0;
        out[13] = 0;
        out[14] = d;
        out[15] = 0;
        return out;
    }

    /**
     * Generates a orthogonal projection matrix with the given bounds
     */
    export function ortho(out: Mat4, left: number, right: number, top: number, bottom: number, near: number, far: number) {
        const w = 1.0 / (right - left);
        const h = 1.0 / (top - bottom);
        const p = 1.0 / (far - near);

        const x = (right + left) * w;
        const y = (top + bottom) * h;
        const z = (far + near) * p;

        out[0] = 2 * w;
        out[1] = 0;
        out[2] = 0;
        out[3] = 0;
        out[4] = 0;
        out[5] = 2 * h;
        out[6] = 0;
        out[7] = 0;
        out[8] = 0;
        out[9] = 0;
        out[10] = -2 * p;
        out[11] = 0;
        out[12] = -x;
        out[13] = -y;
        out[14] = -z;
        out[15] = 1;
        return out;
    }

    /**
     * Generates a look-at matrix with the given eye position, focal point, and up axis
     */
    export function lookAt(out: Mat4, eye: Vec3, center: Vec3, up: Vec3) {
        let x0, x1, x2, y0, y1, y2, z0, z1, z2, len;
        const eyex = eye[0];
        const eyey = eye[1];
        const eyez = eye[2];
        const upx = up[0];
        const upy = up[1];
        const upz = up[2];
        const centerx = center[0];
        const centery = center[1];
        const centerz = center[2];

        if (Math.abs(eyex - centerx) < EPSILON &&
            Math.abs(eyey - centery) < EPSILON &&
            Math.abs(eyez - centerz) < EPSILON
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

    /**
     * Generates a matrix that makes something look at something else.
     */
    export function targetTo(out: Mat4, eye: Vec3, target: Vec3, up: Vec3) {
        const eyex = eye[0],
            eyey = eye[1],
            eyez = eye[2],
            upx = up[0],
            upy = up[1],
            upz = up[2];

        let z0 = eyex - target[0],
            z1 = eyey - target[1],
            z2 = eyez - target[2];

        let len = z0 * z0 + z1 * z1 + z2 * z2;
        if (len > 0) {
            len = 1 / Math.sqrt(len);
            z0 *= len;
            z1 *= len;
            z2 *= len;
        }

        let x0 = upy * z2 - upz * z1,
            x1 = upz * z0 - upx * z2,
            x2 = upx * z1 - upy * z0;

        len = x0 * x0 + x1 * x1 + x2 * x2;
        if (len > 0) {
            len = 1 / Math.sqrt(len);
            x0 *= len;
            x1 *= len;
            x2 *= len;
        }

        out[0] = x0;
        out[1] = x1;
        out[2] = x2;
        out[3] = 0;
        out[4] = z1 * x2 - z2 * x1;
        out[5] = z2 * x0 - z0 * x2;
        out[6] = z0 * x1 - z1 * x0;
        out[7] = 0;
        out[8] = z0;
        out[9] = z1;
        out[10] = z2;
        out[11] = 0;
        out[12] = eyex;
        out[13] = eyey;
        out[14] = eyez;
        out[15] = 1;
        return out;
    }

    /**
     * Perm is 0-indexed permutation
     */
    export function fromPermutation(out: Mat4, perm: number[]) {
        setZero(out);
        for (let i = 0; i < 4; i++) {
            const p = perm[i];
            setValue(out, i, p, 1);
        }
        return out;
    }

    export function getMaxScaleOnAxis(m: Mat4) {
        const scaleXSq = m[0] * m[0] + m[1] * m[1] + m[2] * m[2];
        const scaleYSq = m[4] * m[4] + m[5] * m[5] + m[6] * m[6];
        const scaleZSq = m[8] * m[8] + m[9] * m[9] + m[10] * m[10];
        return Math.sqrt(Math.max(scaleXSq, scaleYSq, scaleZSq));
    }

    const xAxis = Vec3.create(1, 0, 0);
    const yAxis = Vec3.create(0, 1, 0);
    const zAxis = Vec3.create(0, 0, 1);

    /** Rotation matrix for 90deg around x-axis */
    export const rotX90: ReadonlyMat4 = Mat4.fromRotation(Mat4(), degToRad(90), xAxis);
    /** Rotation matrix for 180deg around x-axis */
    export const rotX180: ReadonlyMat4 = Mat4.fromRotation(Mat4(), degToRad(180), xAxis);
    /** Rotation matrix for 90deg around y-axis */
    export const rotY90: ReadonlyMat4 = Mat4.fromRotation(Mat4(), degToRad(90), yAxis);
    /** Rotation matrix for 180deg around y-axis */
    export const rotY180: ReadonlyMat4 = Mat4.fromRotation(Mat4(), degToRad(180), yAxis);
    /** Rotation matrix for 270deg around y-axis */
    export const rotY270: ReadonlyMat4 = Mat4.fromRotation(Mat4(), degToRad(270), yAxis);
    /** Rotation matrix for 90deg around z-axis */
    export const rotZ90: ReadonlyMat4 = Mat4.fromRotation(Mat4(), degToRad(90), zAxis);
    /** Rotation matrix for 180deg around z-axis */
    export const rotZ180: ReadonlyMat4 = Mat4.fromRotation(Mat4(), degToRad(180), zAxis);
    /** Rotation matrix for 90deg around first x-axis and then y-axis */
    export const rotXY90: ReadonlyMat4 = Mat4.mul(Mat4(), rotX90, rotY90);
    /** Rotation matrix for 90deg around first z-axis and then y-axis */
    export const rotZY90: ReadonlyMat4 = Mat4.mul(Mat4(), rotZ90, rotY90);
    /** Rotation matrix for 90deg around first z-axis and then y-axis and then z-axis */
    export const rotZYZ90: ReadonlyMat4 = Mat4.mul(Mat4(), rotZY90, rotZ90);
    /** Rotation matrix for 90deg around first z-axis and then 180deg around x-axis */
    export const rotZ90X180: ReadonlyMat4 = Mat4.mul(Mat4(), rotZ90, rotX180);
    /** Rotation matrix for 90deg around first y-axis and then 180deg around z-axis */
    export const rotY90Z180: ReadonlyMat4 = Mat4.mul(Mat4(), rotY90, rotZ180);

    /** Identity matrix */
    export const id: ReadonlyMat4 = Mat4.identity();
}

export default Mat4;