/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
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

export type Mat4 = number[]
export type Vec3 = number[]
export type Vec4 = number[]

const enum EPSILON { Value = 0.000001 }

export function Mat4() {
    return Mat4.zero();
}

/**
 * Stores a 4x4 matrix in a column major (j * 4 + i indexing) format.
 */
export namespace Mat4 {
    export function zero(): number[] {
        // force double backing array by 0.1.
        const ret = [0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        ret[0] = 0.0;
        return ret;
    }

    export function identity(): number[] {
        let out = zero();
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

    export function fromIdentity(mat: number[]): number[] {
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

    export function ofRows(rows: number[][]): number[] {
        let out = zero(), i: number, j: number, r: number[];
        for (i = 0; i < 4; i++) {
            r = rows[i];
            for (j = 0; j < 4; j++) {
                out[4 * j + i] = r[j];
            }
        }
        return out;
    }

    export function areEqual(a: number[], b: number[], eps: number) {
        for (let i = 0; i < 16; i++) {
            if (Math.abs(a[i] - b[i]) > eps) {
                return false;
            }
        }
        return true;
    }

    export function setValue(a: number[], i: number, j: number, value: number) {
        a[4 * j + i] = value;
    }

    export function copy(out: number[], a: number[]) {
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

    export function clone(a: number[]) {
        return Mat4.copy(Mat4.zero(), a);
    }

    export function invert(out: number[], a: number[]) {
        let a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
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
            b11 = a22 * a33 - a23 * a32,

            // Calculate the determinant
            det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

        if (!det) {
            return null;
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

    export function mul(out: number[], a: number[], b: number[]) {
        let a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
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

    export function mul3(out: number[], a: number[], b: number[], c: number[]) {
        return mul(out, mul(out, a, b), c);
    }

    export function translate(out: number[], a: number[], v: number[]) {
        let x = v[0], y = v[1], z = v[2],
            a00: number, a01: number, a02: number, a03: number,
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

    export function fromTranslation(out: number[], v: number[]) {
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

    export function rotate(out: number[], a: number[], rad: number, axis: number[]) {
        let x = axis[0], y = axis[1], z = axis[2],
            len = Math.sqrt(x * x + y * y + z * z),
            s, c, t,
            a00, a01, a02, a03,
            a10, a11, a12, a13,
            a20, a21, a22, a23,
            b00, b01, b02,
            b10, b11, b12,
            b20, b21, b22;

        if (Math.abs(len) < EPSILON.Value) { return null; }

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

    export function fromRotation(out: number[], rad: number, axis: number[]) {
        let x = axis[0], y = axis[1], z = axis[2],
            len = Math.sqrt(x * x + y * y + z * z),
            s, c, t;

        if (Math.abs(len) < EPSILON.Value) { return fromIdentity(out); }

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

    export function scale(out: number[], a: number[], v: number[]) {
        let x = v[0], y = v[1], z = v[2];

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

    export function fromScaling(out: number[], v: number[]) {
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

    export function makeTable(m: number[]) {
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

    export function determinant(a: number[]) {
        let a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
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
}

export function Vec3(x?: number, y?: number, z?: number) {
    return Vec3.fromValues(x || 0, y || 0, z || 0);
}

export namespace Vec3 {
    export function zero() {
        let out = [0.1, 0.0, 0.0];
        out[0] = 0;
        return out;
    }

    export function clone(a: number[]) {
        let out = zero();
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        return out;
    }

    export function fromObj(v: { x: number, y: number, z: number }) {
        return fromValues(v.x, v.y, v.z);
    }

    export function toObj(v: number[]) {
        return { x: v[0], y: v[1], z: v[2] };
    }

    export function fromValues(x: number, y: number, z: number) {
        let out = zero();
        out[0] = x;
        out[1] = y;
        out[2] = z;
        return out;
    }

    export function set(out: number[], x: number, y: number, z: number) {
        out[0] = x;
        out[1] = y;
        out[2] = z;
        return out;
    }

    export function copy(out: number[], a: number[]) {
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        return out;
    }

    export function add(out: number[], a: number[], b: number[]) {
        out[0] = a[0] + b[0];
        out[1] = a[1] + b[1];
        out[2] = a[2] + b[2];
        return out;
    }

    export function sub(out: number[], a: number[], b: number[]) {
        out[0] = a[0] - b[0];
        out[1] = a[1] - b[1];
        out[2] = a[2] - b[2];
        return out;
    }

    export function scale(out: number[], a: number[], b: number) {
        out[0] = a[0] * b;
        out[1] = a[1] * b;
        out[2] = a[2] * b;
        return out;
    }

    export function scaleAndAdd(out: number[], a: number[], b: number[], scale: number) {
        out[0] = a[0] + (b[0] * scale);
        out[1] = a[1] + (b[1] * scale);
        out[2] = a[2] + (b[2] * scale);
        return out;
    }

    export function distance(a: number[], b: number[]) {
        let x = b[0] - a[0],
            y = b[1] - a[1],
            z = b[2] - a[2];
        return Math.sqrt(x * x + y * y + z * z);
    }

    export function squaredDistance(a: number[], b: number[]) {
        let x = b[0] - a[0],
            y = b[1] - a[1],
            z = b[2] - a[2];
        return x * x + y * y + z * z;
    }

    export function magnitude(a: number[]) {
        let x = a[0],
            y = a[1],
            z = a[2];
        return Math.sqrt(x * x + y * y + z * z);
    }

    export function squaredMagnitude(a: number[]) {
        let x = a[0],
            y = a[1],
            z = a[2];
        return x * x + y * y + z * z;
    }

    export function normalize(out: number[], a: number[]) {
        let x = a[0],
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

    export function dot(a: number[], b: number[]) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    export function cross(out: number[], a: number[], b: number[]) {
        let ax = a[0], ay = a[1], az = a[2],
            bx = b[0], by = b[1], bz = b[2];

        out[0] = ay * bz - az * by;
        out[1] = az * bx - ax * bz;
        out[2] = ax * by - ay * bx;
        return out;
    }

    export function lerp(out: number[], a: number[], b: number[], t: number) {
        let ax = a[0],
            ay = a[1],
            az = a[2];
        out[0] = ax + t * (b[0] - ax);
        out[1] = ay + t * (b[1] - ay);
        out[2] = az + t * (b[2] - az);
        return out;
    }

    export function transformMat4(out: number[], a: number[], m: number[]) {
        let x = a[0], y = a[1], z = a[2],
            w = (m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0;
        out[0] = (m[0] * x + m[4] * y + m[8] * z + m[12]) / w;
        out[1] = (m[1] * x + m[5] * y + m[9] * z + m[13]) / w;
        out[2] = (m[2] * x + m[6] * y + m[10] * z + m[14]) / w;
        return out;
    }

    const angleTempA = zero(), angleTempB = zero();
    export function angle(a: number[], b: number[]) {
        copy(angleTempA, a);
        copy(angleTempB, b);

        normalize(angleTempA, angleTempA);
        normalize(angleTempB, angleTempB);

        let cosine = dot(angleTempA, angleTempB);

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
        if (Math.abs(by) < 0.0001) return Mat4.fromIdentity(mat);
        const axis = cross(rotTemp, a, b);
        return Mat4.fromRotation(mat, by, axis);
    }
}

export function Vec4(x?: number, y?: number, z?: number, w?: number) {
    return Vec4.fromValues(x || 0, y || 0, z || 0, w || 0);
}

export namespace Vec4 {
    export function zero(): number[] {
        // force double backing array by 0.1.
        const ret = [0.1, 0, 0, 0];
        ret[0] = 0.0;
        return ret;
    }

    export function clone(a: number[]) {
        let out = zero();
        out[0] = a[0];
        out[1] = a[1];
        out[2] = a[2];
        out[3] = a[3];
        return out;
    }

    export function fromValues(x: number, y: number, z: number, w: number) {
        let out = zero();
        out[0] = x;
        out[1] = y;
        out[2] = z;
        out[3] = w;
        return out;
    }

    export function set(out: number[], x: number, y: number, z: number, w: number) {
        out[0] = x;
        out[1] = y;
        out[2] = z;
        out[3] = w;
        return out;
    }

    export function distance(a: number[], b: number[]) {
        let x = b[0] - a[0],
            y = b[1] - a[1],
            z = b[2] - a[2],
            w = b[3] - a[3];
        return Math.sqrt(x * x + y * y + z * z + w * w);
    }

    export function squaredDistance(a: number[], b: number[]) {
        let x = b[0] - a[0],
            y = b[1] - a[1],
            z = b[2] - a[2],
            w = b[3] - a[3];
        return x * x + y * y + z * z + w * w;
    }

    export function norm(a: number[]) {
        let x = a[0],
            y = a[1],
            z = a[2],
            w = a[3];
        return Math.sqrt(x * x + y * y + z * z + w * w);
    }

    export function squaredNorm(a: number[]) {
        let x = a[0],
            y = a[1],
            z = a[2],
            w = a[3];
        return x * x + y * y + z * z + w * w;
    }

    export function transform(out: number[], a: number[], m: number[]) {
        let x = a[0], y = a[1], z = a[2], w = a[3];
        out[0] = m[0] * x + m[4] * y + m[8] * z + m[12] * w;
        out[1] = m[1] * x + m[5] * y + m[9] * z + m[13] * w;
        out[2] = m[2] * x + m[6] * y + m[10] * z + m[14] * w;
        out[3] = m[3] * x + m[7] * y + m[11] * z + m[15] * w;
        return out;
    }
}