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

import Mat4 from './mat4';
import { Quat, Mat3, EPSILON } from '../3d';
import { spline as _spline, quadraticBezier as _quadraticBezier, clamp } from '../../interpolate';
import { NumberArray } from '../../../mol-util/type-helpers';

export { ReadonlyVec3 };

interface Vec3 extends Array<number> { [d: number]: number, '@type': 'vec3', length: 3 }
interface ReadonlyVec3 extends Array<number> { readonly [d: number]: number, '@type': 'vec3', length: 3 }

function Vec3() {
    return Vec3.zero();
}

namespace Vec3 {
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

    export function hasNaN(a: Vec3) {
        return isNaN(a[0]) || isNaN(a[1]) || isNaN(a[2]);
    }

    export function setNaN(out: Vec3) {
        out[0] = NaN;
        out[1] = NaN;
        out[2] = NaN;
        return out;
    }

    export function fromObj(v: { x: number, y: number, z: number }): Vec3 {
        return create(v.x, v.y, v.z);
    }

    export function toObj(v: Vec3) {
        return { x: v[0], y: v[1], z: v[2] };
    }

    export function fromArray(v: Vec3, array: ArrayLike<number>, offset: number) {
        v[0] = array[offset + 0];
        v[1] = array[offset + 1];
        v[2] = array[offset + 2];
        return v;
    }

    export function toArray(v: Vec3, out: NumberArray, offset: number) {
        out[offset + 0] = v[0];
        out[offset + 1] = v[1];
        out[offset + 2] = v[2];
        return out;
    }

    export function create(x: number, y: number, z: number): Vec3 {
        const out = zero();
        out[0] = x;
        out[1] = y;
        out[2] = z;
        return out;
    }

    export function ofArray(array: ArrayLike<number>) {
        const out = zero();
        out[0] = array[0];
        out[1] = array[1];
        out[2] = array[2];
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

    export function mul(out: Vec3, a: Vec3, b: Vec3) {
        out[0] = a[0] * b[0];
        out[1] = a[1] * b[1];
        out[2] = a[2] * b[2];
        return out;
    }

    export function div(out: Vec3, a: Vec3, b: Vec3) {
        out[0] = a[0] / b[0];
        out[1] = a[1] / b[1];
        out[2] = a[2] / b[2];
        return out;
    }

    export function scale(out: Vec3, a: Vec3, b: number) {
        out[0] = a[0] * b;
        out[1] = a[1] * b;
        out[2] = a[2] * b;
        return out;
    }

    /** Scales b, then adds a and b together */
    export function scaleAndAdd(out: Vec3, a: Vec3, b: Vec3, scale: number) {
        out[0] = a[0] + (b[0] * scale);
        out[1] = a[1] + (b[1] * scale);
        out[2] = a[2] + (b[2] * scale);
        return out;
    }

    /** Scales b, then subtracts b from a */
    export function scaleAndSub(out: Vec3, a: Vec3, b: Vec3, scale: number) {
        out[0] = a[0] - (b[0] * scale);
        out[1] = a[1] - (b[1] * scale);
        out[2] = a[2] - (b[2] * scale);
        return out;
    }

    /**
     * Math.round the components of a Vec3
     */
    export function round(out: Vec3, a: Vec3) {
        out[0] = Math.round(a[0]);
        out[1] = Math.round(a[1]);
        out[2] = Math.round(a[2]);
        return out;
    }

    /**
     * Math.ceil the components of a Vec3
     */
    export function ceil(out: Vec3, a: Vec3) {
        out[0] = Math.ceil(a[0]);
        out[1] = Math.ceil(a[1]);
        out[2] = Math.ceil(a[2]);
        return out;
    }

    /**
     * Math.floor the components of a Vec3
     */
    export function floor(out: Vec3, a: Vec3) {
        out[0] = Math.floor(a[0]);
        out[1] = Math.floor(a[1]);
        out[2] = Math.floor(a[2]);
        return out;
    }

    /**
     * Returns the minimum of two Vec3's
     */
    export function min(out: Vec3, a: Vec3, b: Vec3) {
        out[0] = Math.min(a[0], b[0]);
        out[1] = Math.min(a[1], b[1]);
        out[2] = Math.min(a[2], b[2]);
        return out;
    }

    /**
     * Returns the maximum of two Vec3's
     */
    export function max(out: Vec3, a: Vec3, b: Vec3) {
        out[0] = Math.max(a[0], b[0]);
        out[1] = Math.max(a[1], b[1]);
        out[2] = Math.max(a[2], b[2]);
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

    export function setMagnitude(out: Vec3, a: Vec3, l: number) {
        return Vec3.scale(out, Vec3.normalize(out, a), l);
    }

    /**
     * Negates the components of a vec3
     */
    export function negate(out: Vec3, a: Vec3) {
        out[0] = -a[0];
        out[1] = -a[1];
        out[2] = -a[2];
        return out;
    }

    /**
     * Returns the inverse of the components of a Vec3
     */
    export function inverse(out: Vec3, a: Vec3) {
        out[0] = 1.0 / a[0];
        out[1] = 1.0 / a[1];
        out[2] = 1.0 / a[2];
        return out;
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

    /**
     * Performs a linear interpolation between two Vec3's
     */
    export function lerp(out: Vec3, a: Vec3, b: Vec3, t: number) {
        const ax = a[0],
            ay = a[1],
            az = a[2];
        out[0] = ax + t * (b[0] - ax);
        out[1] = ay + t * (b[1] - ay);
        out[2] = az + t * (b[2] - az);
        return out;
    }

    const slerpRelVec = Vec3.zero();
    export function slerp(out: Vec3, a: Vec3, b: Vec3, t: number) {
        const dot = clamp(Vec3.dot(a, b), -1, 1);
        const theta = Math.acos(dot) * t;
        Vec3.scaleAndAdd(slerpRelVec, b, a, -dot);
        Vec3.normalize(slerpRelVec, slerpRelVec);
        return Vec3.add(out, Vec3.scale(out, a, Math.cos(theta)), Vec3.scale(slerpRelVec, slerpRelVec, Math.sin(theta)));
    }

    /**
     * Performs a hermite interpolation with two control points
     */
    export function hermite(out: Vec3, a: Vec3, b: Vec3, c: Vec3, d: Vec3, t: number) {
        const factorTimes2 = t * t;
        const factor1 = factorTimes2 * (2 * t - 3) + 1;
        const factor2 = factorTimes2 * (t - 2) + t;
        const factor3 = factorTimes2 * (t - 1);
        const factor4 = factorTimes2 * (3 - 2 * t);

        out[0] = a[0] * factor1 + b[0] * factor2 + c[0] * factor3 + d[0] * factor4;
        out[1] = a[1] * factor1 + b[1] * factor2 + c[1] * factor3 + d[1] * factor4;
        out[2] = a[2] * factor1 + b[2] * factor2 + c[2] * factor3 + d[2] * factor4;

        return out;
    }

    /**
     * Performs a bezier interpolation with two control points
     */
    export function bezier(out: Vec3, a: Vec3, b: Vec3, c: Vec3, d: Vec3, t: number) {
        const inverseFactor = 1 - t;
        const inverseFactorTimesTwo = inverseFactor * inverseFactor;
        const factorTimes2 = t * t;
        const factor1 = inverseFactorTimesTwo * inverseFactor;
        const factor2 = 3 * t * inverseFactorTimesTwo;
        const factor3 = 3 * factorTimes2 * inverseFactor;
        const factor4 = factorTimes2 * t;

        out[0] = a[0] * factor1 + b[0] * factor2 + c[0] * factor3 + d[0] * factor4;
        out[1] = a[1] * factor1 + b[1] * factor2 + c[1] * factor3 + d[1] * factor4;
        out[2] = a[2] * factor1 + b[2] * factor2 + c[2] * factor3 + d[2] * factor4;

        return out;
    }

    export function quadraticBezier(out: Vec3, a: Vec3, b: Vec3, c: Vec3, t: number) {
        out[0] = _quadraticBezier(a[0], b[0], c[0], t);
        out[1] = _quadraticBezier(a[1], b[1], c[1], t);
        out[2] = _quadraticBezier(a[2], b[2], c[2], t);

        return out;
    }

    /**
     * Performs a spline interpolation with two control points and a tension parameter
     */
    export function spline(out: Vec3, a: Vec3, b: Vec3, c: Vec3, d: Vec3, t: number, tension: number) {
        out[0] = _spline(a[0], b[0], c[0], d[0], t, tension);
        out[1] = _spline(a[1], b[1], c[1], d[1], t, tension);
        out[2] = _spline(a[2], b[2], c[2], d[2], t, tension);

        return out;
    }

    /**
     * Generates a random vector with the given scale
     */
    export function random(out: Vec3, scale: number) {
        const r = Math.random() * 2.0 * Math.PI;
        const z = (Math.random() * 2.0) - 1.0;
        const zScale = Math.sqrt(1.0 - z * z) * scale;

        out[0] = Math.cos(r) * zScale;
        out[1] = Math.sin(r) * zScale;
        out[2] = z * scale;
        return out;
    }

    /**
     * Transforms the Vec3 with a Mat4. 4th vector component is implicitly '1'
     */
    export function transformMat4(out: Vec3, a: Vec3, m: Mat4) {
        const x = a[0], y = a[1], z = a[2],
            w = 1 / ((m[3] * x + m[7] * y + m[11] * z + m[15]) || 1.0);
        out[0] = (m[0] * x + m[4] * y + m[8] * z + m[12]) * w;
        out[1] = (m[1] * x + m[5] * y + m[9] * z + m[13]) * w;
        out[2] = (m[2] * x + m[6] * y + m[10] * z + m[14]) * w;
        return out;
    }

    /**
     * Like `transformMat4` but with offsets into arrays
     */
    export function transformMat4Offset(out: NumberArray, a: NumberArray, m: NumberArray, outO: number, aO: number, oM: number) {
        const x = a[0 + aO], y = a[1 + aO], z = a[2 + aO],
            w = 1 / ((m[3 + oM] * x + m[7 + oM] * y + m[11 + oM] * z + m[15 + oM]) || 1.0);
        out[0 + outO] = (m[0 + oM] * x + m[4 + oM] * y + m[8 + oM] * z + m[12 + oM]) * w;
        out[1 + outO] = (m[1 + oM] * x + m[5 + oM] * y + m[9 + oM] * z + m[13 + oM]) * w;
        out[2 + outO] = (m[2 + oM] * x + m[6 + oM] * y + m[10 + oM] * z + m[14 + oM]) * w;
        return out;
    }

    /**
     * Transforms the Vec3 with a Mat3.
     */
    export function transformMat3(out: Vec3, a: Vec3, m: Mat3) {
        const x = a[0], y = a[1], z = a[2];
        out[0] = x * m[0] + y * m[3] + z * m[6];
        out[1] = x * m[1] + y * m[4] + z * m[7];
        out[2] = x * m[2] + y * m[5] + z * m[8];
        return out;
    }

    /** Transforms the Vec3 with a quat */
    export function transformQuat(out: Vec3, a: Vec3, q: Quat) {
        // benchmarks: http://jsperf.com/quaternion-transform-vec3-implementations

        const x = a[0], y = a[1], z = a[2];
        const qx = q[0], qy = q[1], qz = q[2], qw = q[3];

        // calculate quat * vec
        const ix = qw * x + qy * z - qz * y;
        const iy = qw * y + qz * x - qx * z;
        const iz = qw * z + qx * y - qy * x;
        const iw = -qx * x - qy * y - qz * z;

        // calculate result * inverse quat
        out[0] = ix * qw + iw * -qx + iy * -qz - iz * -qy;
        out[1] = iy * qw + iw * -qy + iz * -qx - ix * -qz;
        out[2] = iz * qw + iw * -qz + ix * -qy - iy * -qx;
        return out;
    }

    /** Computes the angle between 2 vectors, reports in radians. */
    export function angle(a: Vec3, b: Vec3) {
        const theta = dot(a, b) / Math.sqrt(squaredMagnitude(a) * squaredMagnitude(b));
        return Math.acos(clamp(theta, -1, 1)); // clamp to avoid numerical problems
    }

    const tmp_dh_ab = zero();
    const tmp_dh_cb = zero();
    const tmp_dh_bc = zero();
    const tmp_dh_dc = zero();
    const tmp_dh_abc = zero();
    const tmp_dh_bcd = zero();
    const tmp_dh_cross = zero();
    /**
     * Computes the dihedral angles of 4 points, reports in radians.
     */
    export function dihedralAngle(a: Vec3, b: Vec3, c: Vec3, d: Vec3): number {
        sub(tmp_dh_ab, a, b);
        sub(tmp_dh_cb, c, b);
        sub(tmp_dh_bc, b, c);
        sub(tmp_dh_dc, d, c);

        cross(tmp_dh_abc, tmp_dh_ab, tmp_dh_cb);
        cross(tmp_dh_bcd, tmp_dh_bc, tmp_dh_dc);

        const _angle = angle(tmp_dh_abc, tmp_dh_bcd);
        cross(tmp_dh_cross, tmp_dh_abc, tmp_dh_bcd);
        return dot(tmp_dh_cb, tmp_dh_cross) > 0 ? _angle : -_angle;
    }

    /**
     * Returns whether or not the vectors have exactly the same elements in the same position (when compared with ===)
     */
    export function exactEquals(a: Vec3, b: Vec3) {
        return a[0] === b[0] && a[1] === b[1] && a[2] === b[2];
    }

    /**
     * Returns whether or not the vectors have approximately the same elements in the same position.
     */
    export function equals(a: Vec3, b: Vec3) {
        const a0 = a[0], a1 = a[1], a2 = a[2];
        const b0 = b[0], b1 = b[1], b2 = b[2];
        return (Math.abs(a0 - b0) <= EPSILON * Math.max(1.0, Math.abs(a0), Math.abs(b0)) &&
                Math.abs(a1 - b1) <= EPSILON * Math.max(1.0, Math.abs(a1), Math.abs(b1)) &&
                Math.abs(a2 - b2) <= EPSILON * Math.max(1.0, Math.abs(a2), Math.abs(b2)));
    }

    const rotTemp = zero();
    export function makeRotation(mat: Mat4, a: Vec3, b: Vec3): Mat4 {
        const by = angle(a, b);
        if (Math.abs(by) < 0.0001) return Mat4.setIdentity(mat);
        if (Math.abs(by - Math.PI) < EPSILON) {
            // here, axis can be [0,0,0] but the rotation is a simple flip
            return Mat4.fromScaling(mat, negUnit);
        }
        const axis = cross(rotTemp, a, b);
        return Mat4.fromRotation(mat, by, axis);
    }

    export function isZero(v: Vec3) {
        return v[0] === 0 && v[1] === 0 && v[2] === 0;
    }

    /** Project `point` onto `vector` starting from `origin` */
    export function projectPointOnVector(out: Vec3, point: Vec3, vector: Vec3, origin: Vec3) {
        sub(out, copy(out, point), origin);
        const scalar = dot(vector, out) / squaredMagnitude(vector);
        return add(out, scale(out, copy(out, vector), scalar), origin);
    }

    export function projectOnVector(out: Vec3, p: Vec3, vector: Vec3 ) {
        const scalar = dot(vector, p) / squaredMagnitude(vector);
        return scale(out, vector, scalar);
    }

    const tmpProject = Vec3();
    export function projectOnPlane(out: Vec3, p: Vec3, normal: Vec3) {
        projectOnVector(tmpProject, p, normal);
        return sub(out, p, tmpProject);
    }

    /** Get a vector that is similar to `b` but orthogonal to `a` */
    export function orthogonalize(out: Vec3, a: Vec3, b: Vec3) {
        return normalize(out, cross(out, cross(out, a, b), a));
    }

    /**
     * Get a vector like `a` that point into the same general direction as `b`,
     * i.e. where the dot product is > 0
     */
    export function matchDirection(out: Vec3, a: Vec3, b: Vec3) {
        if (Vec3.dot(a, b) > 0) Vec3.copy(out, a);
        else Vec3.negate(out, Vec3.copy(out, a));
        return out;
    }

    const triangleNormalTmpAB = zero();
    const triangleNormalTmpAC = zero();
    /** Calculate normal for the triangle defined by `a`, `b` and `c` */
    export function triangleNormal(out: Vec3, a: Vec3, b: Vec3, c: Vec3) {
        sub(triangleNormalTmpAB, b, a);
        sub(triangleNormalTmpAC, c, a);
        return normalize(out, cross(out, triangleNormalTmpAB, triangleNormalTmpAC));
    }

    export function toString(a: Vec3, precision?: number) {
        return `[${a[0].toPrecision(precision)} ${a[1].toPrecision(precision)} ${a[2].toPrecision(precision)}]`;
    }

    export const origin: ReadonlyVec3 = Vec3.create(0, 0, 0);

    export const unit: ReadonlyVec3 = Vec3.create(1, 1, 1);
    export const negUnit: ReadonlyVec3 = Vec3.create(-1, -1, -1);

    export const unitX: ReadonlyVec3 = Vec3.create(1, 0, 0);
    export const unitY: ReadonlyVec3 = Vec3.create(0, 1, 0);
    export const unitZ: ReadonlyVec3 = Vec3.create(0, 0, 1);
}

export default Vec3;