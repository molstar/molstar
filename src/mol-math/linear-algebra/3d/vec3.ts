/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { Quat } from '../3d';

interface Vec3 extends Array<number> { [d: number]: number, '@type': 'vec3', length: 3 }

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
        return v
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

    /**
     * Generates a random vector with the given scale
     */
    export function random(out: Vec3, scale: number) {
        const r = Math.random() * 2.0 * Math.PI;
        const z = (Math.random() * 2.0) - 1.0;
        const zScale = Math.sqrt(1.0-z*z) * scale;

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

    /** Transforms the vec3 with a quat */
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

    export function isZero(v: Vec3) {
        return v[0] === 0 && v[1] === 0 && v[2] === 0
    }
}

export default Vec3