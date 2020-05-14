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
import { EPSILON } from '../3d';
import { NumberArray } from '../../../mol-util/type-helpers';
import { Sphere3D } from '../../geometry/primitives/sphere3d';

interface Vec4 extends Array<number> { [d: number]: number, '@type': 'vec4', length: 4 }

function Vec4() {
    return Vec4.zero();
}

namespace Vec4 {
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

    export function fromSphere(out: Vec4, sphere: Sphere3D) {
        out[0] = sphere.center[0];
        out[1] = sphere.center[1];
        out[2] = sphere.center[2];
        out[3] = sphere.radius;
        return out;
    }

    export function ofSphere(sphere: Sphere3D) {
        return fromSphere(zero(), sphere);
    }

    export function hasNaN(a: Vec4) {
        return isNaN(a[0]) || isNaN(a[1]) || isNaN(a[2]) || isNaN(a[3]);
    }

    export function toArray(a: Vec4, out: NumberArray, offset: number) {
        out[offset + 0] = a[0];
        out[offset + 1] = a[1];
        out[offset + 2] = a[2];
        out[offset + 3] = a[3];
        return out;
    }

    export function fromArray(a: Vec4, array: NumberArray, offset: number) {
        a[0] = array[offset + 0];
        a[1] = array[offset + 1];
        a[2] = array[offset + 2];
        a[3] = array[offset + 3];
        return a;
    }

    export function toVec3Array(a: Vec4, out: NumberArray, offset: number) {
        out[offset + 0] = a[0];
        out[offset + 1] = a[1];
        out[offset + 2] = a[2];
    }

    export function fromVec3Array(a: Vec4, array: NumberArray, offset: number) {
        a[0] = array[offset + 0];
        a[1] = array[offset + 1];
        a[2] = array[offset + 2];
        a[3] = 0;
        return a;
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

    export function add(out: Vec4, a: Vec4, b: Vec4) {
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

    export function scale(out: Vec4, a: Vec4, b: number) {
        out[0] = a[0] * b;
        out[1] = a[1] * b;
        out[2] = a[2] * b;
        out[4] = a[4] * b;
        return out;
    }

    /**
     * Math.round the components of a Vec4
     */
    export function round(out: Vec4, a: Vec4) {
        out[0] = Math.round(a[0]);
        out[1] = Math.round(a[1]);
        out[2] = Math.round(a[2]);
        out[3] = Math.round(a[3]);
        return out;
    }

    /**
     * Math.ceil the components of a Vec4
     */
    export function ceil(out: Vec4, a: Vec4) {
        out[0] = Math.ceil(a[0]);
        out[1] = Math.ceil(a[1]);
        out[2] = Math.ceil(a[2]);
        out[3] = Math.ceil(a[3]);
        return out;
    }

    /**
     * Math.floor the components of a Vec3
     */
    export function floor(out: Vec4, a: Vec4) {
        out[0] = Math.floor(a[0]);
        out[1] = Math.floor(a[1]);
        out[2] = Math.floor(a[2]);
        out[3] = Math.floor(a[3]);
        return out;
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

    export function transformMat4(out: Vec4, a: Vec4, m: Mat4) {
        const x = a[0], y = a[1], z = a[2], w = a[3];
        out[0] = m[0] * x + m[4] * y + m[8] * z + m[12] * w;
        out[1] = m[1] * x + m[5] * y + m[9] * z + m[13] * w;
        out[2] = m[2] * x + m[6] * y + m[10] * z + m[14] * w;
        out[3] = m[3] * x + m[7] * y + m[11] * z + m[15] * w;
        return out;
    }

    export function dot(a: Vec4, b: Vec4) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    }

    /**
     * Returns the inverse of the components of a Vec4
     */
    export function inverse(out: Vec4, a: Vec4) {
        out[0] = 1.0 / a[0];
        out[1] = 1.0 / a[1];
        out[2] = 1.0 / a[2];
        out[3] = 1.0 / a[3];
        return out;
    }

    /**
     * Returns whether or not the vectors have exactly the same elements in the same position (when compared with ===)
     */
    export function exactEquals(a: Vec4, b: Vec4) {
        return a[0] === b[0] && a[1] === b[1] && a[2] === b[2] && a[3] === b[3];
    }

    /**
     * Returns whether or not the vectors have approximately the same elements in the same position.
     */
    export function equals(a: Vec4, b: Vec4) {
        const a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
        const b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
        return (Math.abs(a0 - b0) <= EPSILON * Math.max(1.0, Math.abs(a0), Math.abs(b0)) &&
                Math.abs(a1 - b1) <= EPSILON * Math.max(1.0, Math.abs(a1), Math.abs(b1)) &&
                Math.abs(a2 - b2) <= EPSILON * Math.max(1.0, Math.abs(a2), Math.abs(b2)) &&
                Math.abs(a3 - b3) <= EPSILON * Math.max(1.0, Math.abs(a3), Math.abs(b3)));
    }

    export function toString(a: Vec4, precision?: number) {
        return `[${a[0].toPrecision(precision)} ${a[1].toPrecision(precision)} ${a[2].toPrecision(precision)}  ${a[3].toPrecision(precision)}]`;
    }
}

export default Vec4;