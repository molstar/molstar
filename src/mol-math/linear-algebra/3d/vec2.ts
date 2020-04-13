/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from '../../../mol-util/type-helpers';

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

interface Vec2 extends Array<number> { [d: number]: number, '@type': 'vec2', length: 2 }

function Vec2() {
    return Vec2.zero();
}

namespace Vec2 {
    export function zero(): Vec2 {
        // force double backing array by 0.1.
        const ret = [0.1, 0];
        ret[0] = 0.0;
        return ret as any;
    }

    export function clone(a: Vec2) {
        const out = zero();
        out[0] = a[0];
        out[1] = a[1];
        return out;
    }

    export function create(x: number, y: number) {
        const out = zero();
        out[0] = x;
        out[1] = y;
        return out;
    }

    export function hasNaN(a: Vec2) {
        return isNaN(a[0]) || isNaN(a[1]);
    }

    export function toArray(a: Vec2, out: NumberArray, offset: number) {
        out[offset + 0] = a[0];
        out[offset + 1] = a[1];
        return out;
    }

    export function fromArray(a: Vec2, array: NumberArray, offset: number) {
        a[0] = array[offset + 0];
        a[1] = array[offset + 1];
        return a;
    }

    export function copy(out: Vec2, a: Vec2) {
        out[0] = a[0];
        out[1] = a[1];
        return out;
    }

    export function set(out: Vec2, x: number, y: number) {
        out[0] = x;
        out[1] = y;
        return out;
    }

    export function add(out: Vec2, a: Vec2, b: Vec2) {
        out[0] = a[0] + b[0];
        out[1] = a[1] + b[1];
        return out;
    }

    export function sub(out: Vec2, a: Vec2, b: Vec2) {
        out[0] = a[0] - b[0];
        out[1] = a[1] - b[1];
        return out;
    }

    export function mul(out: Vec2, a: Vec2, b: Vec2) {
        out[0] = a[0] * b[0];
        out[1] = a[1] * b[1];
        return out;
    }

    export function div(out: Vec2, a: Vec2, b: Vec2) {
        out[0] = a[0] / b[0];
        out[1] = a[1] / b[1];
        return out;
    }

    export function scale(out: Vec2, a: Vec2, b: number) {
        out[0] = a[0] * b;
        out[1] = a[1] * b;
        return out;
    }

    /**
     * Math.round the components of a Vec2
     */
    export function round(out: Vec2, a: Vec2) {
        out[0] = Math.round(a[0]);
        out[1] = Math.round(a[1]);
        return out;
    }

    /**
     * Math.ceil the components of a Vec2
     */
    export function ceil(out: Vec2, a: Vec2) {
        out[0] = Math.ceil(a[0]);
        out[1] = Math.ceil(a[1]);
        return out;
    }

    /**
     * Math.floor the components of a Vec2
     */
    export function floor(out: Vec2, a: Vec2) {
        out[0] = Math.floor(a[0]);
        out[1] = Math.floor(a[1]);
        return out;
    }

    export function distance(a: Vec2, b: Vec2) {
        const x = b[0] - a[0],
            y = b[1] - a[1];
        return Math.sqrt(x * x + y * y);
    }

    export function squaredDistance(a: Vec2, b: Vec2) {
        const x = b[0] - a[0],
            y = b[1] - a[1];
        return x * x + y * y;
    }

    export function magnitude(a: Vec2) {
        const x = a[0],
            y = a[1];
        return Math.sqrt(x * x + y * y);
    }

    export function squaredMagnitude(a: Vec2) {
        const x = a[0],
            y = a[1];
        return x * x + y * y;
    }

    /**
     * Returns the inverse of the components of a Vec2
     */
    export function inverse(out: Vec2, a: Vec2) {
        out[0] = 1.0 / a[0];
        out[1] = 1.0 / a[1];
        return out;
    }

    export function areEqual(a: Vec2, b: Vec2) {
        return a[0] === b[0] && a[1] === b[1];
    }

    export function toString(a: Vec2, precision?: number) {
        return `[${a[0].toPrecision(precision)} ${a[1].toPrecision(precision)}}]`;
    }
}

export default Vec2;