/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
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

interface Vec2 extends Array<number> { [d: number]: number, '@type': 'vec2', length: 2 }

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

    export function toArray(a: Vec2, out: Helpers.NumberArray, offset: number) {
        out[offset + 0] = a[0];
        out[offset + 1] = a[1];
    }

    export function fromArray(a: Vec2, array: Helpers.NumberArray, offset: number) {
        a[0] = array[offset + 0]
        a[1] = array[offset + 1]
        return a
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
}

export default Vec2