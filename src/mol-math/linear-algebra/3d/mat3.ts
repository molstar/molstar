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

interface Mat3 extends Array<number> { [d: number]: number, '@type': 'mat3', length: 9 }

namespace Mat3 {
    export function zero(): Mat3 {
        // force double backing array by 0.1.
        const ret = [0.1, 0, 0, 0, 0, 0, 0, 0, 0];
        ret[0] = 0.0;
        return ret as any;
    }

    export function toArray(a: Mat3, out: Helpers.NumberArray, offset: number) {
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
    }

    export function fromArray(a: Mat3, array: Helpers.NumberArray, offset: number) {
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
        return a
    }
}

export default Mat3