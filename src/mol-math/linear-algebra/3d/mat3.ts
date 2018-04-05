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
}

export default Mat3