/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from '../../mol-util/type-helpers';

/**
 * 2D Euclidean distance transform by Felzenszwalb & Huttenlocher https://cs.brown.edu/~pff/papers/dt-final.pdf
 */
export function edt(data: NumberArray, width: number, height: number, f: NumberArray, d: NumberArray, v: NumberArray, z: NumberArray) {
    for (let x = 0; x < width; x++) {
        for (let y = 0; y < height; y++) {
            f[y] = data[y * width + x];
        }
        edt1d(f, d, v, z, height);
        for (let y = 0; y < height; y++) {
            data[y * width + x] = d[y];
        }
    }
    for (let y = 0; y < height; y++) {
        for (let x = 0; x < width; x++) {
            f[x] = data[y * width + x];
        }
        edt1d(f, d, v, z, width);
        for (let x = 0; x < width; x++) {
            data[y * width + x] = Math.sqrt(d[x]);
        }
    }
}

/**
 * 1D squared distance transform
 */
function edt1d(f: NumberArray, d: NumberArray, v: NumberArray, z: NumberArray, n: number) {
    v[0] = 0;
    z[0] = Number.MIN_SAFE_INTEGER;
    z[1] = Number.MAX_SAFE_INTEGER;

    for (let q = 1, k = 0; q < n; q++) {
        let s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
        while (s <= z[k]) {
            k--;
            s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
        }
        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = Number.MAX_SAFE_INTEGER;
    }

    for (let q = 0, k = 0; q < n; q++) {
        while (z[k + 1] < q) k++;
        d[q] = (q - v[k]) * (q - v[k]) + f[v[k]];
    }
}
