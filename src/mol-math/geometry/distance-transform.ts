/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Tadej Satler <tadej.satler@gmail.com>
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

const Inf = 1e20;

/**
 * Squared 3D Euclidean distance from every cell to the nearest cell with `inside[i] !== 0`,
 * in grid units (Felzenszwalb & Huttenlocher, one separable pass per axis). `inside` is in
 * x-fastest order (`x + y * nx + z * nx * ny`). Cells inside get 0; if nothing is inside,
 * all values are `>= 1e20`.
 */
export function squaredDistanceTransform3D(inside: ArrayLike<number>, nx: number, ny: number, nz: number, out?: Float32Array): Float32Array {
    const n = nx * ny * nz;
    const dist = out && out.length === n ? out : new Float32Array(n);
    for (let i = 0; i < n; i++) dist[i] = inside[i] ? 0 : Inf;

    const maxDim = Math.max(nx, ny, nz);
    const f = new Float64Array(maxDim);
    const d = new Float64Array(maxDim);
    const v = new Int32Array(maxDim);
    const z = new Float64Array(maxDim + 1);
    const nxy = nx * ny;

    for (let k = 0; k < nz; k++) {
        for (let j = 0; j < ny; j++) {
            const base = j * nx + k * nxy;
            for (let i = 0; i < nx; i++) f[i] = dist[base + i];
            edt1d(f, d, v, z, nx);
            for (let i = 0; i < nx; i++) dist[base + i] = d[i];
        }
    }
    for (let k = 0; k < nz; k++) {
        for (let i = 0; i < nx; i++) {
            const base = i + k * nxy;
            for (let j = 0; j < ny; j++) f[j] = dist[base + j * nx];
            edt1d(f, d, v, z, ny);
            for (let j = 0; j < ny; j++) dist[base + j * nx] = d[j];
        }
    }
    for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
            const base = i + j * nx;
            for (let k = 0; k < nz; k++) f[k] = dist[base + k * nxy];
            edt1d(f, d, v, z, nz);
            for (let k = 0; k < nz; k++) dist[base + k * nxy] = d[k];
        }
    }
    return dist;
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
