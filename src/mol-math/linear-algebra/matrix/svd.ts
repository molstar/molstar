/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Matrix from './matrix';
import { NumberArray } from '../../../mol-util/type-helpers';

// svd method adapted from http://inspirit.github.io/jsfeat/ MIT Eugene Zatepyakin

export function swap(A: NumberArray, i0: number, i1: number, t: number) {
    t = A[i0];
    A[i0] = A[i1];
    A[i1] = t;
}

export function hypot(a: number, b: number) {
    a = Math.abs(a);
    b = Math.abs(b);
    if (a > b) {
        b /= a;
        return a * Math.sqrt(1.0 + b * b);
    }
    if (b > 0) {
        a /= b;
        return b * Math.sqrt(1.0 + a * a);
    }
    return 0.0;
}

const EPSILON = 0.0000001192092896;
const FLT_MIN = 1E-37;

export function JacobiSVDImpl(At: NumberArray, astep: number, _W: NumberArray, Vt: NumberArray, vstep: number, m: number, n: number, n1: number) {
    const eps = EPSILON * 2.0;
    const minval = FLT_MIN;
    let i = 0;
    let j = 0;
    let k = 0;
    let iter = 0;
    const maxIter = Math.max(m, 30);
    let Ai = 0;
    let Aj = 0;
    let Vi = 0;
    let Vj = 0;
    let changed = 0;
    let c = 0.0;
    let s = 0.0;
    let t = 0.0;
    let t0 = 0.0;
    let t1 = 0.0;
    let sd = 0.0;
    let beta = 0.0;
    let gamma = 0.0;
    let delta = 0.0;
    let a = 0.0;
    let p = 0.0;
    let b = 0.0;
    let seed = 0x1234;
    let val = 0.0;
    let val0 = 0.0;
    let asum = 0.0;

    const W = new Float64Array(n << 3);

    for (; i < n; i++) {
        for (k = 0, sd = 0; k < m; k++) {
            t = At[i * astep + k];
            sd += t * t;
        }
        W[i] = sd;

        if (Vt) {
            for (k = 0; k < n; k++) {
                Vt[i * vstep + k] = 0;
            }
            Vt[i * vstep + i] = 1;
        }
    }

    for (; iter < maxIter; iter++) {
        changed = 0;

        for (i = 0; i < n - 1; i++) {
            for (j = i + 1; j < n; j++) {
                Ai = (i * astep) | 0;
                Aj = (j * astep) | 0;
                a = W[i];
                p = 0;
                b = W[j];

                k = 2;
                p += At[Ai] * At[Aj];
                p += At[Ai + 1] * At[Aj + 1];

                for (; k < m; k++) { p += At[Ai + k] * At[Aj + k]; }

                if (Math.abs(p) <= eps * Math.sqrt(a * b)) continue;

                p *= 2.0;
                beta = a - b;
                gamma = hypot(p, beta);
                if (beta < 0) {
                    delta = (gamma - beta) * 0.5;
                    s = Math.sqrt(delta / gamma);
                    c = (p / (gamma * s * 2.0));
                } else {
                    c = Math.sqrt((gamma + beta) / (gamma * 2.0));
                    s = (p / (gamma * c * 2.0));
                }

                a = 0.0;
                b = 0.0;

                k = 2; // unroll
                t0 = c * At[Ai] + s * At[Aj];
                t1 = -s * At[Ai] + c * At[Aj];
                At[Ai] = t0; At[Aj] = t1;
                a += t0 * t0; b += t1 * t1;

                t0 = c * At[Ai + 1] + s * At[Aj + 1];
                t1 = -s * At[Ai + 1] + c * At[Aj + 1];
                At[Ai + 1] = t0; At[Aj + 1] = t1;
                a += t0 * t0; b += t1 * t1;

                for (; k < m; k++) {
                    t0 = c * At[Ai + k] + s * At[Aj + k];
                    t1 = -s * At[Ai + k] + c * At[Aj + k];
                    At[Ai + k] = t0; At[Aj + k] = t1;

                    a += t0 * t0; b += t1 * t1;
                }

                W[i] = a;
                W[j] = b;

                changed = 1;

                if (Vt) {
                    Vi = (i * vstep) | 0;
                    Vj = (j * vstep) | 0;

                    k = 2;
                    t0 = c * Vt[Vi] + s * Vt[Vj];
                    t1 = -s * Vt[Vi] + c * Vt[Vj];
                    Vt[Vi] = t0; Vt[Vj] = t1;

                    t0 = c * Vt[Vi + 1] + s * Vt[Vj + 1];
                    t1 = -s * Vt[Vi + 1] + c * Vt[Vj + 1];
                    Vt[Vi + 1] = t0; Vt[Vj + 1] = t1;

                    for (; k < n; k++) {
                        t0 = c * Vt[Vi + k] + s * Vt[Vj + k];
                        t1 = -s * Vt[Vi + k] + c * Vt[Vj + k];
                        Vt[Vi + k] = t0; Vt[Vj + k] = t1;
                    }
                }
            }
        }
        if (changed === 0) break;
    }

    for (i = 0; i < n; i++) {
        for (k = 0, sd = 0; k < m; k++) {
            t = At[i * astep + k];
            sd += t * t;
        }
        W[i] = Math.sqrt(sd);
    }

    for (i = 0; i < n - 1; i++) {
        j = i;
        for (k = i + 1; k < n; k++) {
            if (W[j] < W[k]) { j = k; }
        }
        if (i !== j) {
            swap(W, i, j, sd);
            if (Vt) {
                for (k = 0; k < m; k++) {
                    swap(At, i * astep + k, j * astep + k, t);
                }

                for (k = 0; k < n; k++) {
                    swap(Vt, i * vstep + k, j * vstep + k, t);
                }
            }
        }
    }

    for (i = 0; i < n; i++) {
        _W[i] = W[i];
    }

    if (!Vt) {
        return;
    }

    for (i = 0; i < n1; i++) {
        sd = i < n ? W[i] : 0;

        while (sd <= minval) {
            // if we got a zero singular value, then in order to get the corresponding left singular vector
            // we generate a random vector, project it to the previously computed left singular vectors,
            // subtract the projection and normalize the difference.
            val0 = (1.0 / m);
            for (k = 0; k < m; k++) {
                seed = (seed * 214013 + 2531011);
                val = (((seed >> 16) & 0x7fff) & 256) !== 0 ? val0 : -val0;
                At[i * astep + k] = val;
            }
            for (iter = 0; iter < 2; iter++) {
                for (j = 0; j < i; j++) {
                    sd = 0;
                    for (k = 0; k < m; k++) {
                        sd += At[i * astep + k] * At[j * astep + k];
                    }
                    asum = 0.0;
                    for (k = 0; k < m; k++) {
                        t = (At[i * astep + k] - sd * At[j * astep + k]);
                        At[i * astep + k] = t;
                        asum += Math.abs(t);
                    }
                    asum = asum ? 1.0 / asum : 0;
                    for (k = 0; k < m; k++) {
                        At[i * astep + k] *= asum;
                    }
                }
            }
            sd = 0;
            for (k = 0; k < m; k++) {
                t = At[i * astep + k];
                sd += t * t;
            }
            sd = Math.sqrt(sd);
        }

        s = (1.0 / sd);
        for (k = 0; k < m; k++) {
            At[i * astep + k] *= s;
        }
    }
}

export function svd(A: Matrix, W: Matrix, U: Matrix, V: Matrix) {
    let at = 0;
    let i = 0;
    const _m = A.rows;
    const _n = A.cols;
    let m = _m;
    let n = _n;

    if (m < n) {
        at = 1;
        i = m;
        m = n;
        n = i;
    }

    const amt = Matrix.create(m, m);
    const wmt = Matrix.create(1, n);
    const vmt = Matrix.create(n, n);

    if (at === 0) {
        Matrix.transpose(amt, A);
    } else {
        for (i = 0; i < _n * _m; i++) {
            amt.data[i] = A.data[i];
        }
        for (; i < n * m; i++) {
            amt.data[i] = 0;
        }
    }

    JacobiSVDImpl(amt.data, m, wmt.data, vmt.data, n, m, n, m);

    if (W) {
        for (i = 0; i < n; i++) {
            W.data[i] = wmt.data[i];
        }
        for (; i < _n; i++) {
            W.data[i] = 0;
        }
    }

    if (at === 0) {
        if (U) Matrix.transpose(U, amt);
        if (V) Matrix.transpose(V, vmt);
    } else {
        if (U) Matrix.transpose(U, vmt);
        if (V) Matrix.transpose(V, amt);
    }
}