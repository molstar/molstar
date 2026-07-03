/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as B from 'benchmark';
import { fft1d, fft3d, fft3dHermitian } from '../mol-math/fourier-transform/fft';

// Compact radix-2 (Float64) baseline, kept here only for a head-to-head 1D
// comparison against the radix-4 fft1d. Each benched call copies fresh input
// first (fft is in-place); the O(n) copy is the same for both variants.
function fft1dRadix2(re: Float64Array, im: Float64Array, n: number) {
    for (let i = 1, j = 0; i < n; i++) {
        let bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            const tr = re[i]; re[i] = re[j]; re[j] = tr;
            const ti = im[i]; im[i] = im[j]; im[j] = ti;
        }
    }
    for (let len = 2; len <= n; len <<= 1) {
        const half = len >> 1;
        const ang = -2 * Math.PI / len;
        const wr = Math.cos(ang), wi = Math.sin(ang);
        for (let s = 0; s < n; s += len) {
            let cr = 1, ci = 0;
            for (let k = 0; k < half; k++) {
                const a = s + k, b = a + half;
                const tr = re[b] * cr - im[b] * ci;
                const ti = re[b] * ci + im[b] * cr;
                re[b] = re[a] - tr; im[b] = im[a] - ti;
                re[a] += tr; im[a] += ti;
                const ncr = cr * wr - ci * wi;
                ci = cr * wi + ci * wr; cr = ncr;
            }
        }
    }
}

function randomData(n: number) {
    const re = new Float64Array(n), im = new Float64Array(n);
    for (let i = 0; i < n; i++) { re[i] = Math.random() * 2 - 1; im[i] = Math.random() * 2 - 1; }
    return { re, im };
}

export function run1d(n: number) {
    const { re, im } = randomData(n);
    const wr = new Float64Array(n), wi = new Float64Array(n);
    new B.Suite()
        .add(`radix-2  n=${n}`, () => { wr.set(re); wi.set(im); fft1dRadix2(wr, wi, n); })
        .add(`radix-4  n=${n}`, () => { wr.set(re); wi.set(im); fft1d(wr, wi, 0, n, 1); })
        .on('cycle', (e: any) => console.log(String(e.target)))
        .run();
    console.log('---------------------');
}

export function run3d(N0: number, N1: number, N2: number) {
    const total = N0 * N1 * N2;
    const { re, im } = randomData(total);
    const wr = new Float32Array(total), wi = new Float32Array(total);
    new B.Suite()
        .add(`fft3d ${N0}x${N1}x${N2}`, () => { wr.set(re); wi.set(im); fft3d(wr, wi, N0, N1, N2); })
        .on('cycle', (e: any) => console.log(String(e.target)))
        .run();
    console.log('---------------------');
}

// Random Hermitian grid: F[-h,-k,-l] = conj(F[h,k,l]); self-conjugate points real.
function randomHermitian(N0: number, N1: number, N2: number) {
    const total = N0 * N1 * N2;
    const re = new Float32Array(total), im = new Float32Array(total);
    const wrap = (x: number, N: number) => ((x % N) + N) % N;
    for (let h = 0; h < N0; h++) {
        for (let k = 0; k < N1; k++) {
            for (let l = 0; l < N2; l++) {
                const i = h * N1 * N2 + k * N2 + l;
                const j = wrap(-h, N0) * N1 * N2 + wrap(-k, N1) * N2 + wrap(-l, N2);
                if (i < j) {
                    const r = Math.random() * 2 - 1, m = Math.random() * 2 - 1;
                    re[i] = r; im[i] = m; re[j] = r; im[j] = -m;
                } else if (i === j) {
                    re[i] = Math.random() * 2 - 1; im[i] = 0;
                }
            }
        }
    }
    return { re, im };
}

// Hermitian input + real output: fft3dHermitian vs the general complex fft3d.
export function run3dHermitian(N0: number, N1: number, N2: number) {
    const total = N0 * N1 * N2;
    const P = (N2 >> 1) + 1;
    const half = N0 * N1 * P;
    const { re, im } = randomHermitian(N0, N1, N2);
    // Compact lower half-spectrum (planes l = 0..N2/2) for the specialized transform.
    const hr = new Float32Array(half), hi = new Float32Array(half);
    for (let h = 0; h < N0; h++) {
        for (let k = 0; k < N1; k++) {
            for (let l = 0; l < P; l++) {
                const s = h * N1 * N2 + k * N2 + l;
                const d = h * N1 * P + k * P + l;
                hr[d] = re[s]; hi[d] = im[s];
            }
        }
    }
    const wr = new Float32Array(total), wi = new Float32Array(total);
    const whr = new Float32Array(half), whi = new Float32Array(half), out = new Float32Array(total);
    new B.Suite()
        .add(`fft3d          ${N0}x${N1}x${N2}`, () => { wr.set(re); wi.set(im); fft3d(wr, wi, N0, N1, N2); })
        .add(`fft3dHermitian ${N0}x${N1}x${N2}`, () => { whr.set(hr); whi.set(hi); fft3dHermitian(whr, whi, N0, N1, N2, out); })
        .on('cycle', (e: any) => console.log(String(e.target)))
        .run();
    console.log('---------------------');
}

run1d(1024);
run1d(4096);
run1d(8192);
run3d(64, 64, 64);
run3d(128, 128, 128);
run3dHermitian(64, 64, 64);
run3dHermitian(128, 128, 128);

// radix-2  n=1024 x 91,184 ops/sec ±0.71% (96 runs sampled)
// radix-4  n=1024 x 99,081 ops/sec ±0.38% (97 runs sampled)
// ---------------------
// radix-2  n=4096 x 14,877 ops/sec ±0.33% (99 runs sampled)
// radix-4  n=4096 x 20,120 ops/sec ±0.29% (98 runs sampled)
// ---------------------
// radix-2  n=8192 x 8,182 ops/sec ±2.75% (92 runs sampled)
// radix-4  n=8192 x 9,683 ops/sec ±0.34% (100 runs sampled)
// ---------------------
// fft3d 64x64x64 x 82.82 ops/sec ±0.23% (72 runs sampled)
// ---------------------
// fft3d 128x128x128 x 9.37 ops/sec ±0.36% (28 runs sampled)
// ---------------------
// fft3d          64x64x64 x 84.76 ops/sec ±0.23% (74 runs sampled)
// fft3dHermitian 64x64x64 x 169 ops/sec ±0.32% (87 runs sampled)
// ---------------------
// fft3d          128x128x128 x 9.54 ops/sec ±0.30% (28 runs sampled)
// fft3dHermitian 128x128x128 x 18.82 ops/sec ±0.53% (51 runs sampled)
