/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { fft1d, fft3d, fft3dHermitian, nextPow2, FftArray } from '../fft';

// Small deterministic PRNG so the tests are reproducible.
function mulberry32(seed: number) {
    return () => {
        seed |= 0; seed = (seed + 0x6D2B79F5) | 0;
        let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

function makeArray(Ctor: Float32ArrayConstructor | Float64ArrayConstructor, data: number[]): FftArray {
    const a = new Ctor(data.length);
    a.set(data);
    return a;
}

// Naive O(n^2) reference DFT (forward, -2πi sign).
function naiveDft1d(re: number[], im: number[]) {
    const n = re.length;
    const outRe = new Float64Array(n);
    const outIm = new Float64Array(n);
    for (let k = 0; k < n; k++) {
        let sRe = 0, sIm = 0;
        for (let j = 0; j < n; j++) {
            const angle = -2 * Math.PI * j * k / n;
            const c = Math.cos(angle), s = Math.sin(angle);
            sRe += re[j] * c - im[j] * s;
            sIm += re[j] * s + im[j] * c;
        }
        outRe[k] = sRe; outIm[k] = sIm;
    }
    return { re: outRe, im: outIm };
}

// Naive O((N0*N1*N2)^2) reference 3D DFT (forward, -2πi sign).
function naiveDft3d(re: number[], im: number[], N0: number, N1: number, N2: number) {
    const total = N0 * N1 * N2;
    const outRe = new Float64Array(total);
    const outIm = new Float64Array(total);
    for (let a0 = 0; a0 < N0; a0++) {
        for (let a1 = 0; a1 < N1; a1++) {
            for (let a2 = 0; a2 < N2; a2++) {
                let sRe = 0, sIm = 0;
                for (let h = 0; h < N0; h++) {
                    for (let k = 0; k < N1; k++) {
                        for (let l = 0; l < N2; l++) {
                            const angle = -2 * Math.PI * (h * a0 / N0 + k * a1 / N1 + l * a2 / N2);
                            const c = Math.cos(angle), s = Math.sin(angle);
                            const idx = h * N1 * N2 + k * N2 + l;
                            sRe += re[idx] * c - im[idx] * s;
                            sIm += re[idx] * s + im[idx] * c;
                        }
                    }
                }
                const o = a0 * N1 * N2 + a1 * N2 + a2;
                outRe[o] = sRe; outIm[o] = sIm;
            }
        }
    }
    return { re: outRe, im: outIm };
}

describe('fft1d', () => {
    // Includes powers of 4 (4, 16, 64) and non-powers of 4 (2, 8, 32, 128)
    // to exercise both the pure radix-4 path and the leading radix-2 stage.
    const sizes = [1, 2, 4, 8, 16, 32, 64, 128];

    for (const Ctor of [Float64Array, Float32Array] as const) {
        const isF32 = Ctor === Float32Array;
        for (const n of sizes) {
            it(`matches naive DFT, n=${n} (${Ctor.name})`, () => {
                const rng = mulberry32(1234 + n);
                const re: number[] = [], im: number[] = [];
                for (let i = 0; i < n; i++) { re.push(rng() * 2 - 1); im.push(rng() * 2 - 1); }
                const ref = naiveDft1d(re, im);

                const fr = makeArray(Ctor, re);
                const fi = makeArray(Ctor, im);
                fft1d(fr, fi, 0, n, 1);

                const tol = (isF32 ? 1e-3 : 1e-9) * n;
                let maxErr = 0;
                for (let k = 0; k < n; k++) {
                    maxErr = Math.max(maxErr, Math.abs(fr[k] - ref.re[k]), Math.abs(fi[k] - ref.im[k]));
                }
                expect(maxErr).toBeLessThan(tol);
            });
        }
    }

    it('respects offset and stride', () => {
        const n = 16, stride = 3, offset = 5, total = offset + n * stride + 2;
        const rng = mulberry32(99);
        const re: number[] = [], im: number[] = [];
        for (let i = 0; i < n; i++) { re.push(rng() * 2 - 1); im.push(rng() * 2 - 1); }
        const ref = naiveDft1d(re, im);

        const fr = new Float64Array(total);
        const fi = new Float64Array(total);
        for (let i = 0; i < n; i++) { fr[offset + i * stride] = re[i]; fi[offset + i * stride] = im[i]; }
        fft1d(fr, fi, offset, n, stride);

        let maxErr = 0;
        for (let k = 0; k < n; k++) {
            maxErr = Math.max(maxErr, Math.abs(fr[offset + k * stride] - ref.re[k]), Math.abs(fi[offset + k * stride] - ref.im[k]));
        }
        expect(maxErr).toBeLessThan(1e-9);
        // Positions outside the strided sequence must be untouched.
        expect(fr[0]).toBe(0);
        expect(fr[offset + 1]).toBe(0);
    });

    it('transforms a unit impulse to a flat spectrum', () => {
        const n = 32;
        const re = new Float64Array(n), im = new Float64Array(n);
        re[0] = 1;
        fft1d(re, im, 0, n, 1);
        for (let k = 0; k < n; k++) {
            expect(re[k]).toBeCloseTo(1, 10);
            expect(im[k]).toBeCloseTo(0, 10);
        }
    });
});

describe('fft3d', () => {
    const grids: [number, number, number][] = [[2, 4, 8], [4, 8, 16], [8, 8, 8]];

    for (const Ctor of [Float64Array, Float32Array] as const) {
        const isF32 = Ctor === Float32Array;
        for (const [N0, N1, N2] of grids) {
            it(`matches naive 3D DFT, ${N0}x${N1}x${N2} (${Ctor.name})`, () => {
                const total = N0 * N1 * N2;
                const rng = mulberry32(777 + total);
                const re: number[] = [], im: number[] = [];
                for (let i = 0; i < total; i++) { re.push(rng() * 2 - 1); im.push(rng() * 2 - 1); }
                const ref = naiveDft3d(re, im, N0, N1, N2);

                const fr = makeArray(Ctor, re);
                const fi = makeArray(Ctor, im);
                fft3d(fr, fi, N0, N1, N2);

                const tol = isF32 ? 1e-2 : 1e-7;
                let maxErr = 0;
                for (let i = 0; i < total; i++) {
                    maxErr = Math.max(maxErr, Math.abs(fr[i] - ref.re[i]), Math.abs(fi[i] - ref.im[i]));
                }
                expect(maxErr).toBeLessThan(tol);
            });
        }
    }
});

// Build a random Hermitian-symmetric grid: F[-h,-k,-l] = conj(F[h,k,l]).
// Self-conjugate points (where (-h,-k,-l) wraps onto (h,k,l)) are forced real.
function makeHermitianGrid(N0: number, N1: number, N2: number, seed: number) {
    const total = N0 * N1 * N2;
    const re = new Array<number>(total).fill(0);
    const im = new Array<number>(total).fill(0);
    const rng = mulberry32(seed);
    const wrap = (x: number, N: number) => ((x % N) + N) % N;
    for (let h = 0; h < N0; h++) {
        for (let k = 0; k < N1; k++) {
            for (let l = 0; l < N2; l++) {
                const i = h * N1 * N2 + k * N2 + l;
                const j = wrap(-h, N0) * N1 * N2 + wrap(-k, N1) * N2 + wrap(-l, N2);
                if (i < j) {
                    const r = rng() * 2 - 1, m = rng() * 2 - 1;
                    re[i] = r; im[i] = m;
                    re[j] = r; im[j] = -m;
                } else if (i === j) {
                    re[i] = rng() * 2 - 1; im[i] = 0;
                }
            }
        }
    }
    return { re, im };
}

// Extract the lower half-spectrum (planes l = 0..N2/2) of a full grid into the
// compact (N0, N1, P) layout, P = N2/2 + 1, that fft3dHermitian expects.
function halfSpectrum(re: number[], im: number[], N0: number, N1: number, N2: number, Ctor: Float32ArrayConstructor | Float64ArrayConstructor) {
    const P = (N2 >> 1) + 1;
    const hr = new Ctor(N0 * N1 * P);
    const hi = new Ctor(N0 * N1 * P);
    for (let h = 0; h < N0; h++) {
        for (let k = 0; k < N1; k++) {
            for (let l = 0; l < P; l++) {
                const s = h * N1 * N2 + k * N2 + l;
                const d = h * N1 * P + k * P + l;
                hr[d] = re[s]; hi[d] = im[s];
            }
        }
    }
    return { hr: hr as FftArray, hi: hi as FftArray };
}

describe('fft3dHermitian', () => {
    const grids: [number, number, number][] = [[2, 4, 8], [4, 8, 16], [8, 8, 8], [8, 16, 4]];

    for (const Ctor of [Float64Array, Float32Array] as const) {
        const isF32 = Ctor === Float32Array;
        for (const [N0, N1, N2] of grids) {
            it(`matches the real part of the naive 3D DFT, ${N0}x${N1}x${N2} (${Ctor.name})`, () => {
                const total = N0 * N1 * N2;
                const { re, im } = makeHermitianGrid(N0, N1, N2, 555 + total);
                const ref = naiveDft3d(re, im, N0, N1, N2);

                // A Hermitian input must transform to a purely real output.
                let maxImag = 0;
                for (let i = 0; i < total; i++) maxImag = Math.max(maxImag, Math.abs(ref.im[i]));
                expect(maxImag).toBeLessThan(1e-9);

                const { hr, hi } = halfSpectrum(re, im, N0, N1, N2, Ctor);
                const out = makeArray(Ctor, new Array<number>(total).fill(0));
                fft3dHermitian(hr, hi, N0, N1, N2, out);

                const tol = isF32 ? 1e-2 : 1e-7;
                let maxErr = 0;
                for (let i = 0; i < total; i++) maxErr = Math.max(maxErr, Math.abs(out[i] - ref.re[i]));
                expect(maxErr).toBeLessThan(tol);
            });
        }
    }

    it('agrees with the general fft3d real output', () => {
        const N0 = 8, N1 = 8, N2 = 8, total = N0 * N1 * N2;
        const { re, im } = makeHermitianGrid(N0, N1, N2, 4242);

        const gr = makeArray(Float64Array, re), gi = makeArray(Float64Array, im);
        fft3d(gr, gi, N0, N1, N2);

        const { hr, hi } = halfSpectrum(re, im, N0, N1, N2, Float64Array);
        const out = makeArray(Float64Array, new Array<number>(total).fill(0));
        fft3dHermitian(hr, hi, N0, N1, N2, out);

        let maxErr = 0;
        for (let i = 0; i < total; i++) maxErr = Math.max(maxErr, Math.abs(out[i] - gr[i]));
        expect(maxErr).toBeLessThan(1e-9);
    });
});

describe('nextPow2', () => {
    it('rounds up to powers of two', () => {
        expect(nextPow2(1)).toBe(1);
        expect(nextPow2(2)).toBe(2);
        expect(nextPow2(3)).toBe(4);
        expect(nextPow2(5)).toBe(8);
        expect(nextPow2(8)).toBe(8);
        expect(nextPow2(17)).toBe(32);
        expect(nextPow2(1024)).toBe(1024);
        expect(nextPow2(1025)).toBe(2048);
    });

    it('handles n <= 1', () => {
        expect(nextPow2(0)).toBe(1);
        expect(nextPow2(-5)).toBe(1);
    });
});
