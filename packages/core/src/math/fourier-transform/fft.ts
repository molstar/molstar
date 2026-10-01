/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Float32Array or Float64Array for FFT data storage. */
export type FftArray = Float32Array | Float64Array;

/**
 * Returns the smallest integer >= n that is a power of 2.
 * Used to determine FFT-friendly grid dimensions.
 */
export function nextPow2(n: number): number {
    if (n <= 1) return 1;
    let p = 1;
    while (p < n) p <<= 1;
    return p;
}

// ---------------------------------------------------------------------------
// Precomputed twiddle-factor cache, keyed by transform size n.
// Always Float64 for trig precision; JS arithmetic is float64 regardless of
// the working array type.
// Entry: [cos(-2πk/n), sin(-2πk/n)] for k = 0..n-1 (full length so radix-4 can
// look up W^k, W^2k, W^3k directly), plus the bit-reversal permutation table.
// ---------------------------------------------------------------------------
const _twiddleCache = new Map<number, [Float64Array, Float64Array, Uint32Array]>();

function getTwiddle(n: number): [Float64Array, Float64Array, Uint32Array] {
    let entry = _twiddleCache.get(n);
    if (!entry) {
        const tr = new Float64Array(n);
        const ti = new Float64Array(n);
        for (let k = 0; k < n; k++) {
            const angle = -2 * Math.PI * k / n;
            tr[k] = Math.cos(angle);
            ti[k] = Math.sin(angle);
        }
        // Precompute bit-reversal permutation
        const perm = new Uint32Array(n);
        let j = 0;
        for (let i = 1; i < n; i++) {
            let bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            perm[i] = j;
        }
        entry = [tr, ti, perm];
        _twiddleCache.set(n, entry);
    }
    return entry;
}

// ---------------------------------------------------------------------------
// Persistent scratch buffers used by fft3d transpose passes.
// One pair per precision; grown on demand, never freed.
// ---------------------------------------------------------------------------
let _scratchRe32 = new Float32Array(0);
let _scratchIm32 = new Float32Array(0);
let _scratchRe64 = new Float64Array(0);
let _scratchIm64 = new Float64Array(0);

function getScratch<T extends FftArray>(n: number, ref: T): [T, T] {
    if (ref instanceof Float32Array) {
        if (_scratchRe32.length < n) {
            _scratchRe32 = new Float32Array(n);
            _scratchIm32 = new Float32Array(n);
        }
        return [_scratchRe32 as unknown as T, _scratchIm32 as unknown as T];
    } else {
        if (_scratchRe64.length < n) {
            _scratchRe64 = new Float64Array(n);
            _scratchIm64 = new Float64Array(n);
        }
        return [_scratchRe64 as unknown as T, _scratchIm64 as unknown as T];
    }
}

/**
 * In-place radix-4 (radix-2²) Cooley-Tukey FFT (DIT), with a single radix-2
 * stage when log2(n) is odd. Reuses the radix-2 bit-reversal ordering.
 * Computes the forward DFT: X[k] = Σ_{n=0}^{N-1} x[n] * exp(-2πi*n*k/N)
 *
 * @param re - real parts (modified in-place), length must be a power of 2
 * @param im - imaginary parts (modified in-place), length must be a power of 2
 * @param offset - starting index within the arrays
 * @param n - transform length (must be a power of 2)
 * @param stride - element stride in the arrays
 */
export function fft1d<T extends FftArray>(re: T, im: T, offset: number, n: number, stride: number): void {
    // Bit-reversal permutation (precomputed table)
    const [twRe, twIm, perm] = getTwiddle(n);
    for (let i = 1; i < n; i++) {
        const j = perm[i];
        if (i < j) {
            const ii = offset + i * stride;
            const jj = offset + j * stride;
            let tmp = re[ii]; re[ii] = re[jj]; re[jj] = tmp;
            tmp = im[ii]; im[ii] = im[jj]; im[jj] = tmp;
        }
    }

    // Butterfly passes (forward transform: -2πi sign).
    // Radix-4 (radix-2²): each stage quadruples the sub-transform size, halving
    // the number of memory passes vs. radix-2. When log2(n) is odd, one radix-2
    // stage runs first so the remaining stages pair up cleanly into radix-4.
    let log2n = 0;
    for (let t = n; t > 1; t >>= 1) log2n++;

    let m = 1; // current sub-transform size
    if ((log2n & 1) === 1) {
        // Single radix-2 stage: combine size-1 sub-transforms into size-2.
        for (let start = 0; start < n; start += 2) {
            const i0 = offset + start * stride;
            const i1 = offset + (start + 1) * stride;
            const aRe = re[i0], aIm = im[i0];
            const bRe = re[i1], bIm = im[i1];
            re[i0] = aRe + bRe; im[i0] = aIm + bIm;
            re[i1] = aRe - bRe; im[i1] = aIm - bIm;
        }
        m = 2;
    }

    // Radix-4 stages. For a block of size L = 4m the four contiguous quarters
    // o0..o3 (size m each) are combined. Twiddles: o1←W^{2k}, o2←W^{k}, o3←W^{3k}
    // (exponents scaled by step = n/L). The four-point DFT then folds in the
    // multiply-by-(±i) structure that arises from fusing two radix-2 stages.
    for (; m < n; m <<= 2) {
        const L = m << 2;
        const step = (n / L) | 0;
        const step2 = step << 1;
        const step3 = step2 + step;
        for (let start = 0; start < n; start += L) {
            let e1 = 0, e2 = 0, e3 = 0; // twiddle exponents k·step, 2k·step, 3k·step
            for (let k = 0; k < m; k++) {
                const o0 = offset + (start + k) * stride;
                const o1 = o0 + m * stride;
                const o2 = o1 + m * stride;
                const o3 = o2 + m * stride;

                const aRe = re[o0], aIm = im[o0];
                const bRe = re[o1], bIm = im[o1];
                const cRe = re[o2], cIm = im[o2];
                const dRe = re[o3], dIm = im[o3];

                const w1Re = twRe[e1], w1Im = twIm[e1]; // W^{k}   for o2
                const w2Re = twRe[e2], w2Im = twIm[e2]; // W^{2k}  for o1
                const w3Re = twRe[e3], w3Im = twIm[e3]; // W^{3k}  for o3

                // Twiddle-multiply the three "odd" quarters.
                const pRe = bRe * w2Re - bIm * w2Im;
                const pIm = bRe * w2Im + bIm * w2Re;
                const uRe = cRe * w1Re - cIm * w1Im;
                const uIm = cRe * w1Im + cIm * w1Re;
                const sRe = dRe * w3Re - dIm * w3Im;
                const sIm = dRe * w3Im + dIm * w3Re;

                const g0Re = aRe + pRe, g0Im = aIm + pIm;
                const g1Re = aRe - pRe, g1Im = aIm - pIm;
                const h0Re = uRe + sRe, h0Im = uIm + sIm;
                const h1Re = uRe - sRe, h1Im = uIm - sIm;

                re[o0] = g0Re + h0Re; im[o0] = g0Im + h0Im;
                re[o2] = g0Re - h0Re; im[o2] = g0Im - h0Im;
                re[o1] = g1Re + h1Im; im[o1] = g1Im - h1Re;
                re[o3] = g1Re - h1Im; im[o3] = g1Im + h1Re;

                e1 += step; e2 += step2; e3 += step3;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Transpose helpers (out-of-place, into pre-allocated buffers)
// ---------------------------------------------------------------------------

// Tile size for cache-blocked transposes. A tile of re+im stays within L1 for
// both Float32 and Float64; tunable.
const TRANSPOSE_BLOCK = 16;

/**
 * Cache-blocked 2D transpose: src (R×C, row-major) → dst (C×R, row-major).
 * srcOff / dstOff select a sub-block within larger buffers.
 * Tiling keeps the strided destination writes within cache.
 */
function transpose2dBlocked<T extends FftArray>(
    srcRe: T, srcIm: T,
    dstRe: T, dstIm: T,
    R: number, C: number,
    srcOff: number, dstOff: number
): void {
    const BLK = TRANSPOSE_BLOCK;
    for (let r0 = 0; r0 < R; r0 += BLK) {
        const rMax = r0 + BLK < R ? r0 + BLK : R;
        for (let c0 = 0; c0 < C; c0 += BLK) {
            const cMax = c0 + BLK < C ? c0 + BLK : C;
            for (let r = r0; r < rMax; r++) {
                const srcRow = srcOff + r * C;
                for (let c = c0; c < cMax; c++) {
                    const s = srcRow + c;
                    const d = dstOff + c * R + r;
                    dstRe[d] = srcRe[s];
                    dstIm[d] = srcIm[s];
                }
            }
        }
    }
}

/** Shape (A, B, C) → (A, C, B): swap last two axes. */
function transposeAxes12<T extends FftArray>(
    srcRe: T, srcIm: T,
    dstRe: T, dstIm: T,
    A: number, B: number, C: number
): void {
    const BC = B * C;
    for (let a = 0; a < A; a++) {
        transpose2dBlocked(srcRe, srcIm, dstRe, dstIm, B, C, a * BC, a * BC);
    }
}

/** Shape (N0, N1, N2) → (N1, N2, N0): axis 0 moved to last. */
function transposeAxis0ToLast<T extends FftArray>(
    srcRe: T, srcIm: T,
    dstRe: T, dstIm: T,
    N0: number, N1: number, N2: number
): void {
    transpose2dBlocked(srcRe, srcIm, dstRe, dstIm, N0, N1 * N2, 0, 0);
}

/** Shape (N1, N2, N0) → (N0, N1, N2): last axis moved to first. */
function transposeAxis0FromLast<T extends FftArray>(
    srcRe: T, srcIm: T,
    dstRe: T, dstIm: T,
    N0: number, N1: number, N2: number
): void {
    transpose2dBlocked(srcRe, srcIm, dstRe, dstIm, N1 * N2, N0, 0, 0);
}

/**
 * In-place 3D forward FFT.
 * Computes the forward DFT along all three dimensions.
 * The data is stored in row-major order: data[h * N1 * N2 + k * N2 + l].
 *
 * The forward DFT directly gives the electron density when the input is
 * populated with structure factors: ρ[n1,n2,n3] = Σ_{h,k,l} F[h,k,l] * exp(-2πi*(h*n1/N0+k*n2/N1+l*n3/N2))
 *
 * @param re - real parts (in-place), length N0*N1*N2
 * @param im - imaginary parts (in-place), length N0*N1*N2
 * @param N0 - size along axis 0 (must be power of 2)
 * @param N1 - size along axis 1 (must be power of 2)
 * @param N2 - size along axis 2 (must be power of 2)
 */
export function fft3d<T extends FftArray>(re: T, im: T, N0: number, N1: number, N2: number): void {
    const [scrRe, scrIm] = getScratch(N0 * N1 * N2, re);

    // Axis 2 — already contiguous (stride = 1), no transpose needed
    for (let i0 = 0; i0 < N0; i0++) {
        for (let i1 = 0; i1 < N1; i1++) {
            fft1d(re, im, i0 * N1 * N2 + i1 * N2, N2, 1);
        }
    }

    // Axis 1 — transpose (N0,N1,N2)→(N0,N2,N1), FFT at stride=1, transpose back
    transposeAxes12(re, im, scrRe, scrIm, N0, N1, N2);
    for (let i0 = 0; i0 < N0; i0++) {
        for (let i2 = 0; i2 < N2; i2++) {
            fft1d(scrRe, scrIm, i0 * N2 * N1 + i2 * N1, N1, 1);
        }
    }
    transposeAxes12(scrRe, scrIm, re, im, N0, N2, N1);

    // Axis 0 — transpose (N0,N1,N2)→(N1,N2,N0), FFT at stride=1, transpose back
    transposeAxis0ToLast(re, im, scrRe, scrIm, N0, N1, N2);
    for (let i1 = 0; i1 < N1; i1++) {
        for (let i2 = 0; i2 < N2; i2++) {
            fft1d(scrRe, scrIm, i1 * N2 * N0 + i2 * N0, N0, 1);
        }
    }
    transposeAxis0FromLast(scrRe, scrIm, re, im, N0, N1, N2);
}

/**
 * 3D forward FFT specialized for a Hermitian-symmetric input, producing a purely
 * real output. Roughly 2x faster than {@link fft3d} for this case and uses about
 * half the memory: only the non-redundant lower half-spectrum is stored.
 *
 * Exploits two facts that hold for crystallographic structure-factor synthesis:
 *  - the input grid F is Hermitian: F[-h,-k,-l] = conj(F[h,k,l]); hence
 *  - the forward DFT ρ[n] = Σ F·exp(-2πi·…) is real-valued.
 *
 * The input is the compact lower half-spectrum: planes l = 0..N2/2 only, in
 * row-major (N0, N1, P) layout with P = N2/2 + 1, i.e. index h*N1*P + k*P + l.
 * The caller fills these planes directly; the upper planes are implied by
 * Hermitian symmetry and never materialized.
 *
 * After transforming axes 0 and 1, every axis-2 line is conjugate-symmetric, so
 * the axis-2 pass uses a "two-for-one" real transform: two Hermitian lines are
 * packed into one complex length-N2 FFT whose real/imaginary outputs are the two
 * real result lines.
 *
 * The real result is written to `out`. The input arrays `re`/`im` are used as
 * scratch and their contents are destroyed.
 *
 * @param re - real parts of the half-spectrum (destroyed), length N0*N1*(N2/2+1)
 * @param im - imaginary parts of the half-spectrum (destroyed), length N0*N1*(N2/2+1)
 * @param N0 - size along axis 0 (must be power of 2)
 * @param N1 - size along axis 1 (must be power of 2)
 * @param N2 - size along axis 2 (must be power of 2)
 * @param out - real output (in-place), length N0*N1*N2, row-major
 */
export function fft3dHermitian<T extends FftArray>(re: T, im: T, N0: number, N1: number, N2: number, out: T): void {
    const M = N2 >> 1; // Nyquist index
    const P = M + 1; // number of stored l-planes (0..M); re/im hold the (N0,N1,P) half-spectrum
    const lines = N0 * N1;

    // Transpose scratch sized for the compact (N0,N1,P) half-spectrum. The input
    // re/im already hold that half-spectrum, so no compaction copy is needed.
    const [sRe, sIm] = getScratch(lines * P, re);

    // 2D FFT over axes 0 and 1 on the half-spectrum held in re/im, ping-ponging
    // with the scratch exactly like fft3d does.

    // Axis 1: (N0,N1,P)→(N0,P,N1), FFT along N1 at stride 1, transpose back.
    transposeAxes12(re, im, sRe, sIm, N0, N1, P);
    for (let i0 = 0; i0 < N0; i0++) {
        for (let p = 0; p < P; p++) {
            fft1d(sRe, sIm, i0 * P * N1 + p * N1, N1, 1);
        }
    }
    transposeAxes12(sRe, sIm, re, im, N0, P, N1);

    // Axis 0: (N0,N1,P)→(N1,P,N0), FFT along N0 at stride 1, transpose back.
    transposeAxis0ToLast(re, im, sRe, sIm, N0, N1, P);
    for (let i1 = 0; i1 < N1; i1++) {
        for (let p = 0; p < P; p++) {
            fft1d(sRe, sIm, i1 * P * N0 + p * N0, N0, 1);
        }
    }
    transposeAxis0FromLast(sRe, sIm, re, im, N0, N1, P);

    // Axis 2 — two-for-one real transform. re/im now hold G[a0,a1,l] for
    // l = 0..M in layout [q*P + l] with q = a0*N1 + a1, conjugate-symmetric along
    // axis 2. Pairs of lines q, q+1 are packed into one complex length-N2 FFT:
    // y = DFT(gA + i·gB) = ρA + i·ρB with ρA, ρB real.
    const yr = (re instanceof Float32Array ? new Float32Array(N2) : new Float64Array(N2)) as T;
    const yi = (re instanceof Float32Array ? new Float32Array(N2) : new Float64Array(N2)) as T;

    let q = 0;
    for (; q + 1 < lines; q += 2) {
        const baseA = q * P;
        const baseB = (q + 1) * P;

        // l = 0 (self-conjugate)
        yr[0] = re[baseA] - im[baseB];
        yi[0] = im[baseA] + re[baseB];
        // l = 1..M-1, with the mirror frequency N2-l from conjugate symmetry
        for (let l = 1; l < M; l++) {
            const aR = re[baseA + l], aI = im[baseA + l];
            const bR = re[baseB + l], bI = im[baseB + l];
            yr[l] = aR - bI; yi[l] = aI + bR;
            const ln = N2 - l;
            yr[ln] = aR + bI; yi[ln] = bR - aI;
        }
        // l = M (Nyquist, self-conjugate)
        yr[M] = re[baseA + M] - im[baseB + M];
        yi[M] = im[baseA + M] + re[baseB + M];

        fft1d(yr, yi, 0, N2, 1);

        const outA = q * N2;
        const outB = (q + 1) * N2;
        for (let a2 = 0; a2 < N2; a2++) {
            out[outA + a2] = yr[a2];
            out[outB + a2] = yi[a2];
        }
    }

    // Leftover single line (only when N0*N1 is odd, i.e. N0 = N1 = 1).
    if (q < lines) {
        const baseA = q * P;
        yr[0] = re[baseA]; yi[0] = im[baseA];
        for (let l = 1; l < M; l++) {
            const aR = re[baseA + l], aI = im[baseA + l];
            yr[l] = aR; yi[l] = aI;
            const ln = N2 - l;
            yr[ln] = aR; yi[ln] = -aI;
        }
        yr[M] = re[baseA + M]; yi[M] = im[baseA + M];

        fft1d(yr, yi, 0, N2, 1);

        const outA = q * N2;
        for (let a2 = 0; a2 < N2; a2++) out[outA + a2] = yr[a2];
    }
}
