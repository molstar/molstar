/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 *
 * Self-contained real spherical-harmonics utilities: associated Legendre
 * polynomials, the real SH basis Y_l^m, least-squares fitting of a radial
 * function r(theta, phi) and reconstruction from the fitted coefficients.
 *
 * The math is implemented from scratch; only Mol*'s generic dense `svd`
 * (mol-math/linear-algebra/matrix) is reused for the least-squares solve.
 *
 * Coefficient/basis ordering is flat by (l, m): index(l, m) = l*l + l + m,
 * so an expansion up to degree L has (L + 1)^2 terms.
 */

import { Matrix } from '../linear-algebra/matrix/matrix';
import { svd } from '../linear-algebra/matrix/svd';

/** Number of (l, m) terms for an expansion up to (and including) degree L. */
export function shTermCount(L: number): number {
    return (L + 1) * (L + 1);
}

/** Flat index of the (l, m) term, m in [-l, l]. */
export function shIndex(l: number, m: number): number {
    return l * l + l + m;
}

/** Packed lower-triangular index for 0 <= m <= l. */
function legendreIndex(l: number, m: number): number {
    return (l * (l + 1)) / 2 + m;
}

/**
 * Associated Legendre polynomials P_l^m(x), m >= 0, computed with the standard
 * stable recurrences. Results are written into `out` indexed by `legendreIndex`.
 * `x` is expected in [-1, 1] (typically cos(theta)).
 */
export function assocLegendre(L: number, x: number, out?: Float64Array): Float64Array {
    const size = ((L + 1) * (L + 2)) / 2;
    const p = out && out.length >= size ? out : new Float64Array(size);

    // P_0^0 = 1
    p[legendreIndex(0, 0)] = 1;
    if (L === 0) return p;

    const somx2 = Math.sqrt(Math.max(0, 1 - x * x)); // sin(theta)

    // Diagonal terms P_m^m = (-1)^m (2m-1)!! (1-x^2)^(m/2)
    let pmm = 1;
    let fact = 1;
    for (let m = 1; m <= L; ++m) {
        pmm *= -fact * somx2;
        fact += 2;
        p[legendreIndex(m, m)] = pmm;
    }

    // P_{m+1}^m = x (2m+1) P_m^m
    for (let m = 0; m < L; ++m) {
        p[legendreIndex(m + 1, m)] = x * (2 * m + 1) * p[legendreIndex(m, m)];
    }

    // Upward recurrence in l:
    // (l - m) P_l^m = x (2l - 1) P_{l-1}^m - (l + m - 1) P_{l-2}^m
    for (let m = 0; m <= L; ++m) {
        for (let l = m + 2; l <= L; ++l) {
            p[legendreIndex(l, m)] = (x * (2 * l - 1) * p[legendreIndex(l - 1, m)] -
                (l + m - 1) * p[legendreIndex(l - 2, m)]) / (l - m);
        }
    }

    return p;
}

// Precomputed normalization factors K_l^m = sqrt((2l+1)/(4pi) * (l-|m|)!/(l+|m|)!)
const normCache = new Map<number, Float64Array>();
function getNormFactors(L: number): Float64Array {
    const cached = normCache.get(L);
    if (cached) return cached;

    const size = ((L + 1) * (L + 2)) / 2;
    const norm = new Float64Array(size);
    const fourPi = 4 * Math.PI;
    for (let l = 0; l <= L; ++l) {
        for (let m = 0; m <= l; ++m) {
            // (l-m)!/(l+m)! computed as a running product to avoid large factorials
            let ratio = 1;
            for (let k = l - m + 1; k <= l + m; ++k) ratio /= k;
            norm[legendreIndex(l, m)] = Math.sqrt(((2 * l + 1) / fourPi) * ratio);
        }
    }
    normCache.set(L, norm);
    return norm;
}

/**
 * Evaluate the real SH basis Y_l^m(theta, phi) for all 0 <= l <= L, -l <= m <= l.
 * theta is the polar angle [0, pi], phi the azimuth [-pi, pi].
 * Results are written into `out` (length shTermCount(L)).
 *
 * Real form:
 *   m > 0:  sqrt(2) K_l^m cos(m phi) P_l^m(cos theta)
 *   m = 0:  K_l^0 P_l^0(cos theta)
 *   m < 0:  sqrt(2) K_l^|m| sin(|m| phi) P_l^|m|(cos theta)
 */
export function realSph(L: number, theta: number, phi: number, out?: Float64Array, legendreScratch?: Float64Array): Float64Array {
    const y = out && out.length >= shTermCount(L) ? out : new Float64Array(shTermCount(L));
    const p = assocLegendre(L, Math.cos(theta), legendreScratch);
    const norm = getNormFactors(L);
    const sqrt2 = Math.SQRT2;

    for (let l = 0; l <= L; ++l) {
        // m = 0
        y[shIndex(l, 0)] = norm[legendreIndex(l, 0)] * p[legendreIndex(l, 0)];
        for (let m = 1; m <= l; ++m) {
            const k = norm[legendreIndex(l, m)] * p[legendreIndex(l, m)] * sqrt2;
            y[shIndex(l, m)] = k * Math.cos(m * phi);
            y[shIndex(l, -m)] = k * Math.sin(m * phi);
        }
    }
    return y;
}

/** Spherical coordinates of a point relative to a center. */
export interface SphericalCoord { r: number, theta: number, phi: number }
export function toSpherical(x: number, y: number, z: number, out?: SphericalCoord): SphericalCoord {
    const o = out ?? { r: 0, theta: 0, phi: 0 };
    const r = Math.sqrt(x * x + y * y + z * z);
    o.r = r;
    o.theta = r > 1e-12 ? Math.acos(Math.min(1, Math.max(-1, z / r))) : 0;
    o.phi = Math.atan2(y, x);
    return o;
}

export interface SphericalHarmonicFit {
    /** Fitted coefficients, length shTermCount(L). */
    coeffs: Float64Array
    L: number
    /** Largest observed sample radius about the center; used to clamp reconstruction overshoot. */
    rMax: number
}

/**
 * Least-squares fit of the radial function r(theta, phi) of `points` (flat xyz
 * array) about `center` to a real spherical-harmonic expansion up to degree L.
 *
 * Builds the symmetric normal-equations system A = B^T B (size K x K,
 * K = (L+1)^2) and rhs = B^T r without materializing the full design matrix,
 * then solves via a truncated SVD pseudo-inverse for stability.
 *
 * `maxPoints` caps the number of samples used for fitting (stride subsampling)
 * to keep high-L fits affordable; reconstruction quality is unaffected.
 *
 * `regularization` adds a scale-invariant Tikhonov term to the normal-equations
 * diagonal with per-band `l(l+1)` damping (a smoothness prior). It is the main
 * defense against the ill-conditioned, exploding fits that arise from sparse or
 * clustered samples (e.g. trace-only or low-resolution clouds): high-degree bands
 * are penalized so an under-determined system relaxes toward the mean-radius
 * sphere instead of diverging. 0 disables it.
 */
export function fitSphericalHarmonics(points: ArrayLike<number>, center: ArrayLike<number>, L: number, maxPoints = 8192, regularization = 0): SphericalHarmonicFit {
    const K = shTermCount(L);
    const A = Matrix.create(K, K, Float64Array); // K x K, symmetric
    const rhs = new Float64Array(K);

    const cx = center[0], cy = center[1], cz = center[2];
    const pointCount = Math.floor(points.length / 3);
    const stride = Math.max(1, Math.ceil(pointCount / maxPoints));

    const basis = new Float64Array(K);
    const legendreScratch = new Float64Array(((L + 1) * (L + 2)) / 2);
    const sc: SphericalCoord = { r: 0, theta: 0, phi: 0 };
    let rMax = 0;

    for (let i = 0; i < pointCount; i += stride) {
        toSpherical(points[i * 3] - cx, points[i * 3 + 1] - cy, points[i * 3 + 2] - cz, sc);
        if (sc.r <= 1e-12) continue;
        if (sc.r > rMax) rMax = sc.r;
        realSph(L, sc.theta, sc.phi, basis, legendreScratch);

        // accumulate A += basis basis^T (upper triangle) and rhs += basis * r
        for (let a = 0; a < K; ++a) {
            const ba = basis[a];
            if (ba === 0) continue;
            rhs[a] += ba * sc.r;
            for (let b = a; b < K; ++b) {
                A.data[a * K + b] += ba * basis[b];
            }
        }
    }
    // mirror to full symmetric matrix
    for (let a = 0; a < K; ++a) {
        for (let b = a + 1; b < K; ++b) {
            A.data[b * K + a] = A.data[a * K + b];
        }
    }

    // scale-invariant Tikhonov regularization with per-band l(l+1) damping:
    // A[a][a] += regularization * meanDiag * (1 + l(l+1)), l = band of flat index a
    if (regularization > 0) {
        let trace = 0;
        for (let a = 0; a < K; ++a) trace += A.data[a * K + a];
        const meanDiag = trace / K;
        if (meanDiag > 0) {
            for (let a = 0; a < K; ++a) {
                const l = Math.floor(Math.sqrt(a)); // indices [l*l, (l+1)*(l+1)-1] belong to band l
                A.data[a * K + a] += regularization * meanDiag * (1 + l * (l + 1));
            }
        }
    }

    // With Tikhonov regularization the normal matrix is SPD, so solve it exactly via Cholesky -
    // several times faster than the SVD pseudo-inverse, and since the fit runs once per unit the
    // solve dominates build/update time on large multi-chain structures. Fall back to the truncated
    // SVD pseudo-inverse when unregularized (and possibly rank-deficient) or if Cholesky fails.
    if (regularization > 0) {
        const chol = choleskySolveSymmetric(A.data, K, rhs);
        if (chol) return { coeffs: chol, L, rMax };
    }

    // Solve A x = rhs via truncated SVD pseudo-inverse: x = V diag(1/w) U^T rhs
    const W = Matrix.create(1, K, Float64Array);
    const U = Matrix.create(K, K, Float64Array);
    const V = Matrix.create(K, K, Float64Array);
    svd(A, W, U, V);

    // singular-value truncation relative to the largest
    let wMax = 0;
    for (let j = 0; j < K; ++j) wMax = Math.max(wMax, W.data[j]);
    const tol = wMax * 1e-8;

    // y = U^T rhs (column j of U dotted with rhs), z = y / w
    const z = new Float64Array(K);
    for (let j = 0; j < K; ++j) {
        const w = W.data[j];
        if (w <= tol) { z[j] = 0; continue; }
        let y = 0;
        for (let r = 0; r < K; ++r) y += U.data[r * K + j] * rhs[r];
        z[j] = y / w;
    }
    // x = V z
    const coeffs = new Float64Array(K);
    for (let i = 0; i < K; ++i) {
        let s = 0;
        for (let j = 0; j < K; ++j) s += V.data[i * K + j] * z[j];
        coeffs[i] = s;
    }

    return { coeffs, L, rMax };
}

/**
 * Solve a symmetric positive-definite system `A x = rhs` (A row-major, K x K) by Cholesky
 * factorization. Returns `undefined` if A is not positive definite (a non-positive pivot from
 * round-off, or a rank-deficient/unregularized system), letting the caller fall back to SVD.
 *
 * For the regularized normal equations this is exact and several times faster than the SVD
 * pseudo-inverse - and the fit runs once per unit, so the solve dominates representation
 * build/update time on large multi-chain structures.
 */
function choleskySolveSymmetric(A: ArrayLike<number>, K: number, rhs: ArrayLike<number>): Float64Array | undefined {
    const Lm = new Float64Array(K * K); // lower-triangular Cholesky factor
    for (let i = 0; i < K; ++i) {
        for (let j = 0; j <= i; ++j) {
            let sum = A[i * K + j];
            for (let k = 0; k < j; ++k) sum -= Lm[i * K + k] * Lm[j * K + k];
            if (i === j) {
                if (!(sum > 0)) return undefined; // not positive definite (also catches NaN)
                Lm[i * K + i] = Math.sqrt(sum);
            } else {
                Lm[i * K + j] = sum / Lm[j * K + j];
            }
        }
    }
    const x = new Float64Array(K);
    // forward substitution L y = rhs (y accumulated in x)
    for (let i = 0; i < K; ++i) {
        let sum = rhs[i];
        for (let k = 0; k < i; ++k) sum -= Lm[i * K + k] * x[k];
        x[i] = sum / Lm[i * K + i];
    }
    // back substitution L^T x = y
    for (let i = K - 1; i >= 0; --i) {
        let sum = x[i];
        for (let k = i + 1; k < K; ++k) sum -= Lm[k * K + i] * x[k];
        x[i] = sum / Lm[i * K + i];
    }
    return x;
}

/** Reconstruct the radius from fitted coefficients at the given direction. */
export function reconstructRadius(coeffs: ArrayLike<number>, L: number, theta: number, phi: number, basisScratch?: Float64Array, legendreScratch?: Float64Array): number {
    const K = shTermCount(L);
    const basis = realSph(L, theta, phi, basisScratch, legendreScratch);
    let r = 0;
    for (let i = 0; i < K; ++i) r += coeffs[i] * basis[i];
    return r;
}

/** One star-shaped lobe: a radial SH expansion about its own center. */
export interface SphericalHarmonicLobe { center: number[], coeffs: Float64Array, rMax: number }
export interface SphericalHarmonicLobesFit { lobes: SphericalHarmonicLobe[], L: number }

/** Gather a subset of an interleaved xyz cloud into a flat array and its centroid. */
function gatherSubset(points: ArrayLike<number>, idx: ArrayLike<number>): { flat: Float64Array, center: number[] } {
    const m = idx.length;
    const flat = new Float64Array(m * 3);
    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < m; ++i) {
        const p = idx[i] * 3;
        const x = points[p], y = points[p + 1], z = points[p + 2];
        flat[i * 3] = x; flat[i * 3 + 1] = y; flat[i * 3 + 2] = z;
        cx += x; cy += y; cz += z;
    }
    return { flat, center: m > 0 ? [cx / m, cy / m, cz / m] : [0, 0, 0] };
}

/** Small deterministic PRNG (mulberry32) for seeded, reproducible k-means++ initialization. */
function mulberry32(seed: number): () => number {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

/**
 * Partition an interleaved xyz cloud into `k` spatial clusters by k-means (k-means++ seeded init,
 * Lloyd iterations) and return a per-point cluster label. Deterministic for a given `seed`, so the
 * same cloud yields the same partition (the fit caches on it). Spatial clusters give locally compact,
 * star-convex lobes - for a threaded chain (e.g. rRNA) far more compact per cluster than contiguous
 * sequence runs.
 *
 * The Lloyd iterations (the O(points * k * iterations) cost) run on a strided subsample of at most
 * `maxSamples` points to find the cluster centers; every point is then labelled by its nearest center
 * in a single final pass. Subsampling a large cloud changes the centers negligibly but cuts the
 * dominant cost several-fold (e.g. a 63k-atom rRNA chain).
 */
export function kmeansLabels(points: ArrayLike<number>, k: number, options: { maxIterations?: number, seed?: number, maxSamples?: number } = {}): Int32Array {
    const n = Math.floor(points.length / 3);
    const labels = new Int32Array(n);
    const K = Math.min(Math.max(1, Math.round(k)), n);
    if (K <= 1 || n === 0) return labels;

    const maxIterations = options.maxIterations ?? 20;
    const maxSamples = options.maxSamples ?? 8192;
    const rand = mulberry32((options.seed ?? 0x9e3779b9) >>> 0);

    // strided subsample that the iterative center search runs on (all points are labelled at the end)
    const stride = Math.max(1, Math.floor(n / maxSamples));
    const sn = Math.floor((n + stride - 1) / stride);
    const sample = new Float64Array(sn * 3);
    for (let s = 0; s < sn; ++s) {
        const i = s * stride;
        sample[s * 3] = points[i * 3]; sample[s * 3 + 1] = points[i * 3 + 1]; sample[s * 3 + 2] = points[i * 3 + 2];
    }

    // k-means++ seeding: spread initial centers by squared-distance probability
    const centers = new Float64Array(K * 3);
    const d2 = new Float64Array(sn).fill(Infinity);
    let pick = Math.floor(rand() * sn);
    for (let c = 0; c < K; ++c) {
        centers[c * 3] = sample[pick * 3]; centers[c * 3 + 1] = sample[pick * 3 + 1]; centers[c * 3 + 2] = sample[pick * 3 + 2];
        if (c === K - 1) break;
        let sum = 0;
        for (let i = 0; i < sn; ++i) {
            const dx = sample[i * 3] - centers[c * 3], dy = sample[i * 3 + 1] - centers[c * 3 + 1], dz = sample[i * 3 + 2] - centers[c * 3 + 2];
            const dd = dx * dx + dy * dy + dz * dz;
            if (dd < d2[i]) d2[i] = dd;
            sum += d2[i];
        }
        let target = rand() * sum;
        pick = sn - 1;
        for (let i = 0; i < sn; ++i) { target -= d2[i]; if (target <= 0) { pick = i; break; } }
    }

    // Lloyd iterations on the subsample
    const sLabels = new Int32Array(sn);
    const sums = new Float64Array(K * 3);
    const counts = new Int32Array(K);
    for (let it = 0; it < maxIterations; ++it) {
        let changed = false;
        for (let i = 0; i < sn; ++i) {
            let best = Infinity, bi = 0;
            for (let c = 0; c < K; ++c) {
                const dx = sample[i * 3] - centers[c * 3], dy = sample[i * 3 + 1] - centers[c * 3 + 1], dz = sample[i * 3 + 2] - centers[c * 3 + 2];
                const dd = dx * dx + dy * dy + dz * dz;
                if (dd < best) { best = dd; bi = c; }
            }
            if (sLabels[i] !== bi) { sLabels[i] = bi; changed = true; }
        }
        if (!changed && it > 0) break;
        sums.fill(0); counts.fill(0);
        for (let i = 0; i < sn; ++i) {
            const c = sLabels[i];
            sums[c * 3] += sample[i * 3]; sums[c * 3 + 1] += sample[i * 3 + 1]; sums[c * 3 + 2] += sample[i * 3 + 2];
            counts[c]++;
        }
        for (let c = 0; c < K; ++c) {
            if (counts[c] > 0) { centers[c * 3] = sums[c * 3] / counts[c]; centers[c * 3 + 1] = sums[c * 3 + 1] / counts[c]; centers[c * 3 + 2] = sums[c * 3 + 2] / counts[c]; }
        }
    }

    // label every point by its nearest final center
    for (let i = 0; i < n; ++i) {
        let best = Infinity, bi = 0;
        for (let c = 0; c < K; ++c) {
            const dx = points[i * 3] - centers[c * 3], dy = points[i * 3 + 1] - centers[c * 3 + 1], dz = points[i * 3 + 2] - centers[c * 3 + 2];
            const dd = dx * dx + dy * dy + dz * dz;
            if (dd < best) { best = dd; bi = c; }
        }
        labels[i] = bi;
    }
    return labels;
}

/**
 * Fit one star-shaped radial SH lobe per distinct label: points sharing the same `labels[i]`
 * form one lobe. The caller decides the partition (e.g. contiguous sequence runs, k-means clusters).
 *
 * Each lobe is centered on the centroid of its (original) points. When `radii` is given, every point
 * is inflated outward from that centroid by its radius *after* the centroid is computed, so the fitted
 * envelope tracks the outer surface (atom size + probe) while the center stays put. Inflating the
 * points in 3D *before* the centroid would shift the center proportionally to the offset on an
 * asymmetric cluster, drifting the lobe off the atoms - so the inflation must follow the centroid.
 */
export function fitSphericalHarmonicLobesByLabel(points: ArrayLike<number>, labels: ArrayLike<number>, L: number, options: { maxPoints?: number, regularization?: number, radii?: ArrayLike<number> } = {}): SphericalHarmonicLobesFit {
    const maxPoints = options.maxPoints ?? 8192;
    const regularization = options.regularization ?? 0;
    const radii = options.radii;
    const n = Math.floor(points.length / 3);

    const byLabel = new Map<number, number[]>();
    for (let i = 0; i < n; ++i) {
        const lbl = labels[i];
        let arr = byLabel.get(lbl);
        if (!arr) { arr = []; byLabel.set(lbl, arr); }
        arr.push(i);
    }

    const lobes: SphericalHarmonicLobe[] = [];
    for (const idxArr of byLabel.values()) {
        const idx = Int32Array.from(idxArr);
        const { flat, center } = gatherSubset(points, idx);
        if (radii) {
            const cx = center[0], cy = center[1], cz = center[2];
            for (let t = 0; t < idx.length; ++t) {
                const p = t * 3;
                const dx = flat[p] - cx, dy = flat[p + 1] - cy, dz = flat[p + 2] - cz;
                const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
                if (len > 1e-3) {
                    const s = radii[idx[t]] / len; // displace by `radius` along the centroid->point ray
                    flat[p] += dx * s; flat[p + 1] += dy * s; flat[p + 2] += dz * s;
                }
            }
        }
        const { coeffs, rMax } = fitSphericalHarmonics(flat, center, L, maxPoints, regularization);
        lobes.push({ center, coeffs, rMax });
    }
    return { lobes, L };
}
