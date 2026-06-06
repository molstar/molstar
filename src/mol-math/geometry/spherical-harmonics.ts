/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

/**
 * Associated Legendre polynomials P_l^m(x), m >= 0, computed with the standard
 * stable recurrences. Results are written into `out` indexed by `legendreIndex`.
 * `x` is expected in [-1, 1] (typically cos(theta)).
 */
function legendreIndex(l: number, m: number): number {
    // packed lower-triangular storage for 0 <= m <= l <= L
    return (l * (l + 1)) / 2 + m;
}

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
export function realSphAdvanced(L: number, theta: number, phi: number, out?: Float64Array, legendreScratch?: Float64Array): Float64Array {
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
 */
export function fitSphericalHarmonics(points: ArrayLike<number>, center: ArrayLike<number>, L: number, maxPoints = 8192): SphericalHarmonicFit {
    const K = shTermCount(L);
    const A = Matrix.create(K, K, Float64Array); // K x K, symmetric
    const rhs = new Float64Array(K);

    const cx = center[0], cy = center[1], cz = center[2];
    const pointCount = Math.floor(points.length / 3);
    const stride = Math.max(1, Math.ceil(pointCount / maxPoints));

    const basis = new Float64Array(K);
    const legendreScratch = new Float64Array(((L + 1) * (L + 2)) / 2);
    const sc: SphericalCoord = { r: 0, theta: 0, phi: 0 };

    for (let i = 0; i < pointCount; i += stride) {
        toSpherical(points[i * 3] - cx, points[i * 3 + 1] - cy, points[i * 3 + 2] - cz, sc);
        if (sc.r <= 1e-12) continue;
        realSphAdvanced(L, sc.theta, sc.phi, basis, legendreScratch);

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

    return { coeffs, L };
}

/** Reconstruct the radius from fitted coefficients at the given direction. */
export function reconstructRadius(coeffs: ArrayLike<number>, L: number, theta: number, phi: number, basisScratch?: Float64Array, legendreScratch?: Float64Array): number {
    const K = shTermCount(L);
    const basis = realSphAdvanced(L, theta, phi, basisScratch, legendreScratch);
    let r = 0;
    for (let i = 0; i < K; ++i) r += coeffs[i] * basis[i];
    return r;
}
