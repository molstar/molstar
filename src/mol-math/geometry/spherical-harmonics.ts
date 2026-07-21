/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
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

import { fastAcos, fastAtan2, fastCos, fastSin } from '../approx';
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
// Indexed directly by L (a small non-negative integer) rather than a `Map`: a plain array lookup
// avoids the hashing/bucket overhead of `Map.get`, which matters here since this is looked up on
// every `realSph` call unless the caller supplies a precomputed `norm` (see below) - `realSph` is
// called once per sample point during fitting and once per grid voxel during reconstruction, i.e.
// this is one of the hottest call sites in the whole module.
const normCache: (Float64Array | undefined)[] = [];
export function getNormFactors(L: number): Float64Array {
    const cached = normCache[L];
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
    normCache[L] = norm;
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
export function realSph(L: number, theta: number, phi: number, out: Float64Array, legendreScratch: Float64Array, norm: Float64Array): Float64Array {
    const y = out && out.length >= shTermCount(L) ? out : new Float64Array(shTermCount(L));
    const p = assocLegendre(L, fastCos(theta), legendreScratch);
    const sqrt2 = Math.SQRT2;

    // m = 0 band: no azimuthal (phi) dependence
    for (let l = 0; l <= L; ++l) {
        y[shIndex(l, 0)] = norm[legendreIndex(l, 0)] * p[legendreIndex(l, 0)];
    }

    // cos(m phi)/sin(m phi) via the angle-addition recurrence from a single cos(phi)/sin(phi)
    // evaluation, instead of calling Math.cos/Math.sin(m * phi) directly for every (l, m) pair -
    // the original loop nested m inside l, so it recomputed the SAME cos(m phi)/sin(m phi)
    // (which does not depend on l at all) once per l >= m, i.e. up to O(L^2) redundant
    // transcendental calls. This brings it down to 2 `Math.cos`/`Math.sin` calls total plus O(L)
    // cheap multiply-adds.
    if (L > 0) {
        const cos1 = fastCos(phi), sin1 = fastSin(phi);
        let cm = 1, sm = 0; // cos(0 * phi), sin(0 * phi)
        for (let m = 1; m <= L; ++m) {
            const cmNext = cos1 * cm - sin1 * sm;
            const smNext = sin1 * cm + cos1 * sm;
            cm = cmNext; sm = smNext;

            for (let l = m; l <= L; ++l) {
                const k = norm[legendreIndex(l, m)] * p[legendreIndex(l, m)] * sqrt2;
                y[shIndex(l, m)] = k * cm;
                y[shIndex(l, -m)] = k * sm;
            }
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
    o.theta = r > 1e-12 ? fastAcos(Math.min(1, Math.max(-1, z / r))) : 0;
    o.phi = fastAtan2(y, x);
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
    // fetched once and reused for every sample - `realSph` would otherwise redo this lookup
    // (previously a Map.get, now an array index) on every single one of the up to `maxPoints`
    // calls below, for a value that never changes within one fit
    const norm = getNormFactors(L);

    for (let i = 0; i < pointCount; i += stride) {
        toSpherical(points[i * 3] - cx, points[i * 3 + 1] - cy, points[i * 3 + 2] - cz, sc);
        if (sc.r <= 1e-12) continue;
        if (sc.r > rMax) rMax = sc.r;
        realSph(L, sc.theta, sc.phi, basis, legendreScratch, norm);

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

/** Reconstruct the radius from fitted coefficients at the given direction. `norm` can be a
 * precomputed `getNormFactors(L)` result, saving a lookup when called repeatedly for the same
 * `L` (e.g. once per grid voxel while reconstructing a blob's surface). */
export function reconstructRadius(coeffs: ArrayLike<number>, L: number, theta: number, phi: number, basisScratch: Float64Array, legendreScratch: Float64Array, norm: Float64Array): number {
    const K = shTermCount(L);
    const basis = realSph(L, theta, phi, basisScratch, legendreScratch, norm);
    let r = 0;
    for (let i = 0; i < K; ++i) r += coeffs[i] * basis[i];
    return r;
}

/** One star-shaped lobe: a radial SH expansion about its own center. */
export interface SphericalHarmonicLobe { center: number[], coeffs: Float64Array, rMax: number }
export interface SphericalHarmonicLobesFit { lobes: SphericalHarmonicLobe[], L: number }

/**
 * A bilinear-interpolatable lookup table of `R(theta, phi)` over a uniform `(nTheta+1) x nPhi`
 * grid spanning `theta` in `[0, pi]` (both poles included) and `phi` in `[-pi, pi)` (wrapping).
 */
export interface RadiusLUT { nTheta: number, nPhi: number, values: Float64Array }

/**
 * Precomputes a `RadiusLUT` for the fitted radial boundary `R(theta, phi) = sum coeffs[k] *
 * Y_k(theta, phi)`.
 *
 * ALGORITHMIC MOTIVATION: evaluating `R` directly (`reconstructRadius`) costs O(L^2) per call (a
 * fresh associated-Legendre table plus O(L) trig recurrence, since both depend on the specific
 * direction). Caching a scalar or a few intermediate values does not help when `R` must be
 * evaluated at a very large number of DIFFERENT directions - which is exactly what
 * `computeBlobSurface`'s marching-cubes splat loop does, calling it once per grid voxel within a
 * blob's padded box (easily thousands to millions of directions for one blob). Building this
 * table ONCE per blob - a fixed O(tableSize * L^2) cost, independent of how many voxels will
 * query it - and then looking values up via `sampleRadiusLUT` (O(1): 4 array reads + a bilinear
 * blend, no trigonometry or Legendre recursion at all) changes the TOTAL complexity of
 * reconstructing one blob's surface from O(voxels * L^2) to O(tableSize * L^2 + voxels), which is
 * a genuine asymptotic win whenever (as is typical) voxels >> tableSize.
 *
 * `minR`/`maxR` clamp each table entry once at build time (rather than after every lookup).
 * `out`, if given and the right size, has its `values` array reused instead of reallocating -
 * useful when rebuilding the table once per blob but the table dimensions (which only depend on
 * the shared `shDegree` prop, not on any individual blob) stay constant across all of them.
 */
export function buildRadiusLUT(coeffs: ArrayLike<number>, L: number, nTheta: number, nPhi: number, minR = -Infinity, maxR = Infinity, out?: RadiusLUT): RadiusLUT {
    const size = (nTheta + 1) * nPhi;
    const values = out && out.nTheta === nTheta && out.nPhi === nPhi && out.values.length === size ? out.values : new Float64Array(size);

    const K = shTermCount(L);
    const basis = new Float64Array(K);
    const legendreScratch = new Float64Array(((L + 1) * (L + 2)) / 2);
    const norm = getNormFactors(L);

    for (let it = 0; it <= nTheta; ++it) {
        const theta = (it / nTheta) * Math.PI;
        for (let ip = 0; ip < nPhi; ++ip) {
            const phi = (ip / nPhi) * 2 * Math.PI - Math.PI;
            realSph(L, theta, phi, basis, legendreScratch, norm);

            let r = 0;
            for (let k = 0; k < K; ++k) r += coeffs[k] * basis[k];
            if (r < minR) r = minR; else if (r > maxR) r = maxR;
            values[it * nPhi + ip] = r;
        }
    }
    return { nTheta, nPhi, values };
}

/**
 * Bilinearly samples a `RadiusLUT` (built by `buildRadiusLUT`) at the given direction - O(1): 4
 * array reads plus a bilinear blend, no trigonometry or Legendre recursion.
 */
export function sampleRadiusLUT(lut: RadiusLUT, theta: number, phi: number): number {
    const { nTheta, nPhi, values } = lut;

    let tt = (theta / Math.PI) * nTheta;
    if (tt < 0) tt = 0; else if (tt > nTheta) tt = nTheta;
    let it0 = Math.floor(tt);
    if (it0 >= nTheta) it0 = nTheta - 1;
    const it1 = it0 + 1;
    const ft = tt - it0;

    // phi wraps around the full circle: no separate last column is stored for phi = pi, it is
    // the same direction as phi = -pi (column 0), so ip1 wraps back to 0 past the last column
    let pp = ((phi + Math.PI) / (2 * Math.PI)) * nPhi;
    if (pp < 0) pp = 0; else if (pp >= nPhi) pp -= nPhi;
    let ip0 = Math.floor(pp);
    if (ip0 >= nPhi) ip0 = nPhi - 1;
    const fp = pp - ip0;
    let ip1 = ip0 + 1;
    if (ip1 >= nPhi) ip1 -= nPhi;

    const row0 = it0 * nPhi, row1 = it1 * nPhi;
    const v00 = values[row0 + ip0], v01 = values[row0 + ip1];
    const v10 = values[row1 + ip0], v11 = values[row1 + ip1];

    const v0 = v00 + (v01 - v00) * fp;
    const v1 = v10 + (v11 - v10) * fp;
    return v0 + (v1 - v0) * ft;
}
