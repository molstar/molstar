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
import { PrincipalAxes } from '../linear-algebra/matrix/principal-axes';

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

/** Reconstruct the radius from fitted coefficients at the given direction. */
export function reconstructRadius(coeffs: ArrayLike<number>, L: number, theta: number, phi: number, basisScratch?: Float64Array, legendreScratch?: Float64Array): number {
    const K = shTermCount(L);
    const basis = realSph(L, theta, phi, basisScratch, legendreScratch);
    let r = 0;
    for (let i = 0; i < K; ++i) r += coeffs[i] * basis[i];
    return r;
}

/**
 * Non-star-shapedness of a point cloud about `center`: the fraction of populated
 * direction bins that contain a radial GAP wider than `thickness` — an empty
 * shell of radius along the ray, i.e. the surface is crossed, left, and re-entered.
 *
 * A radial gap (not the min-max spread) is the defect that actually breaks a
 * single-center radial expansion. It works on both a thin surface shell and a
 * solid atom cloud: a star-shaped blob is radially filled from the center out
 * (no gap), whereas a concavity, C-shape, or a lobe seen past an empty region
 * leaves a gap. Using the spread instead would flag every solid cloud (interior
 * points reach r≈0) and every bumpy-but-star-shaped blob, over-splitting both.
 * It is also independent of the expansion degree L.
 *
 * Bins use an equal-area layout: latitude by cos(theta) in [-1, 1] and longitude
 * by phi in [-pi, pi].
 */
export function starShapeViolation(points: ArrayLike<number>, center: ArrayLike<number>, thickness: number, nLat = 12, nLon = 24): number {
    const cx = center[0], cy = center[1], cz = center[2];
    const nBins = nLat * nLon;
    const binRadii: number[][] = [];
    for (let b = 0; b < nBins; ++b) binRadii.push([]);

    const n = Math.floor(points.length / 3);
    const twoPi = 2 * Math.PI;
    for (let i = 0; i < n; ++i) {
        const dx = points[i * 3] - cx, dy = points[i * 3 + 1] - cy, dz = points[i * 3 + 2] - cz;
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r <= 1e-12) continue;
        const ct = Math.min(1, Math.max(-1, dz / r)); // cos(theta)
        const phi = Math.atan2(dy, dx); // [-pi, pi]
        let li = Math.floor(((ct + 1) / 2) * nLat); if (li >= nLat) li = nLat - 1; else if (li < 0) li = 0;
        let oi = Math.floor(((phi + Math.PI) / twoPi) * nLon); if (oi >= nLon) oi = nLon - 1; else if (oi < 0) oi = 0;
        binRadii[li * nLon + oi].push(r);
    }

    // A bin needs enough samples for a radial gap to mean "empty region" rather than
    // "sparse sampling"; sparse/solid clouds (e.g. atom positions, ~handful per bin)
    // are simply not judged, so they are never spuriously split.
    const minBinCount = 6;
    let populated = 0, violating = 0;
    for (let b = 0; b < nBins; ++b) {
        const radii = binRadii[b];
        if (radii.length < minBinCount) continue;
        ++populated;
        radii.sort((p, q) => p - q);
        let maxGap = 0;
        for (let i = 1; i < radii.length; ++i) {
            const gap = radii[i] - radii[i - 1];
            if (gap > maxGap) maxGap = gap;
        }
        if (maxGap > thickness) ++violating;
    }
    return populated > 0 ? violating / populated : 0;
}

/** One star-shaped lobe: a radial SH expansion about its own center. */
export interface SphericalHarmonicLobe { center: number[], coeffs: Float64Array, rMax: number }
export interface SphericalHarmonicLobesFit { lobes: SphericalHarmonicLobe[], L: number }

export interface FitSphericalHarmonicLobesOptions {
    /** Maximum number of lobes to split into (1 = single envelope). */
    maxLobes?: number
    /** Multi-valued-bin fraction above which a lobe is split (see `starShapeViolation`). */
    tolerance?: number
    /** Radial-extent threshold (length units) for the violation metric. */
    thickness?: number
    /** Sample cap forwarded to `fitSphericalHarmonics`. */
    maxPoints?: number
    /** Tikhonov regularization forwarded to `fitSphericalHarmonics`. */
    regularization?: number
}

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

/** Bisect a point subset by the sign of its projection onto its longest principal axis. */
function bisectSubset(flat: Float64Array, idx: ArrayLike<number>): [Int32Array, Int32Array] | undefined {
    const m = idx.length;
    if (m < 2) return undefined;
    const axes = PrincipalAxes.calculateMomentsAxes(flat);
    const ox = axes.origin[0], oy = axes.origin[1], oz = axes.origin[2];
    let ax = axes.dirA[0], ay = axes.dirA[1], az = axes.dirA[2];
    const len = Math.sqrt(ax * ax + ay * ay + az * az);
    if (len < 1e-9) return undefined;
    ax /= len; ay /= len; az /= len;

    const left: number[] = [], right: number[] = [];
    for (let i = 0; i < m; ++i) {
        const p = i * 3;
        const proj = (flat[p] - ox) * ax + (flat[p + 1] - oy) * ay + (flat[p + 2] - oz) * az;
        if (proj < 0) left.push(idx[i]); else right.push(idx[i]);
    }
    if (left.length === 0 || right.length === 0) return undefined;
    return [Int32Array.from(left), Int32Array.from(right)];
}

/**
 * Fit a point cloud to one or more star-shaped radial SH lobes.
 *
 * A single-center radial expansion only represents star-shaped (radially
 * single-valued) shapes; concave, elongated or multi-domain clouds fold onto the
 * outermost crossing. When the directional multi-valuedness (`starShapeViolation`)
 * of a lobe exceeds `tolerance`, the lobe is bisected along its longest principal
 * axis and each half refit, repeated until every lobe is star-shaped enough or
 * `maxLobes` is reached.
 */
export function fitSphericalHarmonicLobes(points: ArrayLike<number>, L: number, options: FitSphericalHarmonicLobesOptions = {}): SphericalHarmonicLobesFit {
    const maxLobes = Math.max(1, Math.round(options.maxLobes ?? 1));
    const tolerance = options.tolerance ?? 0.15;
    const thickness = options.thickness ?? 2;
    const maxPoints = options.maxPoints ?? 8192;
    const regularization = options.regularization ?? 0;

    const n = Math.floor(points.length / 3);

    type Work = { idx: Int32Array, flat: Float64Array, center: number[], coeffs: Float64Array, rMax: number, violation: number };
    const buildLobe = (idx: Int32Array): Work => {
        const { flat, center } = gatherSubset(points, idx);
        const { coeffs, rMax } = fitSphericalHarmonics(flat, center, L, maxPoints, regularization);
        const violation = starShapeViolation(flat, center, thickness);
        return { idx, flat, center, coeffs, rMax, violation };
    };

    const all = new Int32Array(n);
    for (let i = 0; i < n; ++i) all[i] = i;
    const lobes: Work[] = n > 0 ? [buildLobe(all)] : [];

    while (lobes.length > 0 && lobes.length < maxLobes) {
        // pick the worst star-shape violator that is still worth splitting
        let wi = -1, wv = tolerance;
        for (let i = 0; i < lobes.length; ++i) {
            if (lobes[i].violation > wv && lobes[i].idx.length >= 4) { wv = lobes[i].violation; wi = i; }
        }
        if (wi < 0) break;

        const split = bisectSubset(lobes[wi].flat, lobes[wi].idx);
        if (!split) { lobes[wi].violation = 0; continue; } // can't split this one; leave it
        lobes.splice(wi, 1, buildLobe(split[0]), buildLobe(split[1]));
    }

    return { lobes: lobes.map(l => ({ center: l.center, coeffs: l.coeffs, rMax: l.rMax })), L };
}
