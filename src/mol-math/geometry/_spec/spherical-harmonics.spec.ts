/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { assocLegendre, realSph, shIndex, shTermCount, fitSphericalHarmonics, reconstructRadius, fitSphericalHarmonicLobesByLabel, kmeansLabels } from '../spherical-harmonics';

/** Fibonacci-sphere surface points of radius R about (cx, cy, cz), written into `out` from `offset`. */
function sphereCloud(n: number, R: number, cx: number, cy: number, cz: number, out: Float32Array, offset = 0) {
    const golden = Math.PI * (3 - Math.sqrt(5));
    for (let i = 0; i < n; ++i) {
        const z = 1 - (2 * i + 1) / n;
        const theta = Math.acos(z);
        const phi = i * golden;
        const s = Math.sin(theta);
        const o = offset + i * 3;
        out[o] = cx + R * s * Math.cos(phi);
        out[o + 1] = cy + R * s * Math.sin(phi);
        out[o + 2] = cz + R * Math.cos(theta);
    }
}

describe('spherical-harmonics', () => {
    it('associated Legendre values match known closed forms', () => {
        const x = 0.3;
        const s = Math.sqrt(1 - x * x); // sin(theta)
        const p = assocLegendre(3, x);
        const idx = (l: number, m: number) => (l * (l + 1)) / 2 + m;
        // P_0^0 = 1, P_1^0 = x, P_1^1 = -sin, P_2^0 = (3x^2-1)/2, P_2^2 = 3 sin^2
        expect(p[idx(0, 0)]).toBeCloseTo(1, 10);
        expect(p[idx(1, 0)]).toBeCloseTo(x, 10);
        expect(p[idx(1, 1)]).toBeCloseTo(-s, 10);
        expect(p[idx(2, 0)]).toBeCloseTo((3 * x * x - 1) / 2, 10);
        expect(p[idx(2, 2)]).toBeCloseTo(3 * s * s, 10);
    });

    it('real SH basis is orthonormal over the sphere (low degree)', () => {
        const L = 3;
        const K = shTermCount(L);
        const nTheta = 64, nPhi = 128;
        const acc = new Float64Array(K * K);
        // integrate Y_i Y_j sin(theta) dtheta dphi over the sphere
        for (let it = 0; it < nTheta; ++it) {
            const theta = (it + 0.5) / nTheta * Math.PI;
            const w = Math.sin(theta) * (Math.PI / nTheta) * (2 * Math.PI / nPhi);
            for (let ip = 0; ip < nPhi; ++ip) {
                const phi = (ip + 0.5) / nPhi * 2 * Math.PI;
                const y = realSph(L, theta, phi);
                for (let a = 0; a < K; ++a) {
                    for (let b = 0; b < K; ++b) acc[a * K + b] += y[a] * y[b] * w;
                }
            }
        }
        for (let a = 0; a < K; ++a) {
            for (let b = 0; b < K; ++b) {
                expect(acc[a * K + b]).toBeCloseTo(a === b ? 1 : 0, 2);
            }
        }
    });

    it('round-trips known coefficients through fit (catches SVD-solve layout bugs)', () => {
        const L = 4;
        const K = shTermCount(L);
        // arbitrary non-trivial coefficients across several (l, m) terms
        const truth = new Float64Array(K);
        truth[shIndex(0, 0)] = 10;
        truth[shIndex(1, -1)] = 1.3;
        truth[shIndex(1, 0)] = -0.7;
        truth[shIndex(2, 2)] = 0.9;
        truth[shIndex(3, -2)] = -0.4;
        truth[shIndex(4, 1)] = 0.6;

        // scatter sample directions (deterministic, well spread) and set r = sum c_i Y_i
        const n = 2000;
        const pts = new Float32Array(n * 3);
        const golden = Math.PI * (3 - Math.sqrt(5));
        for (let i = 0; i < n; ++i) {
            const z = 1 - (2 * i + 1) / n;
            const theta = Math.acos(z);
            const phi = i * golden;
            const r = reconstructRadius(truth, L, theta, phi);
            const s = Math.sin(theta);
            pts[i * 3] = r * s * Math.cos(phi);
            pts[i * 3 + 1] = r * s * Math.sin(phi);
            pts[i * 3 + 2] = r * Math.cos(theta);
        }

        const { coeffs } = fitSphericalHarmonics(pts, [0, 0, 0], L);
        for (let i = 0; i < K; ++i) {
            expect(coeffs[i]).toBeCloseTo(truth[i], 4);
        }
    });

    it('regularized fit (Cholesky path) reconstructs a known low-degree shape accurately', () => {
        const L = 4;
        const K = shTermCount(L);
        // a low-degree, well-conditioned truth shape (dominant l=0 keeps the radius strictly positive)
        const truth = new Float64Array(K);
        truth[shIndex(0, 0)] = 10;
        truth[shIndex(1, 0)] = -0.7;
        truth[shIndex(2, 2)] = 0.9;
        truth[shIndex(3, -2)] = -0.4;

        const n = 2000;
        const pts = new Float32Array(n * 3);
        const golden = Math.PI * (3 - Math.sqrt(5));
        for (let i = 0; i < n; ++i) {
            const z = 1 - (2 * i + 1) / n;
            const theta = Math.acos(z);
            const phi = i * golden;
            const r = reconstructRadius(truth, L, theta, phi);
            const s = Math.sin(theta);
            pts[i * 3] = r * s * Math.cos(phi);
            pts[i * 3 + 1] = r * s * Math.sin(phi);
            pts[i * 3 + 2] = r * Math.cos(theta);
        }

        // regularization > 0 routes through the Cholesky solve; with dense samples and light damping
        // the reconstruction tracks the truth radius closely at arbitrary directions
        const { coeffs } = fitSphericalHarmonics(pts, [0, 0, 0], L, undefined, 0.001);
        for (let i = 0; i < 50; ++i) {
            const theta = Math.acos(1 - (2 * i + 1) / 50);
            const phi = i * golden;
            expect(reconstructRadius(coeffs, L, theta, phi)).toBeCloseTo(reconstructRadius(truth, L, theta, phi), 1);
        }
    });

    it('regularization tames an under-determined (sparse, clustered) fit', () => {
        const L = 8; // K = 81 coefficients
        // 20 points strung along x with jitter: far fewer than K, ill-conditioned
        const n = 20;
        const pts = new Float32Array(n * 3);
        let seed = 1;
        const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647 - 0.5; };
        let dataMax = 0;
        for (let i = 0; i < n; ++i) {
            const x = ((i / (n - 1)) * 2 - 1) * 30 + rnd() * 6;
            const y = rnd() * 6, z = rnd() * 6;
            pts[i * 3] = x; pts[i * 3 + 1] = y; pts[i * 3 + 2] = z;
            dataMax = Math.max(dataMax, Math.hypot(x, y, z));
        }

        const maxReconstructed = (coeffs: Float64Array) => {
            let m = -Infinity;
            const golden = Math.PI * (3 - Math.sqrt(5));
            for (let i = 0; i < 2000; ++i) {
                const th = Math.acos(1 - (2 * i + 1) / 2000);
                m = Math.max(m, reconstructRadius(coeffs, L, th, i * golden));
            }
            return m;
        };

        // unregularized: blows up far beyond the data extent
        const bare = fitSphericalHarmonics(pts, [0, 0, 0], L, undefined, 0);
        expect(maxReconstructed(bare.coeffs)).toBeGreaterThan(dataMax * 5);

        // regularized: stays within a small factor of the data extent
        const reg = fitSphericalHarmonics(pts, [0, 0, 0], L, undefined, 0.05);
        expect(maxReconstructed(reg.coeffs)).toBeLessThan(dataMax * 2);
        expect(reg.rMax).toBeCloseTo(dataMax, 5);
    });

    it('kmeansLabels separates a dumbbell into two compact clusters (deterministic)', () => {
        const n = 2000;
        // two R=8 spheres centered at x = +-20
        const dumbbell = new Float32Array(2 * n * 3);
        sphereCloud(n, 8, -20, 0, 0, dumbbell, 0);
        sphereCloud(n, 8, 20, 0, 0, dumbbell, n * 3);

        const labels = kmeansLabels(dumbbell, 2, { seed: 1 });
        // every point of one input sphere shares a label, distinct from the other sphere's label
        const a = labels[0], b = labels[2 * n - 1];
        expect(a).not.toBe(b);
        for (let i = 0; i < n; ++i) expect(labels[i]).toBe(a);
        for (let i = n; i < 2 * n; ++i) expect(labels[i]).toBe(b);

        // deterministic for a fixed seed
        const labels2 = kmeansLabels(dumbbell, 2, { seed: 1 });
        for (let i = 0; i < 2 * n; ++i) expect(labels2[i]).toBe(labels[i]);
    });

    it('fitSphericalHarmonicLobesByLabel on k-means clusters gives compact per-lobe envelopes', () => {
        const n = 2000;
        const L = 6;
        const dumbbell = new Float32Array(2 * n * 3);
        sphereCloud(n, 8, -20, 0, 0, dumbbell, 0);
        sphereCloud(n, 8, 20, 0, 0, dumbbell, n * 3);

        const labels = kmeansLabels(dumbbell, 2, { seed: 1 });
        const { lobes } = fitSphericalHarmonicLobesByLabel(dumbbell, labels, L, { regularization: 0.01 });
        expect(lobes.length).toBe(2);
        // each lobe is the size of one R=8 sphere, not the whole ~48 A span
        for (const lobe of lobes) {
            expect(Math.abs(lobe.center[0])).toBeGreaterThan(10);
            expect(lobe.rMax).toBeLessThan(16);
        }
    });

    it('radii inflation keeps the lobe centered (center independent of offset, rMax grows with it)', () => {
        const n = 2000;
        const L = 6;
        // a single off-origin sphere; an asymmetric cloud would drift its center if inflated in 3D first
        const sphere = new Float32Array(n * 3);
        sphereCloud(n, 10, 30, -5, 12, sphere);
        const labels = new Int32Array(n); // one lobe

        const base = fitSphericalHarmonicLobesByLabel(sphere, labels, L, { regularization: 0.01 }).lobes[0];
        const inflated = fitSphericalHarmonicLobesByLabel(sphere, labels, L, { regularization: 0.01, radii: new Float32Array(n).fill(5) }).lobes[0];
        // center unchanged by the 5 A inflation
        for (let d = 0; d < 3; ++d) expect(inflated.center[d]).toBeCloseTo(base.center[d], 6);
        // surface grows by ~the offset
        expect(inflated.rMax - base.rMax).toBeGreaterThan(4);
        expect(inflated.rMax - base.rMax).toBeLessThan(6);
    });

    it('fitSphericalHarmonicLobesByLabel fits one lobe per label, centered on its subset', () => {
        const n = 1000;
        const L = 6;

        // two separated spheres, one label each
        const points = new Float32Array(2 * n * 3);
        sphereCloud(n, 8, -20, 0, 0, points, 0);
        sphereCloud(n, 8, 20, 0, 0, points, n * 3);
        const labels = new Int32Array(2 * n);
        for (let i = n; i < 2 * n; ++i) labels[i] = 1;

        const fit = fitSphericalHarmonicLobesByLabel(points, labels, L);
        expect(fit.lobes.length).toBe(2);
        // labels are partitioned by Map insertion order: lobe 0 = label 0 (x ~ -20), lobe 1 = label 1 (x ~ +20)
        expect(fit.lobes[0].center[0]).toBeLessThan(-10);
        expect(fit.lobes[1].center[0]).toBeGreaterThan(10);
    });
});
