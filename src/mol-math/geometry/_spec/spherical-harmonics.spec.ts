/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { assocLegendre, realSph, shIndex, shTermCount, fitSphericalHarmonics, reconstructRadius, buildRadiusLUT, sampleRadiusLUT, getNormFactors } from '../spherical-harmonics';

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
        const sphOut = new Float64Array(K);
        const legendreScratch = new Float64Array((L + 1) * (L + 2) / 2);
        const norm = getNormFactors(L);
        // integrate Y_i Y_j sin(theta) dtheta dphi over the sphere
        for (let it = 0; it < nTheta; ++it) {
            const theta = (it + 0.5) / nTheta * Math.PI;
            const w = Math.sin(theta) * (Math.PI / nTheta) * (2 * Math.PI / nPhi);
            for (let ip = 0; ip < nPhi; ++ip) {
                const phi = (ip + 0.5) / nPhi * 2 * Math.PI;
                const y = realSph(L, theta, phi, sphOut, legendreScratch, norm);
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

        const sphOut = new Float64Array(K);
        const legendreScratch = new Float64Array((L + 1) * (L + 2) / 2);
        const norm = getNormFactors(L);

        // scatter sample directions (deterministic, well spread) and set r = sum c_i Y_i
        const n = 2000;
        const pts = new Float32Array(n * 3);
        const golden = Math.PI * (3 - Math.sqrt(5));
        for (let i = 0; i < n; ++i) {
            const z = 1 - (2 * i + 1) / n;
            const theta = Math.acos(z);
            const phi = i * golden;
            const r = reconstructRadius(truth, L, theta, phi, sphOut, legendreScratch, norm);
            const s = Math.sin(theta);
            pts[i * 3] = r * s * Math.cos(phi);
            pts[i * 3 + 1] = r * s * Math.sin(phi);
            pts[i * 3 + 2] = r * Math.cos(theta);
        }

        const { coeffs } = fitSphericalHarmonics(pts, [0, 0, 0], L);
        for (let i = 0; i < K; ++i) {
            expect(coeffs[i]).toBeCloseTo(truth[i], 3);
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

        const sphOut = new Float64Array(K);
        const legendreScratch = new Float64Array((L + 1) * (L + 2) / 2);
        const norm = getNormFactors(L);

        const n = 2000;
        const pts = new Float32Array(n * 3);
        const golden = Math.PI * (3 - Math.sqrt(5));
        for (let i = 0; i < n; ++i) {
            const z = 1 - (2 * i + 1) / n;
            const theta = Math.acos(z);
            const phi = i * golden;
            const r = reconstructRadius(truth, L, theta, phi, sphOut, legendreScratch, norm);
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
            expect(reconstructRadius(coeffs, L, theta, phi, sphOut, legendreScratch, norm)).toBeCloseTo(reconstructRadius(truth, L, theta, phi, sphOut, legendreScratch, norm), 1);
        }
    });

    it('regularization tames an under-determined (sparse, clustered) fit', () => {
        const L = 8; // K = 81 coefficients
        const K = shTermCount(L);
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
            const sphOut = new Float64Array(K);
            const legendreScratch = new Float64Array((L + 1) * (L + 2) / 2);
            const norm = getNormFactors(L);
            for (let i = 0; i < 2000; ++i) {
                const th = Math.acos(1 - (2 * i + 1) / 2000);
                m = Math.max(m, reconstructRadius(coeffs, L, th, i * golden, sphOut, legendreScratch, norm));
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

    it('RadiusLUT bilinearly approximates reconstructRadius closely, including at poles/phi wraparound', () => {
        const L = 4;
        const K = shTermCount(L);
        const truth = new Float64Array(K);
        truth[shIndex(0, 0)] = 10;
        truth[shIndex(1, -1)] = 1.3;
        truth[shIndex(1, 0)] = -0.7;
        truth[shIndex(2, 2)] = 0.9;
        truth[shIndex(3, -2)] = -0.4;
        truth[shIndex(4, 1)] = 0.6;

        const nTheta = 6 * (L + 1), nPhi = 12 * (L + 1);
        const lut = buildRadiusLUT(truth, L, nTheta, nPhi);

        const sphOut = new Float64Array(K);
        const legendreScratch = new Float64Array((L + 1) * (L + 2) / 2);
        const norm = getNormFactors(L);

        const golden = Math.PI * (3 - Math.sqrt(5));
        for (let i = 0; i < 500; ++i) {
            const theta = Math.acos(1 - (2 * i + 1) / 500);
            const phi = ((i * golden + Math.PI) % (2 * Math.PI)) - Math.PI;
            expect(sampleRadiusLUT(lut, theta, phi)).toBeCloseTo(reconstructRadius(truth, L, theta, phi, sphOut, legendreScratch, norm), 1);
        }

        // poles (theta = 0 and theta = pi) and phi wraparound (phi = -pi vs phi = pi) must not throw
        // or produce wildly different values for what are effectively the same/adjacent directions
        expect(sampleRadiusLUT(lut, 0, 0)).toBeCloseTo(reconstructRadius(truth, L, 0, 0, sphOut, legendreScratch, norm), 1);
        expect(sampleRadiusLUT(lut, Math.PI, 0)).toBeCloseTo(reconstructRadius(truth, L, Math.PI, 0, sphOut, legendreScratch, norm), 1);
        expect(sampleRadiusLUT(lut, Math.PI / 2, -Math.PI)).toBeCloseTo(sampleRadiusLUT(lut, Math.PI / 2, Math.PI), 6);
    });

    it('RadiusLUT reuses the provided out buffer when dimensions match', () => {
        const L = 2;
        const K = shTermCount(L);
        const coeffsA = new Float64Array(shTermCount(L)); coeffsA[0] = 5;
        const coeffsB = new Float64Array(shTermCount(L)); coeffsB[0] = 8;
        const nTheta = 8, nPhi = 16;

        const lutA = buildRadiusLUT(coeffsA, L, nTheta, nPhi);
        const lutB = buildRadiusLUT(coeffsB, L, nTheta, nPhi, -Infinity, Infinity, lutA);

        const sphOut = new Float64Array(K);
        const legendreScratch = new Float64Array((L + 1) * (L + 2) / 2);
        const norm = getNormFactors(L);

        expect(lutB.values).toBe(lutA.values); // buffer reused, not reallocated
        expect(sampleRadiusLUT(lutB, Math.PI / 2, 0)).toBeCloseTo(reconstructRadius(coeffsB, L, Math.PI / 2, 0, sphOut, legendreScratch, norm), 6);
    });
});