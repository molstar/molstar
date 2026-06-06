/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 */

import { assocLegendre, realSphAdvanced, shIndex, shTermCount, fitSphericalHarmonics, reconstructRadius } from '../spherical-harmonics';

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
                const y = realSphAdvanced(L, theta, phi);
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
});
