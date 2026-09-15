/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { squaredDistanceTransform3D } from '../distance-transform';

function bruteForce(inside: Uint8Array, nx: number, ny: number, nz: number) {
    const out = new Float32Array(nx * ny * nz);
    const pts: number[][] = [];
    for (let z = 0; z < nz; z++) for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
        if (inside[x + y * nx + z * nx * ny]) pts.push([x, y, z]);
    }
    for (let z = 0; z < nz; z++) for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
        let best = Infinity;
        for (const [px, py, pz] of pts) {
            const d = (x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2;
            if (d < best) best = d;
        }
        out[x + y * nx + z * nx * ny] = best;
    }
    return out;
}

/** Small deterministic PRNG for reproducible random fixtures. */
function lcg(seed: number) {
    let s = seed >>> 0;
    return () => {
        s = (s * 1664525 + 1013904223) >>> 0;
        return s / 0x100000000;
    };
}

describe('squaredDistanceTransform3D', () => {
    it('matches brute force on random masks', () => {
        const [nx, ny, nz] = [10, 12, 9];
        const rnd = lcg(42);
        for (let trial = 0; trial < 5; trial++) {
            const inside = new Uint8Array(nx * ny * nz);
            for (let i = 0; i < inside.length; i++) inside[i] = rnd() < 0.05 ? 1 : 0;
            if (!inside.some(v => v)) inside[7] = 1;
            const expected = bruteForce(inside, nx, ny, nz);
            const actual = squaredDistanceTransform3D(inside, nx, ny, nz);
            for (let i = 0; i < inside.length; i++) expect(actual[i]).toBeCloseTo(expected[i], 5);
        }
    });

    it('is zero inside and huge when nothing is inside', () => {
        const full = new Uint8Array(27).fill(1);
        expect(Array.from(squaredDistanceTransform3D(full, 3, 3, 3))).toEqual(new Array(27).fill(0));
        const empty = new Uint8Array(27);
        for (const v of squaredDistanceTransform3D(empty, 3, 3, 3)) expect(v).toBeGreaterThanOrEqual(1e20);
    });

    it('reuses the provided output buffer', () => {
        const inside = new Uint8Array(8);
        inside[0] = 1;
        const out = new Float32Array(8);
        expect(squaredDistanceTransform3D(inside, 2, 2, 2, out)).toBe(out);
        expect(out[7]).toBe(3);
    });
});
