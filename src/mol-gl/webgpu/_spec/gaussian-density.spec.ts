/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { createGaussianDensityBins } from '../gaussian-density';

describe('WebGPU Gaussian spatial bins', () => {
    const atoms = new Float32Array(32);
    for (let i = 0; i < 4; i++) atoms.set([i * 30 + 5, 5, 5, 1], i * 8);
    const dimensions = [100, 10, 10], origin = [0, 0, 0];

    it('retains every contributing atom in input order while reducing evaluations', () => {
        const bins = createGaussianDensityBins(atoms, 4, dimensions, origin, 1, 1 << 20);
        expect(bins.evaluations).toBeLessThan(100 * 10 * 10 * 4 / 4);
        for (let x = 0; x < 100; x++) for (let y = 0; y < 10; y++) for (let z = 0; z < 10; z++) {
            const cell = (Math.floor(x / bins.span) * bins.dimensions[1] + Math.floor(y / bins.span)) * bins.dimensions[2] + Math.floor(z / bins.span);
            const offset = bins.ranges[cell * 2], length = bins.ranges[cell * 2 + 1];
            const candidates = Array.from(bins.indices.subarray(offset, offset + length));
            expect(candidates).toEqual([...candidates].sort((a, b) => a - b));
            for (let i = 0; i < 4; i++) {
                if ((x - atoms[i * 8]) ** 2 + (y - 5) ** 2 + (z - 5) ** 2 <= 4) expect(candidates).toContain(i);
            }
        }
    });

    it('coarsens bins to respect storage limits without dropping atoms', () => {
        const bins = createGaussianDensityBins(atoms, 4, dimensions, origin, 1, 32);
        expect(bins.ranges.byteLength).toBeLessThanOrEqual(32);
        expect(bins.indices.byteLength).toBeLessThanOrEqual(32);
        expect(Array.from(bins.indices)).toEqual(expect.arrayContaining([0, 1, 2, 3]));
        const full = createGaussianDensityBins(atoms, 4, dimensions, origin, 1, 32, false);
        expect(Array.from(full.indices)).toEqual([0, 1, 2, 3]);
        expect(full.evaluations).toBe(40000);
    });
});
