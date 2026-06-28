/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { calculateHistogram, downsampleHistogram, histogramPercentile, histogramRobustStats } from '../histogram';

describe('histogram', () => {
    it('histogramPercentile is approximately linear for a uniform distribution', () => {
        const n = 10000;
        const data = new Float32Array(n);
        for (let i = 0; i < n; ++i) data[i] = i / n; // uniform on [0, 1)
        const h = calculateHistogram(data, 100, { min: 0, max: 1 });
        expect(histogramPercentile(h, 0.5)).toBeCloseTo(0.5, 1);
        expect(histogramPercentile(h, 0.01)).toBeCloseTo(0.01, 1);
        expect(histogramPercentile(h, 0.99)).toBeCloseTo(0.99, 1);
    });

    it('histogramRobustStats can ignore the zero bin', () => {
        // Half zeros, half values uniformly in [1, 2].
        const data: number[] = [];
        for (let i = 0; i < 5000; ++i) data.push(0);
        for (let i = 0; i < 5000; ++i) data.push(1 + i / 5000);
        const h = calculateHistogram(data, 200, { min: 0, max: 2 });

        const all = histogramRobustStats(h);
        const noZero = histogramRobustStats(h, { ignoreZero: true });
        // Without filtering, mean is dragged toward 0.
        expect(all.mean).toBeLessThan(noZero.mean);
        // With zero excluded, mean lands near 1.5.
        expect(noZero.mean).toBeCloseTo(1.5, 1);
    });

    it('histogramPercentile handles a signed, symmetric distribution', () => {
        const n = 10000;
        const data = new Float32Array(n);
        for (let i = 0; i < n; ++i) data[i] = -1 + 2 * i / n; // uniform on [-1, 1)
        const h = calculateHistogram(data, 200, { min: -1, max: 1 });
        expect(histogramPercentile(h, 0.5)).toBeCloseTo(0, 1);
        expect(histogramPercentile(h, 0.25)).toBeCloseTo(-0.5, 1);
        expect(histogramPercentile(h, 0.75)).toBeCloseTo(0.5, 1);
    });

    it('histogramPercentile is stable for an empty histogram', () => {
        const h = calculateHistogram(new Float32Array(0), 16, { min: 0, max: 1 });
        // No counts: falls back to the lower bound.
        expect(histogramPercentile(h, 0.5)).toBe(h.min);
    });

    it('histogramPercentile ignores an out-of-range excludeBin', () => {
        const n = 1000;
        const data = new Float32Array(n);
        for (let i = 0; i < n; ++i) data[i] = i / n;
        const h = calculateHistogram(data, 50, { min: 0, max: 1 });
        const ref = histogramPercentile(h, 0.5);
        // excludeBin beyond the last index must not change the result or throw.
        expect(histogramPercentile(h, 0.5, 9999)).toBeCloseTo(ref, 6);
    });

    it('histogramRobustStats keeps the zero bin for signed (difference-map-like) data', () => {
        // Symmetric around 0: the zero bin is the mode and must not be dropped.
        const n = 10000;
        const data = new Float32Array(n);
        for (let i = 0; i < n; ++i) data[i] = -1 + 2 * i / n;
        const h = calculateHistogram(data, 200, { min: -1, max: 1 });
        const all = histogramRobustStats(h);
        const noZero = histogramRobustStats(h, { ignoreZero: true });
        // ignoreZero only drops the first bin, so a centred zero is untouched.
        expect(noZero.mean).toBeCloseTo(all.mean, 6);
        expect(noZero.count).toBe(all.count);
    });

    it('downsampleHistogram preserves the total count and reduces bin count', () => {
        const n = 4096;
        const data = new Float32Array(n);
        for (let i = 0; i < n; ++i) data[i] = i / n;
        const fine = calculateHistogram(data, 1024, { min: 0, max: 1 });
        const coarse = downsampleHistogram(fine, 40);
        const sum = (a: ArrayLike<number>) => { let s = 0; for (let i = 0; i < a.length; ++i) s += a[i]; return s; };
        expect(coarse.counts.length).toBe(40);
        expect(sum(coarse.counts)).toBe(sum(fine.counts));
        // Returns the input unchanged when the target is not coarser.
        expect(downsampleHistogram(fine, 1024)).toBe(fine);
    });
});
