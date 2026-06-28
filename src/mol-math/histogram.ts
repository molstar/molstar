/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { arrayMinMax } from '../mol-util/array';

export interface Histogram {
    min: number,
    max: number,
    binWidth: number,
    counts: Int32Array
}

export interface HistogramRobustStats {
    /** Total voxel count contributing to the stats (after optional bin exclusion). */
    count: number,
    /** Minimum and maximum value over the contributing range. */
    min: number,
    max: number,
    /** Mean and population standard deviation, recomputed from bin centers. */
    mean: number,
    sigma: number,
    /** Approximate percentiles (linear interpolation across bin counts). */
    p1: number,
    p50: number,
    p95: number,
    p99: number,
}

export function calculateHistogram(data: ArrayLike<number>, binCount: number, options?: { min: number, max: number, }): Histogram {
    if (!options) {
        const [min, max] = arrayMinMax(data);
        return _calcHistogram(data, binCount, min, max);
    } else {
        return _calcHistogram(data, binCount, options.min, options.max);
    }
}

function _calcHistogram(data: ArrayLike<number>, binCount: number, min: number, max: number): Histogram {
    let binWidth = (max - min) / binCount;
    if (binWidth === 0) binWidth = 1;

    const counts = new Int32Array(binCount);

    for (let i = 0, _i = data.length; i < _i; i++) {
        let bin = Math.floor((data[i] - min) / binWidth);
        if (bin >= binCount) bin = binCount - 1;
        else if (bin < 0) bin = 0;
        counts[bin]++;
    }

    return { min, max, binWidth, counts };
}

/**
 * Approximate `p`-th percentile (p in [0, 1]) of a histogram via linear
 * interpolation across the cumulative bin counts. Bin centers are used as
 * value samples. If `excludeBin` is provided, that bin is treated as having
 * zero counts (useful to ignore an explicit "background" bin, e.g. the bin
 * containing exact zeros for masked density maps).
 */
export function histogramPercentile(h: Histogram, p: number, excludeBin?: number): number {
    const { counts, min, binWidth } = h;
    const n = counts.length;
    if (n === 0) return min;

    let total = 0;
    for (let i = 0; i < n; ++i) {
        if (i === excludeBin) continue;
        total += counts[i];
    }
    if (total === 0) return min;

    const target = Math.max(0, Math.min(1, p)) * total;

    let cum = 0;
    for (let i = 0; i < n; ++i) {
        if (i === excludeBin) continue;
        const c = counts[i];
        if (c === 0) continue;
        const next = cum + c;
        if (next >= target) {
            // Linear interpolation within the bin.
            const t = c > 0 ? (target - cum) / c : 0;
            return min + (i + t) * binWidth;
        }
        cum = next;
    }
    return min + n * binWidth;
}

/**
 * Compute robust summary statistics directly from a `Histogram`.
 * Mean and sigma are recomputed from bin centers so they stay consistent
 * with the percentiles when `ignoreZero` excludes the bin containing 0.
 *
 * Note: "robust" here means robust against a dominant zero/background bin,
 * not against outliers. The histogram is still binned across the full
 * `[min, max]` range, so a single extreme value coarsens every bin.
 */
export function histogramRobustStats(h: Histogram, opts?: { ignoreZero?: boolean }): HistogramRobustStats {
    const { counts, min, binWidth } = h;
    const n = counts.length;

    let excludeBin: number | undefined;
    if (opts?.ignoreZero && binWidth > 0) {
        // Only treat zero as background when it sits in the first bin (masked/
        // padded maps where 0 is the minimum). For signed maps where 0 is the
        // centre of the distribution (e.g. difference maps) the zero bin is the
        // mode and must not be removed, so leave it in place.
        const b = Math.floor((0 - min) / binWidth);
        if (b === 0) excludeBin = 0;
    }

    let count = 0;
    let sum = 0;
    let sumSq = 0;
    let firstBin = -1;
    let lastBin = -1;
    for (let i = 0; i < n; ++i) {
        if (i === excludeBin) continue;
        const c = counts[i];
        if (c === 0) continue;
        const center = min + (i + 0.5) * binWidth;
        count += c;
        sum += c * center;
        sumSq += c * center * center;
        if (firstBin < 0) firstBin = i;
        lastBin = i;
    }

    if (count === 0) {
        return { count: 0, min, max: min, mean: min, sigma: 0, p1: min, p50: min, p95: min, p99: min };
    }

    const mean = sum / count;
    const variance = Math.max(0, sumSq / count - mean * mean);
    const sigma = Math.sqrt(variance);
    const lo = min + firstBin * binWidth;
    const hi = min + (lastBin + 1) * binWidth;

    return {
        count,
        min: lo,
        max: hi,
        mean,
        sigma,
        p1: histogramPercentile(h, 0.01, excludeBin),
        p50: histogramPercentile(h, 0.50, excludeBin),
        p95: histogramPercentile(h, 0.95, excludeBin),
        p99: histogramPercentile(h, 0.99, excludeBin),
    };
}

/**
 * Regroup a fine-grained histogram into a coarser one with `binCount` bins
 * over the same `[min, max]` range. Counts are preserved (each source bin is
 * mapped to a target bin by index), letting callers reuse an existing
 * histogram for display instead of re-scanning the source data. When
 * `binCount >= h.counts.length` the original histogram is returned unchanged.
 */
export function downsampleHistogram(h: Histogram, binCount: number): Histogram {
    const srcN = h.counts.length;
    if (binCount <= 0 || binCount >= srcN) return h;

    const counts = new Int32Array(binCount);
    for (let i = 0; i < srcN; ++i) {
        let bin = Math.floor(i * binCount / srcN);
        if (bin >= binCount) bin = binCount - 1;
        counts[bin] += h.counts[i];
    }

    const binWidth = (h.max - h.min) / binCount || 1;
    return { min: h.min, max: h.max, binWidth, counts };
}