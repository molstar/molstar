/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Histogram, ValueStats } from '../histogram.js';

// A density map's values for the stories, generated (no fetch): an electron density map is mostly solvent (a narrow
// peak around zero) with the molecule's density a long tail to the right. Seeded, so every render is the same.

function random(seed: number) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** The values of a map: `n` voxels, a fraction of them the molecule's (the right tail). */
export function densityValues(n = 120_000, molecule = 0.3, seed = 7): Float64Array {
  const rnd = random(seed);
  const normal = () => Math.sqrt(-2 * Math.log(1 - rnd())) * Math.cos(2 * Math.PI * rnd());
  const out = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    out[i] = rnd() < molecule ? 0.1 + Math.exp(-1.2 + 0.55 * normal()) : 0.08 * normal();
  }
  return out;
}

/** The statistics and the histogram (`bins` equal bins over [min, max]) of a map's values. */
export function densityData(
  values: Float64Array = densityValues(),
  bins = 160,
): { stats: ValueStats; histogram: Histogram } {
  let min = Infinity;
  let max = -Infinity;
  let sum = 0;
  for (const v of values) {
    if (v < min) min = v;
    if (v > max) max = v;
    sum += v;
  }
  const mean = sum / values.length;
  let sq = 0;
  for (const v of values) sq += (v - mean) ** 2;
  const sigma = Math.sqrt(sq / values.length);
  const counts = new Uint32Array(bins);
  const width = (max - min) / bins || 1;
  for (const v of values) counts[Math.min(bins - 1, Math.floor((v - min) / width))]++;
  return { stats: { min, max, mean, sigma }, histogram: { min, max, counts } };
}
