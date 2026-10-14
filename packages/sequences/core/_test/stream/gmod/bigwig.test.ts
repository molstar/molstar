/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// BigWig against the recorded slices of UCSC's gc5Base.bw (GC percent of 5-base windows). The expected values are the
// REST API's (getData/track, gc5BaseBw) over the same regions, recorded in fixtures/expected.json: the raw windows for a
// data request, the same windows folded onto the level's bins for a level request.

import { readFileSync } from 'node:fs';
import { StreamNotServed } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { CASES, makers, req } from './cases.js';
import { sliceReader } from './slices.js';

const expected = JSON.parse(readFileSync(new URL('./fixtures/expected.json', import.meta.url), 'utf8'))
  .bigwig as Record<string, unknown>;
const signal = new AbortController().signal;
const provider = () => makers.bigwig(sliceReader('bigwig'));

describe('BigWig', () => {
  it("lists the chromosomes, the zoom levels as levels and the file's range as the domain", async () => {
    const p = await provider().fetch(req(null), signal);
    if (p.kind !== 'manifest') throw new Error('not a manifest');
    const sizes = Object.fromEntries(p.manifest.systems.map((s) => [s.key, s.extent[1]]));
    expect([sizes.chr1, sizes.chr21, sizes.chrM]).toEqual([248_956_422, 46_709_983, 16_569]);
    expect(p.manifest.systems.length).toBe(595);
    expect(p.manifest.levels).toEqual([
      100, 400, 1_600, 6_400, 25_600, 102_400, 409_600, 1_638_400, 6_553_600, 26_214_400,
    ]);
    expect(p.manifest.domain).toEqual([0, 100]);
  });

  for (const [name, r] of CASES.bigwig.slice(1).filter(([, r]) => r.level < 0))
    it(`${name}: every position is the API's window value`, async () => {
      const p = await provider().fetch(r, signal);
      const value = p.kind === 'items' && (p.items.value as Float64Array);
      if (!value) throw new Error('not items');
      expect(value.length).toBe(r.stop - r.start);
      const want = new Float64Array(r.stop - r.start).fill(Number.NaN);
      for (const [s, e, v] of expected[name] as [number, number, number][]) want.fill(v, s - r.start, e - r.start);
      expect(value).toEqual(want);
      expect(value.every((v) => v >= 0 && v <= 100)).toBe(true); // percent GC, no gaps in these regions
    });

  for (const [name, r] of CASES.bigwig.slice(1).filter(([, r]) => r.level >= 0))
    it(`${name}: bins hold the API's windows, resampled from the zoom records`, async () => {
      const p = await provider().fetch(r, signal);
      if (p.kind !== 'bins') throw new Error('not bins');
      const { min, max, mean } = p.bins as Record<string, Float64Array>;
      const want = expected[name] as ([number, number, number] | null)[];
      expect(min.length).toBe(want.length);
      let [sum, wsum, n] = [0, 0, 0];
      want.forEach((w, j) => {
        if (!w) {
          expect(Number.isNaN(mean[j])).toBe(true);
          return;
        }
        // A zoom record is not on the level's grid, so a bin takes the extremes of the records it touches: never inside the data's.
        expect(min[j]).toBeLessThanOrEqual(w[0]);
        expect(max[j]).toBeGreaterThanOrEqual(w[1]);
        [sum, wsum, n] = [sum + mean[j], wsum + w[2], n + 1];
      });
      expect(n).toBeGreaterThan(0);
      expect(Math.abs(sum / n - wsum / n)).toBeLessThan(1); // the means keep the total
    });

  it('fails clearly for what it does not serve', async () => {
    const p = provider();
    await expect(p.fetch(req('chr21', 0, 256, 10), signal)).rejects.toBeInstanceOf(StreamNotServed);
    await expect(p.fetch(req('chrZ', 0, 256), signal)).rejects.toThrow('BigWig: no chromosome chrZ');
  });

  it('fails naming the range when a region is not in the recorded slices', async () => {
    await expect(provider().fetch(req('chr2', 50_000_000, 50_004_096), signal)).rejects.toThrow(
      /bigwig: bytes \d+\.\.\d+ are not in the recorded slices/,
    );
  });
});
