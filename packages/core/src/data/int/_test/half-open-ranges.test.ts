/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seeded properties at small, fixed counts: normalization, set operations (against a bitset oracle) and overlap.

import { describe, expect, it } from 'vitest';
import { PCG } from '../../util/hash-functions.js';
import { HalfOpenRanges } from '../half-open-ranges.js';

const Size = 48;
const Cases = 60;

/** Random `[start, end, …]` pairs on [0, Size], empty and reversed ones among them. */
function randomPairs(rng: PCG, max = 6): number[] {
  const n = rng.int() % (max + 1);
  const out: number[] = [];
  for (let i = 0; i < n; i++) out.push(rng.int() % (Size + 1), rng.int() % (Size + 1));
  return out;
}

/** The oracle: one bit per position, a pair covering `[s, e)` only when `s < e` (as ranges do). */
function bits(pairs: ArrayLike<number>): boolean[] {
  const b = new Array<boolean>(Size).fill(false);
  for (let i = 0; i < pairs.length; i += 2) for (let p = pairs[i]; p < pairs[i + 1]; p++) b[p] = true;
  return b;
}

describe('HalfOpenRanges', () => {
  it('normalize gives sorted, disjoint, non-adjacent, non-empty ranges covering the same positions', () => {
    const rng = new PCG(7);
    for (let c = 0; c < Cases; c++) {
      const pairs = randomPairs(rng, 10);
      const r = HalfOpenRanges.normalize(pairs);
      for (let i = 0; i < r.length; i += 2) {
        expect(r[i]).toBeLessThan(r[i + 1]);
        if (i) expect(r[i]).toBeGreaterThan(r[i - 1]);
      }
      expect(bits(r)).toEqual(bits(pairs));
    }
  });

  it('union, intersect, subtract and size match the bitset oracle', () => {
    const rng = new PCG(11);
    for (let c = 0; c < Cases; c++) {
      const pa = randomPairs(rng);
      const pb = randomPairs(rng);
      const a = HalfOpenRanges.normalize(pa);
      const b = HalfOpenRanges.normalize(pb);
      const x = bits(pa);
      const y = bits(pb);
      const results = [
        [HalfOpenRanges.union(a, b), x.map((v, i) => v || y[i])],
        [HalfOpenRanges.intersect(a, b), x.map((v, i) => v && y[i])],
        [HalfOpenRanges.subtract(a, b), x.map((v, i) => v && !y[i])],
      ] as const;
      for (const [r, expected] of results) {
        expect(bits(r)).toEqual(expected);
        // normalized: equal to the normalization of itself
        expect(HalfOpenRanges.compare(r, HalfOpenRanges.normalize(r))).toBe(0);
        expect(HalfOpenRanges.size(r)).toBe(expected.filter(Boolean).length);
      }
    }
  });

  it('overlaps matches a scan for items and boundaries, at the end of the domain too', () => {
    const rng = new PCG(17);
    for (let c = 0; c < Cases; c++) {
      const r = HalfOpenRanges.normalize(randomPairs(rng));
      const b = bits(r);
      for (let k = 0; k < 8; k++) {
        const s = rng.int() % (Size + 1);
        const e = Math.min(Size, s + (rng.int() % 6));
        const expected =
          s === e ? (s === Size ? r.length > 0 && r[r.length - 1] === Size : b[s]) : b.slice(s, e).some((v) => v);
        expect(HalfOpenRanges.overlaps(r, s, e, Size), `${[...r]} [${s}, ${e})`).toBe(expected);
      }
    }
  });

  it('compares pairwise, then by length', () => {
    const r = (...v: number[]) => HalfOpenRanges.normalize(v);
    expect(HalfOpenRanges.compare(r(0, 2), r(0, 2))).toBe(0);
    expect(HalfOpenRanges.compare(r(0, 2), r(1, 2))).toBeLessThan(0);
    expect(HalfOpenRanges.compare(r(0, 2), r(0, 2, 4, 5))).toBeLessThan(0);
  });

  it('refuses a value outside int32, never wraps it', () => {
    expect(() => HalfOpenRanges.normalize([0, 2 ** 31])).toThrow(/int32/);
    expect(() => HalfOpenRanges.normalize([-(2 ** 31) - 1, 0])).toThrow(/int32/);
    expect(() => HalfOpenRanges.normalize([1])).toThrow();
    expect([...HalfOpenRanges.normalize([2 ** 31 - 2, 2 ** 31 - 1])]).toEqual([2 ** 31 - 2, 2 ** 31 - 1]);
  });
});
