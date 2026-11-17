/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seeded properties at small, fixed counts: the range-set operations against a bitset oracle, identity however a
// Location is built, and a short fuzz of normalization and overlap.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import {
  hashLocation,
  intersectLocations,
  location,
  locationsEqual,
  normalizeRanges,
  overlapsRanges,
  partOf,
  subtractLocations,
  unionLocations,
} from '../location.js';

const Size = 48;
const Cases = 60;

/** Random `[start, end, …]` pairs on [0, Size], empty and reversed ones among them. */
function randomPairs(rng: PCG, max = 6): number[] {
  const n = rng.int() % (max + 1);
  const out: number[] = [];
  for (let i = 0; i < n; i++) out.push(rng.int() % (Size + 1), rng.int() % (Size + 1));
  return out;
}

/** The oracle: one bit per position, a pair covering `[min(s, e), max(s, e))` only when `s < e` (as ranges do). */
function bits(pairs: ArrayLike<number>): boolean[] {
  const b = new Array<boolean>(Size).fill(false);
  for (let i = 0; i < pairs.length; i += 2) for (let p = pairs[i]; p < pairs[i + 1]; p++) b[p] = true;
  return b;
}

const rangesOf = (l: ReturnType<typeof location>, system = '/s') => partOf(l, system)?.ranges ?? new Int32Array(0);

describe('query: Location properties', () => {
  it('normalizeRanges gives sorted, disjoint, non-adjacent, non-empty ranges covering the same positions', () => {
    const rng = new PCG(7);
    for (let c = 0; c < Cases; c++) {
      const pairs = randomPairs(rng, 10);
      const r = normalizeRanges(pairs);
      for (let i = 0; i < r.length; i += 2) {
        expect(r[i]).toBeLessThan(r[i + 1]);
        if (i) expect(r[i]).toBeGreaterThan(r[i - 1]);
      }
      expect(bits(r)).toEqual(bits(pairs));
    }
  });

  it('union, intersect and subtract match the bitset oracle', () => {
    const rng = new PCG(11);
    for (let c = 0; c < Cases; c++) {
      const pa = randomPairs(rng);
      const pb = randomPairs(rng);
      const a = location({ parts: [{ system: '/s', ranges: pa }] });
      const b = location({ parts: [{ system: '/s', ranges: pb }] });
      const x = bits(pa);
      const y = bits(pb);
      expect(bits(rangesOf(unionLocations(a, b)))).toEqual(x.map((v, i) => v || y[i]));
      expect(bits(rangesOf(intersectLocations(a, b)))).toEqual(x.map((v, i) => v && y[i]));
      expect(bits(rangesOf(subtractLocations(a, b)))).toEqual(x.map((v, i) => v && !y[i]));
      // and the results are normalized: equal to the normalization of themselves
      const u = unionLocations(a, b);
      expect(locationsEqual(u, location({ parts: u.parts }))).toBe(true);
    }
  });

  it('equal Locations hash equal however built: pairs shuffled, split, given per system in any order', () => {
    const rng = new PCG(13);
    for (let c = 0; c < Cases; c++) {
      const ps = randomPairs(rng);
      const pt = randomPairs(rng);
      const a = location({
        parts: [
          { system: '/s', ranges: ps },
          { system: '/t', ranges: pt },
        ],
      });
      // the same positions: one pair per position, reversed, systems swapped, split in two parts
      const singles = (pairs: number[]) => {
        const out: number[] = [];
        bits(pairs).forEach((v, p) => v && out.unshift(p, p + 1));
        return out;
      };
      const st = singles(pt);
      const b = location({
        parts: [
          { system: '/t', ranges: st.slice(0, (st.length >> 2) << 1) },
          { system: '/s', ranges: singles(ps) },
          { system: '/t', ranges: st.slice((st.length >> 2) << 1) },
        ],
      });
      expect(locationsEqual(a, b)).toBe(true);
      expect(hashLocation(a)).toBe(hashLocation(b));
    }
  });

  it('overlapsRanges matches a scan for items and boundaries, at the system end too', () => {
    const rng = new PCG(17);
    for (let c = 0; c < Cases; c++) {
      const r = normalizeRanges(randomPairs(rng));
      const b = bits(r);
      for (let k = 0; k < 8; k++) {
        const s = rng.int() % (Size + 1);
        const e = Math.min(Size, s + (rng.int() % 6));
        const expected =
          s === e ? (s === Size ? r.length > 0 && r[r.length - 1] === Size : b[s]) : b.slice(s, e).some((v) => v);
        expect(overlapsRanges(r, s, e, Size), `${[...r]} [${s}, ${e})`).toBe(expected);
      }
    }
  });

  it('refuses a position outside int32, never wraps it', () => {
    expect(() => normalizeRanges([0, 2 ** 31])).toThrow(/int32/);
    expect(() => normalizeRanges([-(2 ** 31) - 1, 0])).toThrow(/int32/);
    expect([...normalizeRanges([2 ** 31 - 2, 2 ** 31 - 1])]).toEqual([2 ** 31 - 2, 2 ** 31 - 1]);
  });
});
