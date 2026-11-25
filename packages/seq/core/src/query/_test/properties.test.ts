/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seeded properties at small, fixed counts: the range-set operations against a bitset oracle, identity however a
// Location is built; the ranges' own properties are tested in core (`HalfOpenRanges`).

import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import {
  hashLocation,
  intersectLocations,
  location,
  locationsEqual,
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
  it('a rect minus a rect covers exactly the difference of their cells, in at most two rects', () => {
    const rng = new PCG(29);
    const S = 16;
    const cells = (l: ReturnType<typeof location>) => {
      const m = new Array<boolean>(S * S).fill(false);
      for (const r of l.rects)
        for (let i = 0; i < r.a.ranges.length; i += 2)
          for (let x = r.a.ranges[i]; x < r.a.ranges[i + 1]; x++)
            for (let j = 0; j < r.b.ranges.length; j += 2)
              for (let y = r.b.ranges[j]; y < r.b.ranges[j + 1]; y++) m[x * S + y] = true;
      return m;
    };
    const side = () => {
      const s = rng.int() % S;
      return [s, s + 1 + (rng.int() % (S - s))];
    };
    const rect = (a: number[], b: number[]) => ({ a: { system: '/x', ranges: a }, b: { system: '/y', ranges: b } });
    for (let c = 0; c < Cases; c++) {
      const r = location({ rects: [rect(side(), side())] });
      const q = location({ rects: [rect(side(), side())] });
      const d = subtractLocations(r, q);
      const x = cells(r);
      const y = cells(q);
      expect(cells(d)).toEqual(x.map((v, i) => v && !y[i]));
      expect(d.rects.length).toBeLessThanOrEqual(2);
      // a cell is in one rect of the difference at most: no overlap to merge
      let n = 0;
      for (const rr of d.rects) n += cells(location({ rects: [rr] })).filter(Boolean).length;
      expect(n).toBe(cells(d).filter(Boolean).length);
    }
    // over other systems, untouched
    const r = location({ rects: [rect([0, 4], [0, 4])] });
    const other = location({ rects: [{ a: { system: '/x', ranges: [0, 4] }, b: { system: '/z', ranges: [0, 4] } }] });
    expect(locationsEqual(subtractLocations(r, other), r)).toBe(true);
    expect(subtractLocations(r, r).rects).toHaveLength(0);
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
});
