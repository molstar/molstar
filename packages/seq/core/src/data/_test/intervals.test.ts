/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Intervals found by a region (_spec/README.md "Data", Intervals): the index's matches, by a callback and into a buffer
// (the count past it), equal a brute-force overlap over random intervals, long ones among them; a zero-length interval
// meets the position it is at.

import { itWithin, ratioWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { eachMeeting, intervalIndex, meetingInto } from '../intervals.js';

const meets = (s: number, e: number, a: number, b: number) => s < b && Math.max(e, s + 1) > a;

describe('the interval index', () => {
  it('finds what a brute-force overlap finds, by callback and into a buffer with the count past it', () => {
    let seed = 11;
    const rand = () => {
      seed = (seed * 16807) % 2147483647;
      return seed / 2147483647;
    };
    const n = 2000;
    const s = Int32Array.from({ length: n }, () => Math.floor(rand() * 1e6));
    // most short, a few spanning much of the system (the index's long list)
    const e = Int32Array.from(s, (x, i) => x + (i % 97 === 0 ? Math.floor(rand() * 5e5) : Math.floor(rand() * 200)));
    const ix = intervalIndex(s, e);
    expect(ix.long.length).toBeGreaterThan(0);
    const buffer = new Int32Array(4);
    for (let q = 0; q < 200; q++) {
      const a = Math.floor(rand() * 1e6);
      const b = a + Math.floor(rand() * 3000);
      const brute: number[] = [];
      for (let i = 0; i < n; i++) if (meets(s[i], e[i], a, b)) brute.push(i);
      const found: number[] = [];
      eachMeeting(ix, s, e, a, b, (i) => found.push(i));
      expect(found.sort((x, y) => x - y)).toEqual(brute);
      const count = meetingInto(ix, s, e, a, b, buffer);
      expect(count).toBe(brute.length);
      for (let k = 0; k < Math.min(count, buffer.length); k++) expect(brute).toContain(buffer[k]);
    }
  });

  // a tooltip's lookup under a hover: O(log n + k), so a million-range table answers as a ten-thousand one does
  itWithin('the lookup over 1,000,000 ranges within 1.5× of 10,000 (a hover’s cap of matches)', (within) => {
    const table = (n: number) => {
      const s = Int32Array.from({ length: n }, (_, i) => i * 4);
      const e = Int32Array.from(s, (x) => x + 3);
      return { s, e, ix: intervalIndex(s, e) };
    };
    const small = table(10_000);
    const large = table(1_000_000);
    const out = new Int32Array(9);
    const queries = (t: typeof small, n: number) => () => {
      let hits = 0;
      for (let q = 0; q < 10_000; q++) {
        // a hover's path: positions next to each other, from the middle of the table on
        const at = (n >> 1) * 4 + q;
        hits += meetingInto(t.ix, t.s, t.e, at, at + 1, out);
      }
      return hits;
    };
    ratioWithin(within, 1.5, queries(small, 10_000), queries(large, 1_000_000));
  });

  it('a zero-length interval meets the position it is at', () => {
    const ix = intervalIndex([5], [5]);
    expect(meetingInto(ix, [5], [5], 5, 6, new Int32Array(1))).toBe(1);
    expect(meetingInto(ix, [5], [5], 6, 7, new Int32Array(1))).toBe(0);
  });
});
