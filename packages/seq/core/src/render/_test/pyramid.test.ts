/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { ImageCache, ImageTileBytes } from '../images.js';
import { ExactMajorityRows, NoClass, Pyramid, pyramidLevelFor } from '../pyramid.js';

/** A random grid of classes, each column biased to one class so majorities are mostly clear. */
function classes(rng: PCG, rows: number, columns: number, k = 6): Uint8Array {
  const out = new Uint8Array(rows * columns);
  const bias = Array.from({ length: columns }, () => rng.int() % k);
  for (let i = 0; i < rows; i++)
    for (let j = 0; j < columns; j++) out[i * columns + j] = rng.float() < 0.6 ? bias[j] : rng.int() % k;
  return out;
}

/**
 * The exact majority and its share over a block of the base (classes below 8, as `classes` makes them); `shareOf` reads
 * the block's counts until the next call.
 */
const counts = new Float64Array(8);
function exact(base: Uint8Array, columns: number, r0: number, r1: number, c0: number, c1: number) {
  counts.fill(0);
  let n = 0;
  for (let r = r0; r < r1; r++)
    for (let c = c0; c < c1; c++) {
      counts[base[r * columns + c]]++;
      n++;
    }
  let letter = 0;
  for (let v = 1; v < 8; v++) if (counts[v] > counts[letter]) letter = v;
  let second = 0;
  for (let v = 0; v < 8; v++) if (v !== letter && counts[v] > second) second = counts[v];
  return {
    letter,
    share: counts[letter] / n,
    lead: (counts[letter] - second) / n,
    shareOf: (v: number) => (v < 8 ? counts[v] / n : 0),
  };
}

describe('the shared reduction (pyramid)', () => {
  it('majority levels pair independently: zoomed out on columns a row stays a row', () => {
    const p = new Pyramid({ factor: 4, pairing: 'independent', reducer: 'majority' });
    p.setClasses(10, 40, classes(new PCG(1), 10, 40));
    const l = p.level(0, 1);
    expect([l.rows, l.columns]).toEqual([10, 10]);
    const m = p.level(2, 0);
    expect([m.rows, m.columns]).toEqual([1, 40]);
    expect(p.levels('rows')).toBe(3);
    // clamped to the levels there are
    expect(p.level(9, 9).rows).toBe(1);
    expect(p.level(9, 9).columns).toBe(1);
  });

  it('a block under 64 rows is its exact majority; a larger one the majority of majorities within its tolerance: the letter where the winner leads by more than 0.25, the share within 0.05 on average and 0.3 at worst (seeds 1 to 12)', () => {
    // a majority is not mergeable: a majority of majorities misses the letter in small blocks (measured on PCG32,
    // seeds 1 to 12: leads up to 0.41), so a block under `ExactMajorityRows` rows counts its cells; the recorded
    // tolerance (_spec/README.md "Render", the rows layer) holds for the larger ones
    expect(ExactMajorityRows).toBe(64);
    let sum = 0;
    let count = 0;
    let worst = 0;
    let inexact = 0;
    let misses = 0;
    for (let seed = 1; seed <= 12; seed++) {
      const rng = new PCG(seed);
      for (let t = 0; t < 4; t++) {
        const rows = 70 + (rng.int() % 400);
        const columns = 8 + (rng.int() % 40);
        const base = classes(rng, rows, columns);
        const p = new Pyramid({ factor: 4, pairing: 'independent', reducer: 'majority' });
        p.setClasses(rows, columns, base);
        for (const [rl, cl] of [
          [1, 0],
          [2, 0],
          [0, 1],
          [1, 1],
          [2, 1],
          [3, 0],
          [3, 1],
          [4, 0],
          [4, 1],
        ]) {
          const l = p.level(rl, cl);
          const fr = 4 ** l.rowLevel;
          const fc = 4 ** l.columnLevel;
          const exactLevel = fr < ExactMajorityRows;
          for (let i = 0; i < l.rows; i++)
            for (let j = 0; j < l.columns; j++) {
              const e = exact(
                base,
                columns,
                i * fr,
                Math.min(rows, (i + 1) * fr),
                j * fc,
                Math.min(columns, (j + 1) * fc),
              );
              const o = i * l.columns + j;
              // the share against its letter's exact share, a wrong letter's too
              const d = Math.abs(l.shares[o] / 255 - e.shareOf(l.classes[o]));
              if (exactLevel) {
                if (l.classes[o] !== e.letter || d > 0.5 / 255 + 1e-9) inexact++;
                continue;
              }
              if (e.lead > 0.25 && l.classes[o] !== e.letter) misses++;
              worst = Math.max(worst, d);
              sum += d;
              count++;
            }
        }
      }
    }
    expect(count).toBeGreaterThan(500);
    expect(inexact).toBe(0);
    expect(misses).toBe(0);
    expect(sum / count).toBeLessThanOrEqual(0.05);
    expect(worst).toBeLessThanOrEqual(0.3);
  });

  it('square levels of numbers reduce by mean, max and min, NaN not counted', () => {
    const v = new Float32Array([1, 2, Number.NaN, 4, 5, 6, 7, 8, 9]);
    for (const [reducer, want] of [
      ['mean', [3, 6, 7.5, 9]],
      ['max', [5, 6, 8, 9]],
      ['min', [1, 6, 7, 9]],
    ] as const) {
      const p = new Pyramid({ factor: 2, pairing: 'square', reducer });
      p.setValues(3, 3, v);
      const l = p.level(1);
      expect([l.rows, l.columns]).toEqual([2, 2]);
      expect(Array.from(l.values.subarray(0, 4))).toEqual(want);
      // a mean of means (each level from the one below), the extremes exact
      expect(p.level(5).values[0]).toBe(reducer === 'mean' ? 6.375 : reducer === 'max' ? 9 : 1);
    }
    // the one-axis case: a row of bars, its max levels
    const bars = new Pyramid({ factor: 2, pairing: 'independent', reducer: 'max' });
    bars.setValues(1, 5, new Float32Array([1, 3, 2, 0, 4]));
    expect(Array.from(bars.level(0, 1).values.subarray(0, 3))).toEqual([3, 2, 4]);
  });

  it('reused buffers give what fresh ones give (a base after a larger one)', () => {
    const rng = new PCG(3);
    const big = classes(rng, 300, 120);
    const small = classes(rng, 70, 33);
    const reused = new Pyramid({ factor: 4, pairing: 'independent', reducer: 'majority' });
    reused.setClasses(300, 120, big);
    reused.level(2, 2);
    reused.level(1, 0);
    reused.setClasses(70, 33, small);
    const fresh = new Pyramid({ factor: 4, pairing: 'independent', reducer: 'majority' });
    fresh.setClasses(70, 33, small);
    for (const [r, c] of [
      [1, 0],
      [2, 2],
      [0, 2],
    ]) {
      const a = reused.level(r, c);
      const b = fresh.level(r, c);
      const n = a.rows * a.columns;
      expect([a.rows, a.columns]).toEqual([b.rows, b.columns]);
      expect(a.classes.subarray(0, n)).toEqual(b.classes.subarray(0, n));
      expect(a.shares.subarray(0, n)).toEqual(b.shares.subarray(0, n));
    }
  });

  it('an empty sample has no class; the level for a pixel is the first whose sample covers it', () => {
    const p = new Pyramid({ factor: 4, pairing: 'independent', reducer: 'majority' });
    p.setClasses(2, 2, new Uint8Array([NoClass, NoClass, NoClass, 1]));
    expect(p.level(1, 1).classes[0]).toBe(1);
    expect(pyramidLevelFor(0.5, 4)).toBe(0);
    expect(pyramidLevelFor(1, 4)).toBe(0);
    expect(pyramidLevelFor(4, 4)).toBe(1);
    expect(pyramidLevelFor(5, 4)).toBe(2);
  });
});

describe('the shared image cache', () => {
  it('holds its bytes over a tile’s, lends an evicted tile’s bytes to the next, and reads by key', () => {
    const c = new ImageCache(3 * ImageTileBytes);
    expect(c.capacity).toBe(3);
    for (let i = 0; i < 3; i++) c.put(`t${i}`, c.slot(256, 256));
    expect(c.made).toBe(3);
    const t0 = c.get('t0')!;
    // t1 is now the least recently used: evicted, its bytes reused
    c.put('t3', c.slot(256, 100));
    expect(c.get('t1')).toBeUndefined();
    const t4 = c.slot(256, 256);
    expect(c.made).toBe(4);
    c.put('t4', t4);
    // t2 evicted: its bytes go to the next slot
    const t5 = c.slot(10, 10);
    expect(c.made).toBe(4);
    expect(t5.data.length).toBe(400);
    expect(c.get('t0')).toBe(t0);
    c.clear();
    expect(c.size).toBe(0);
  });
});
