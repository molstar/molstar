/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// KS1: row statistics and greedy redundancy clustering (seqspec:format#computed-nodes, seqspec:model#providers).

import { describe, expect, it } from 'vitest';
import {
  identityOf,
  Redundancy,
  RedundancySizeError,
  type RowCodes,
  redundancy,
  rowStats,
} from '../../src/algorithms/index.js';
import { itWithin } from '../model/timing.js';

const GAP = 45;
/** Rows of text as cells; `-` is the gap. */
function cells(rows: string[]): RowCodes {
  const columns = Math.max(...rows.map((r) => r.length));
  const codes = new Uint8Array(rows.length * columns).fill(GAP);
  rows.forEach((r, i) => {
    for (let c = 0; c < r.length; c++) codes[i * columns + c] = r.charCodeAt(c);
  });
  return { codes, rows: rows.length, columns, gap: GAP };
}

describe('rowStats', () => {
  const c = cells(['ACDEFGHIKL', 'ACDEFGHIKM', 'AC--FGHI--', '----------', 'WYWYWYWYWY']);

  it('identity is equal letters over the columns where both rows have one, NaN with none', () => {
    const { identity } = rowStats(c, 0);
    expect(Array.from(identity.slice(0, 3))).toEqual([1, 0.9, 1]);
    expect(identity[3]).toBeNaN();
    expect(identity[4]).toBe(0);
    expect(identityOf(c, 1, 2)).toBe(1);
    expect(identityOf(c, 4, 3)).toBeNaN();
  });

  it('coverage is the share of columns with a letter, length their count', () => {
    const { coverage, length } = rowStats(c, 0);
    expect(Array.from(coverage)).toEqual([1, 1, 0.6, 0, 1]);
    expect(Array.from(length)).toEqual([10, 10, 6, 0, 10]);
  });

  it('is to the reference row asked for, and any alphabet: symbols are codes', () => {
    const { identity } = rowStats(c, 4);
    expect(Array.from(identity)).toEqual([0, 0, 0, Number.NaN, 1].map((x) => x));
    const dna = cells(['ACGT', 'ACGA']);
    expect(rowStats(dna, 0).identity[1]).toBe(0.75);
  });
});

describe('redundancy', () => {
  const c = cells(['ACDEFGHIKL', 'ACDEFGHIKM', 'ACDEFGHIAM', 'ACDEFGHIKV', 'WYWYWYWYWY', 'WYWYWY----']);

  it('a row joins the first representative it is at least threshold identical to, else is one', async () => {
    const r = await redundancy(c, 0.8, null);
    expect(Array.from(r.cluster)).toEqual([0, 0, 0, 0, 4, 4]);
    expect(Array.from(r.representative)).toEqual([1, 0, 0, 0, 1, 0]);
    const hi = await redundancy(c, 0.9, null);
    expect(Array.from(hi.cluster)).toEqual([0, 0, 2, 0, 4, 4]);
  });

  it('threshold 1 clusters only identical rows; 0 puts every row with a shared column in the first', async () => {
    const same = cells(['AAAA', 'AAAA', 'AAAB']);
    expect(Array.from((await redundancy(same, 1, null)).cluster)).toEqual([0, 0, 2]);
    expect(Array.from((await redundancy(same, 0, null)).cluster)).toEqual([0, 0, 0]);
  });

  it('rows with no common column never join, whatever the threshold', async () => {
    const r = await redundancy(cells(['AA--', '--BB']), 0, null);
    expect(Array.from(r.cluster)).toEqual([0, 1]);
  });

  it('takes rows in the given order, the ones it leaves out after, in row order', async () => {
    const r = await redundancy(c, 0.9, [3, 2]);
    // 3 and 2 first (2 is 0.8 to 3: separate), then 0 (0.9 to 3), 1 (0.9 to 2 and 3: the first, 3), 4, 5.
    expect(Array.from(r.cluster)).toEqual([3, 3, 2, 3, 4, 4]);
    expect(Array.from(r.representative)).toEqual([0, 0, 1, 1, 1, 0]);
  });

  it('is bounded: past Redundancy.MaxCells it throws before any work', async () => {
    const big: RowCodes = { codes: new Uint8Array(0), rows: 1300, columns: 1300, gap: GAP };
    expect(1300 ** 3).toBeGreaterThan(Redundancy.MaxCells);
    await expect(redundancy(big, 0.9, null)).rejects.toBeInstanceOf(RedundancySizeError);
  });

  itWithin('honors an abort signal, before and during a long run', async (within) => {
    const ac = new AbortController();
    ac.abort();
    await expect(redundancy(c, 0.9, null, ac.signal)).rejects.toThrow();
    // 600 distinct rows of 1,200 columns: every row compares with every representative, about 2e8 steps.
    const rows = Array.from(
      { length: 600 },
      (_, r) => String.fromCharCode(65 + (r % 26)).repeat(600) + String(r).padStart(600, 'A'),
    );
    const long = cells(rows);
    const during = new AbortController();
    setTimeout(() => during.abort(), 20);
    await within(5000, () => expect(redundancy(long, 1, null, during.signal)).rejects.toThrow());
  });
});
