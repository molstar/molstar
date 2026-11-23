/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { formatGenomic } from '../number.js';
import { _parsedFormatSpecCount, formatValue } from '../string-format.js';

describe('formatGenomic', () => {
  it('bases below a step of 1,000, grouped by thousands', () => {
    expect(formatGenomic(7_668_402, 1)).toBe('7,668,402');
    expect(formatGenomic(7_668_500, 500)).toBe('7,668,500');
    expect(formatGenomic(400, 100)).toBe('400');
    expect(formatGenomic(0, 10)).toBe('0');
  });

  it('kilobases below a step of 10,000', () => {
    expect(formatGenomic(7_668_000, 1000)).toBe('7,668 kb');
    expect(formatGenomic(5000, 5000)).toBe('5 kb');
  });

  it('megabases from a step of 10,000, with the decimals the step needs', () => {
    expect(formatGenomic(7_670_000, 10_000)).toBe('7.67 Mb');
    expect(formatGenomic(7_700_000, 20_000)).toBe('7.70 Mb');
    expect(formatGenomic(7_600_000, 100_000)).toBe('7.6 Mb');
    expect(formatGenomic(7_500_000, 500_000)).toBe('7.5 Mb');
    expect(formatGenomic(8_000_000, 1_000_000)).toBe('8 Mb');
    expect(formatGenomic(80_000_000, 10_000_000)).toBe('80 Mb');
    expect(formatGenomic(2_000_000_000, 1_000_000_000)).toBe('2,000 Mb');
  });

  it('a negative number keeps its sign', () => {
    expect(formatGenomic(-1500, 1)).toBe('-1,500');
    expect(formatGenomic(-2_000_000, 1_000_000)).toBe('-2 Mb');
    // a value that rounds to zero is unsigned
    expect(formatGenomic(-0.3, 1)).toBe('0');
    expect(formatGenomic(-400, 1000)).toBe('0 kb');
    expect(formatGenomic(-1000, 10_000)).toBe('0.00 Mb');
  });
});

describe('formatGenomic parses each spec once', () => {
  const group = (s: string) => s.replace(/^(-?)(\d+)/, (_, sign, d) => sign + d.replace(/\B(?=(\d{3})+(?!\d))/g, ','));

  it('formatValue `,.Nf` (N 0 to 6) is byte-equal to toFixed and grouping, for random numbers', () => {
    let seed = 12345;
    const rnd = () => (seed = (seed * 1664525 + 1013904223) % 4294967296) / 4294967296;
    for (let i = 0; i < 3000; i++) {
      const v = (rnd() - 0.4) * 10 ** Math.floor(rnd() * 12);
      for (let d = 0; d <= 6; d++) {
        expect(formatValue(String(v), `,.${d}f`)).toBe(group(v.toFixed(d)));
      }
    }
  });

  it('formatGenomic equals the unit and decimals it documents, for random numbers and steps', () => {
    let seed = 777;
    const rnd = () => (seed = (seed * 1664525 + 1013904223) % 4294967296) / 4294967296;
    for (let i = 0; i < 3000; i++) {
      const n = Math.round((rnd() - 0.3) * 10 ** Math.floor(rnd() * 11));
      const step = 10 ** Math.floor(rnd() * 10) * [1, 2, 5][Math.floor(rnd() * 3)];
      const unit = step < 1e3 ? 1 : step < 1e4 ? 1e3 : 1e6;
      const decimals = unit === 1e6 ? Math.max(0, 6 - Math.floor(Math.log10(step) + 1e-9)) : 0;
      const v = n / unit;
      const shown = Number(v.toFixed(decimals)) === 0 ? 0 : v;
      const expected = group(shown.toFixed(decimals)) + (unit === 1 ? '' : unit === 1e3 ? ' kb' : ' Mb');
      expect(formatGenomic(n, step)).toBe(expected);
    }
  });

  it('parses no spec per call: the parsed-spec count stays at one a distinct spec', () => {
    for (let d = 0; d <= 6; d++) formatValue('1', `,.${d}f`);
    const count = _parsedFormatSpecCount();
    for (let i = 0; i < 1000; i++) {
      for (let d = 0; d <= 6; d++) formatValue('1234.5', `,.${d}f`);
      formatGenomic(i * 1000, 10_000);
    }
    expect(_parsedFormatSpecCount()).toBe(count);
  });
});
