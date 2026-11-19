/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { formatGenomic } from '../number.js';

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
