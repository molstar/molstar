/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ElementSymbol } from '../types.js';

// The symbols are memoized lazily (perf audit 2, 12/P1): the same strings out as the table built at load gave.
describe('ElementSymbol', () => {
  it('upper-cases what it reads, one string per symbol', () => {
    for (const [s, e] of [
      ['C', 'C'],
      ['fe', 'FE'],
      ['Zn', 'ZN'],
      ['cl', 'CL'],
      ['x9', 'X9'],
      ['Uuo', 'UUO'],
      ['', ''],
      [' n', ' N'],
      ['abcd', 'ABCD'],
    ]) {
      expect(ElementSymbol(s)).toBe(e);
      expect(ElementSymbol(s)).toBe(ElementSymbol(s));
    }
  });

  it('takes a number as its digits, as a MolScript literal may give it', () => {
    expect(ElementSymbol(6 as unknown as string)).toBe('6');
  });
});
