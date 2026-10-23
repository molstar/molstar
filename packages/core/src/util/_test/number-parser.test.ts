/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { parseFloat, parseFloatBytes, parseInt, parseIntBytes } from '../number-parser.js';

/** mulberry32: a seeded PRNG, so a failure reproduces */
function rng(seed: number) {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Number-like tokens: signs, points, exponents, long digit runs, garbage and an empty value. */
function token(r: () => number) {
  const pick = (xs: string) => xs[Math.floor(r() * xs.length)];
  const digits = (n: number) => Array.from({ length: n }, () => pick('0123456789')).join('');
  const parts = [
    r() < 0.3 ? pick('-+') : '',
    digits(Math.floor(r() * (r() < 0.1 ? 25 : 6))),
    r() < 0.6 ? `.${digits(Math.floor(r() * 8))}` : '',
    r() < 0.2 ? `${pick('eE')}${r() < 0.5 ? pick('-+') : ''}${digits(Math.floor(r() * 3))}` : '',
    r() < 0.1 ? pick('x?.é ') : '',
  ];
  return parts.join('');
}

describe('number-parser bytes', () => {
  it('parses bytes bit-identically to the string parsers', () => {
    const r = rng(7);
    const tokens = Array.from({ length: 20000 }, () => token(r));
    tokens.push('', '-', '+', '.', '-.', 'e5', '1e', '1e+', '12abc', '2147483648', '-2147483649', '1.5E-3', '?', '.');
    // as a CIF parser sees them: ranges in one text, a space after each
    const text = tokens.join(' ');
    const bytes = new TextEncoder().encode(text);
    let charStart = 0;
    let byteStart = 0;
    for (const t of tokens) {
      const byteLength = new TextEncoder().encode(t).length;
      const [cs, ce, bs, be] = [charStart, charStart + t.length, byteStart, byteStart + byteLength];
      expect(Object.is(parseFloatBytes(bytes, bs, be), parseFloat(text, cs, ce)), t).toBe(true);
      expect(Object.is(parseIntBytes(bytes, bs, be), parseInt(text, cs, ce)), t).toBe(true);
      charStart = ce + 1;
      byteStart = be + 1;
    }
  });
});
