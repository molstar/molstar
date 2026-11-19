/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Determine the number of digits in a floating point number
 * Find a number M such that round(M * v) - M * v <= delta.
 * If no such M exists, return -1.
 */
export function getMantissaMultiplier(v: number, maxDigits: number, delta: number) {
  let m = 1,
    i;
  for (i = 0; i < maxDigits; i++) {
    const mv = m * v;
    if (Math.abs(Math.round(mv) - mv) <= delta) return i;
    m *= 10;
  }
  return -1;
}

export function integerDigitCount(v: number, delta: number) {
  const f = Math.abs(v);
  if (f < delta) return 0;
  return Math.floor(Math.log10(Math.abs(v))) + 1;
}

/**
 * Determine the maximum number of digits in a floating point array.
 * Find a number M such that round(M * v) - M * v <= delta.
 * If no such M exists, return -1.
 */
export function getArrayDigitCount(xs: ArrayLike<number>, maxDigits: number, delta: number) {
  let mantissaDigits = 1;
  let integerDigits = 0;
  for (let i = 0, _i = xs.length; i < _i; i++) {
    if (mantissaDigits >= 0) {
      const t = getMantissaMultiplier(xs[i], maxDigits, delta);
      if (t < 0) mantissaDigits = -1;
      else if (t > mantissaDigits) mantissaDigits = t;
    }
    const abs = Math.abs(xs[i]);
    if (abs > delta) {
      const d = Math.floor(Math.log10(Math.abs(abs))) + 1;
      if (d > integerDigits) integerDigits = d;
    }
  }
  return { mantissaDigits, integerDigits };
}

export function isInteger(s: string) {
  s = s.trim();
  const n = parseInt(s, 10);
  return Number.isNaN(n) ? false : n.toString() === s;
}

export function getPrecision(v: number) {
  if (!isFinite(v)) return 0;
  let e = 1;
  let p = 0;
  while (Math.round(v * e) / e !== v) {
    e *= 10;
    ++p;
  }
  return p;
}

export function toPrecision(v: number, precision: number) {
  return parseFloat(v.toPrecision(precision));
}

export function toFixed(v: number, fractionDigits: number) {
  return parseFloat(v.toFixed(fractionDigits));
}

/** A string of digits grouped by thousands with commas (`7668402` → `7,668,402`), whatever the locale. */
function groupThousands(digits: string): string {
  let out = '';
  for (let i = 0; i < digits.length; i++) {
    if (i > 0 && (digits.length - i) % 3 === 0) out += ',';
    out += digits[i];
  }
  return out;
}

/**
 * A genomic position or length in bases, in the unit a ruler's tick step chooses: bases below a step of 1,000 (the
 * number grouped, `7,668,402`), kilobases below a step of 10,000 (`7,668 kb`), else megabases with the decimals the
 * step needs (`7.67 Mb` at a step of 10,000, `7.5 Mb` at 500,000, `8 Mb` at 1,000,000). The labels of one step share
 * one unit and one count of decimals, so a ruler reads evenly.
 */
export function formatGenomic(n: number, step: number): string {
  const s = Math.abs(step);
  const sign = n < 0 ? '-' : '';
  const a = Math.abs(n);
  if (s < 1e3) return sign + groupThousands(String(Math.round(a)));
  if (s < 1e4) return `${sign}${groupThousands(String(Math.round(a / 1e3)))} kb`;
  const decimals = Math.max(0, 6 - Math.floor(Math.log10(s) + 1e-9));
  const [int, frac] = (a / 1e6).toFixed(decimals).split('.');
  return `${sign}${groupThousands(int)}${frac ? `.${frac}` : ''} Mb`;
}
