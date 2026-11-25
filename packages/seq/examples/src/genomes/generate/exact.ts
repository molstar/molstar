/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The generator's smooth shapes in arithmetic every engine rounds alike: integers, `+`, `-`, `*` and `/` (IEEE 754
// rounds each correctly), never `Math.sin`, `Math.exp`, `Math.pow` or `**` (an implementation's approximations, not
// required to agree across engines), so the fixtures are byte-identical wherever the generator runs.

/**
 * A wave of `period` positions at the integer `p`, in [-1, 1]: two parabolas a period (up, then down), within 6% of a
 * sine of that period. The phase is exact: `p mod period` on integers.
 */
export function wave(p: number, period: number): number {
  const m = ((p % period) + period) % period;
  const u = (2 * m) / period;
  return u < 1 ? 4 * u * (1 - u) : -4 * (u - 1) * (2 - u);
}

/** `x` to the integer power `n` ≥ 0, by multiplication. */
export function powi(x: number, n: number): number {
  let r = 1;
  for (let k = 0; k < n; k++) r *= x;
  return r;
}

/** ln 10, written out (`Math.LN10`'s value). */
const Ln10 = 2.302585092994046;

/**
 * 10 to the `w` ≥ 0: its integer part by multiplication, its fraction by e^(f ln 10)'s Taylor series to 16 terms (f ln
 * 10 < 2.31: the last term below 1e-8 of the sum).
 */
export function pow10(w: number): number {
  const n = Math.floor(w);
  const y = (w - n) * Ln10;
  let term = 1;
  let sum = 1;
  for (let k = 1; k <= 16; k++) {
    term = (term * y) / k;
    sum += term;
  }
  return sum * powi(10, n);
}
