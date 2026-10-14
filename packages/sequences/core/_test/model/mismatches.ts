/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An oracle loop's checks, for every package: thousands of `expect`s, one a case, were most of such a test's time (an
// `expect` costs microseconds; the case it checks, often less). The loop records each case that disagrees, and the test
// checks once that none did: the same checks, and a failure lists the cases that failed.

import { expect } from "vitest";

/** Numbers, strings and arrays of them as `toEqual` compares them: NaN equals NaN, -0 is not 0, arrays by element. */
function same(a: unknown, b: unknown): boolean {
  if (Object.is(a, b)) return true;
  const arrays = (Array.isArray(a) || ArrayBuffer.isView(a)) && (Array.isArray(b) || ArrayBuffer.isView(b));
  if (!arrays || a.constructor !== b.constructor) return false;
  const [x, y] = [a as ArrayLike<unknown>, b as ArrayLike<unknown>];
  if (x.length !== y.length) return false;
  for (let i = 0; i < x.length; i++) if (!same(x[i], y[i])) return false;
  return true;
}

/** The disagreements of a loop: `check` each case, then `none()` once. */
export function mismatches(shown = 20) {
  const found: unknown[] = [];
  let count = 0;
  return {
    /** `got` against `want` (numbers, strings, booleans and arrays of them); `at` names the case. */
    check(got: unknown, want: unknown, at?: unknown) {
      if (same(got, want)) return;
      count++;
      if (found.length < shown) found.push(at === undefined ? { got, want } : { at, got, want });
    },
    /** Fails, listing the first disagreements and how many there were, unless every case agreed. */
    none() {
      expect({ count, first: found }).toEqual({ count: 0, first: [] });
    },
  };
}
