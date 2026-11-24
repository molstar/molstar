/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kernels of exact and class patterns (_spec/README.md "Search", Kernels): Horspool for an exact pattern
// (sublinear), shift-and over a table of 256 masks for a pattern of sets (a word every 32 positions). Bytes only: what a
// code stands for is the matcher's.

import type { HitVisitor, SearchKernel, StepBudget } from './types.js';

/** Horspool's search for one exact pattern of codes (copied). */
export function exactKernel(pattern: Uint8Array): SearchKernel {
  const pat = pattern.slice();
  const m = pat.length;
  const shift = new Int32Array(256).fill(m);
  for (let i = 0; i < m - 1; i++) shift[pat[i]] = m - 1 - i;
  const last = pat[m - 1];
  return {
    maxLength: m,
    // a shift a letter at most, charged before; a check's comparisons (O(n·m) at worst: `C` and 999 `A` over poly-A)
    // added as they are made
    words: 1,
    find(seq: Uint8Array, n: number, from: number, to: number, _edges: number, budget: StepBudget, visit: HitVisitor) {
      if (m === 0) return true;
      const end = Math.min(to, n - m + 1);
      for (let s = from; s < end; ) {
        const c = seq[s + m - 1];
        if (c === last) {
          let k = m - 2;
          while (k >= 0 && seq[s + k] === pat[k]) k--;
          if (k < 0) visit(s, s + m);
          budget.steps += m - 1 - k;
          if (budget.steps > budget.limit) return false;
        }
        s += shift[c];
      }
      return true;
    },
  };
}

/**
 * Shift-and for a pattern of sets: position `i` accepts code `c` when bit `i` of `masks[c]` is set, `masks` a table of
 * 256 codes times `ceil(m / 32)` words (code-major). Hits in the order of their ends.
 */
export function classKernel(masks: Int32Array, m: number): SearchKernel {
  const w = Math.ceil(m / 32);
  const top = w - 1;
  const accept = 1 << ((m - 1) & 31);
  const d = new Int32Array(w);
  return {
    maxLength: m,
    words: w,
    find(seq: Uint8Array, n: number, from: number, to: number, _edges: number, _budget: StepBudget, visit: HitVisitor) {
      if (m === 0) return true;
      d.fill(0);
      // the last end of a hit that starts before `to`
      const stop = Math.min(n, to + m - 1);
      if (w === 1) {
        let x = 0;
        for (let u = from; u < stop; u++) {
          x = ((x << 1) | 1) & masks[seq[u]];
          if ((x & accept) !== 0) visit(u + 1 - m, u + 1);
        }
        return true;
      }
      for (let u = from; u < stop; u++) {
        const base = seq[u] * w;
        let carry = 1;
        for (let k = 0; k < w; k++) {
          const v = d[k];
          d[k] = ((v << 1) | carry) & masks[base + k];
          carry = v >>> 31;
        }
        if ((d[top] & accept) !== 0) visit(u + 1 - m, u + 1);
      }
      return true;
    },
  };
}

/** A table of masks for `classKernel` from each position's accepting test over the 256 codes. */
export function classMasks(m: number, accepts: (position: number, code: number) => boolean): Int32Array {
  const w = Math.ceil(m / 32);
  const masks = new Int32Array(256 * w);
  for (let c = 0; c < 256; c++) for (let i = 0; i < m; i++) if (accepts(i, c)) masks[c * w + (i >> 5)] |= 1 << (i & 31);
  return masks;
}
