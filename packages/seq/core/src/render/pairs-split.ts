/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PairsData } from '@molstar/seqx-core/resolve';

// The pairs layer's `split` (_spec/README.md "Render", The pairs layer): the largest nested set of pairs on `side`, the
// rest on the other side; exact to `SplitExact` pairs, greedy past it.

/** Pairs a `split` places exactly (the largest nested set: quadratic); past it, greedily (fixed: no figure asks more). */
export const SplitExact = 2000;

/**
 * The pairs drawn on the other side under `split`: the complement of the largest nested (non-crossing) set of the
 * pairs as intervals between their ends' centres; exact by the interval recurrence to `SplitExact` pairs, greedy past
 * it (narrowest first, kept when it crosses none kept).
 */
export function splitCrossing(d: PairsData, shown?: Uint8Array): Uint8Array {
  // the pairs drawn (a filter's): those left out cross nothing
  const at: number[] = [];
  for (let k = 0; k < d.n; k++) if (!shown || shown[k]) at.push(k);
  const n = at.length;
  const lo = new Float64Array(n);
  const hi = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const k = at[i];
    const ca = (d.aStart[k] + d.aEnd[k]) / 2;
    const cb = (d.bStart[k] + d.bEnd[k]) / 2;
    lo[i] = Math.min(ca, cb);
    hi[i] = Math.max(ca, cb);
  }
  const keep = new Uint8Array(n);
  if (n <= SplitExact) nestedExact(n, lo, hi, keep);
  else nestedGreedy(n, lo, hi, keep);
  const other = new Uint8Array(d.n);
  for (let i = 0; i < n; i++) other[at[i]] = keep[i] ? 0 : 1;
  return other;
}

/** Whether two intervals cross: one's ends either side of one end of the other. */
export function crosses(alo: number, ahi: number, blo: number, bhi: number): boolean {
  return (alo < blo && blo < ahi && ahi < bhi) || (blo < alo && alo < bhi && bhi < ahi);
}

/**
 * The largest set of non-crossing intervals (nested or apart): over the 2n ends in order, F(i, j) the most held within
 * ends i..j, F(i, j) = max(F(i + 1, j), 1 + F(i + 1, k − 1) + F(k + 1, j)) where k is i's partner in (i, j]. At one
 * place closing ends come first, then each self-pair's two ends, then opening ends; a wider interval outside a narrower
 * one and equal ones nested, as `crosses` reads them (touching ends do not cross). O(n²) time, a triangle of int16.
 */
export function nestedExact(n: number, lo: Float64Array, hi: Float64Array, keep: Uint8Array) {
  const m = 2 * n;
  const ends = new Int32Array(m);
  for (let i = 0; i < m; i++) ends[i] = i;
  const at = (e: number) => (e < n ? lo[e] : hi[e - n]);
  // at one place: the others' closing ends, each self-pair's two ends together (it crosses nothing), the others'
  // opening ends; equal intervals nest (opened in order, closed in reverse)
  const rank = (e: number) => (lo[e % n] === hi[e % n] ? 1 : e < n ? 2 : 0);
  ends.sort((x, y) => {
    const d = at(x) - at(y);
    if (d) return d;
    const r = rank(x) - rank(y);
    if (r) return r;
    const px = x % n;
    const py = y % n;
    if (r === 0 && rank(x) === 1) return px !== py ? px - py : x < n ? -1 : 1;
    const w = hi[px] - lo[px] - (hi[py] - lo[py]);
    if (w) return x < n ? -w : w;
    return x < n ? px - py : py - px;
  });
  const pos = new Int32Array(m);
  for (let i = 0; i < m; i++) pos[ends[i]] = i;
  const partner = new Int32Array(m);
  for (let i = 0; i < m; i++) {
    const e = ends[i];
    partner[i] = pos[e < n ? e + n : e - n];
  }
  // F over the triangle i ≤ j, row i from `row(i)`
  const row = (i: number) => i * m - (i * (i - 1)) / 2 - i;
  const F = new Int16Array((m * (m + 1)) / 2);
  const get = (i: number, j: number) => (i > j ? 0 : F[row(i) + j]);
  for (let len = 2; len <= m; len++)
    for (let i = 0; i + len - 1 < m; i++) {
      const j = i + len - 1;
      let best = get(i + 1, j);
      const k = partner[i];
      if (k > i && k <= j) {
        const v = 1 + get(i + 1, k - 1) + get(k + 1, j);
        if (v > best) best = v;
      }
      F[row(i) + j] = best;
    }
  // back through the choices
  const stack: number[] = [0, m - 1];
  while (stack.length) {
    const j = stack.pop()!;
    let i = stack.pop()!;
    while (i < j) {
      const k = partner[i];
      if (k > i && k <= j && get(i, j) === 1 + get(i + 1, k - 1) + get(k + 1, j)) {
        keep[ends[i] % n] = 1;
        stack.push(i + 1, k - 1);
        i = k + 1;
      } else i++;
    }
  }
}

/** Narrowest first, a pair kept when it crosses no pair kept. */
export function nestedGreedy(n: number, lo: Float64Array, hi: Float64Array, keep: Uint8Array) {
  const order = Array.from({ length: n }, (_, k) => k).sort((x, y) => hi[x] - lo[x] - (hi[y] - lo[y]));
  const kept: number[] = [];
  for (const k of order) {
    let ok = true;
    for (let i = 0; i < kept.length && ok; i++) if (crosses(lo[k], hi[k], lo[kept[i]], hi[kept[i]])) ok = false;
    if (ok) {
      kept.push(k);
      keep[k] = 1;
    }
  }
}
