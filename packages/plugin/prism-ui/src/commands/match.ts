/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Subsequence match, case-insensitive: letters in order; start, word-start and contiguous hits score more. */
export function fuzzyScore(text: string, query: string): number | undefined {
  const t = text.toLowerCase();
  const q = query.toLowerCase().trim();
  if (!q) return 0;
  let score = 0;
  let j = 0;
  let last = -2;
  for (let i = 0; i < t.length && j < q.length; i++) {
    if (t[i] !== q[j]) continue;
    score += i === last + 1 ? 3 : 1;
    if (i === 0) score += 4;
    else if (/[\s._:-]/.test(t[i - 1] as string)) score += 2;
    last = i;
    j++;
  }
  return j === q.length ? score - t.length * 0.01 : undefined;
}

/** Items that match a query, best first; all of them, in order, for an empty query. */
export function rankBy<T>(items: readonly T[], query: string, text: (item: T) => string): T[] {
  if (!query.trim()) return [...items];
  return items
    .map((item) => ({ item, s: fuzzyScore(text(item), query) }))
    .filter((x): x is { item: T; s: number } => x.s !== undefined)
    .sort((a, b) => b.s - a.s)
    .map((x) => x.item);
}

/** Where `fuzzyScore` found the query's letters in the text (the same greedy pass), for highlighting; empty if none. */
export function matchPositions(text: string, query: string): number[] {
  const t = text.toLowerCase();
  const q = query.toLowerCase().trim();
  const at: number[] = [];
  let j = 0;
  for (let i = 0; i < t.length && j < q.length; i++) {
    if (t[i] !== q[j]) continue;
    at.push(i);
    j++;
  }
  return j === q.length ? at : [];
}
