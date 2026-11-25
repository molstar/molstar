/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The search tests' helpers: targets over strings, a run collected to its end, seeded letters, and the naive oracle the
// matchers are checked against.

import type { SequenceAlphabet } from '@molstar/core/util/sequence-constants';
import type { Issue } from '@molstar/core/repr/tree';
import {
  compileFailed,
  type SearchBatch,
  type SearchBounds,
  type SearchMatcher,
  type SearchOptions,
  type SearchTarget,
  searchBounds,
} from '../index.js';

/** A target over a string's letters (read counted in `reads`). */
export function textTarget(
  letters: string,
  o: { system?: string; alphabet?: SequenceAlphabet; cyclic?: boolean; within?: number[]; reads?: number[][] } = {},
): SearchTarget {
  const codes = Uint8Array.from(letters, (c) => c.charCodeAt(0));
  return {
    system: o.system ?? '/data/s',
    alphabet: o.alphabet ?? 'dna',
    length: codes.length,
    cyclic: o.cyclic ?? false,
    ...(o.within ? { within: Int32Array.from(o.within) } : {}),
    letters: {
      read: async (start, end) => {
        o.reads?.push([start, end]);
        return codes.slice(start, end);
      },
    },
  };
}

export interface Ran {
  readonly hits: [number, number, number, number][];
  readonly counted: number;
  readonly issues: Issue[];
  readonly batches: SearchBatch[];
}

/** Compiles and runs a pattern to its end; throws the compile's issues. */
export async function runSearch(
  matcher: SearchMatcher<any>,
  pattern: string,
  targets: readonly SearchTarget[],
  options: SearchOptions = {},
  bounds: Partial<SearchBounds> = {},
  signal: AbortSignal = new AbortController().signal,
): Promise<Ran> {
  const b = searchBounds(bounds);
  const compiled = matcher.compile(pattern, options, b);
  if (compileFailed(compiled)) throw new Error(compiled.map((i) => i.code).join(','));
  const hits: [number, number, number, number][] = [];
  const issues: Issue[] = [];
  const batches: SearchBatch[] = [];
  let counted = 0;
  await new Promise<void>((resolve, reject) => {
    matcher.run(compiled, targets, { signal, bounds: b }).subscribe({
      next: (batch) => {
        batches.push(batch);
        counted += batch.counted;
        issues.push(...batch.issues);
        for (let i = 0; i < batch.hits.length; i++)
          hits.push([batch.hits.target[i], batch.hits.start[i], batch.hits.end[i], batch.hits.strand[i]]);
      },
      error: reject,
      complete: resolve,
    });
  });
  return { hits, counted, issues, batches };
}

/** Seeded letters (a linear congruential generator): reproducible across runs. */
export function seeded(seed: number): () => number {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
}

export function randomLetters(rand: () => number, n: number, alphabet: string): string {
  let out = '';
  for (let i = 0; i < n; i++) out += alphabet[Math.floor(rand() * alphabet.length)];
  return out;
}

const Comp: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', R: 'Y', Y: 'R', S: 'S', W: 'W' };
export const revcomp = (s: string) =>
  [...s.toUpperCase().replace(/U/g, 'T')]
    .reverse()
    .map((c) => Comp[c] ?? c)
    .join('');

/**
 * The naive oracle: every start of a target, the letters folded (upper case, U as T, gaps kept out of a hit's letters),
 * `accepts(window)` the longest accepted length from that start (0: none), each strand, palindromes once.
 */
export function oracle(
  letters: string,
  longestAt: (text: string, start: number) => number,
  o: { cyclic?: boolean; strands?: 'both' | 'forward'; protein?: boolean } = {},
): [number, number, number, number][] {
  const L = letters.length;
  const text = letters.toUpperCase().replace(/U/g, 'T');
  const unrolled = o.cyclic ? text + text : text;
  const out: [number, number, number, number][] = [];
  const plus = new Map<number, number>();
  for (let s = 0; s < L; s++) {
    const len = longestAt(unrolled, s);
    if (len > 0 && len <= L) plus.set(s, len);
  }
  const minus: [number, number][] = [];
  if (!o.protein && (o.strands ?? 'both') === 'both') {
    const rc = revcomp(unrolled);
    const n = unrolled.length;
    for (let s = 0; s < n; s++) {
      const len = longestAt(rc, s);
      if (len <= 0 || len > L) continue;
      // forward [n - s - len, n - s)
      const fs = n - s - len;
      if (fs >= 0 && fs < L) minus.push([fs, fs + len]);
    }
  }
  const all: [number, number, number, number][] = [];
  for (const [s, len] of plus) all.push([0, s, s + len, o.protein ? 0 : 1]);
  for (const [s, e] of minus) all.push([0, s, e, -1]);
  const rank = (x: number) => (x === 1 ? 0 : x === 0 ? 1 : 2);
  all.sort((a, b) => a[1] - b[1] || a[2] - b[2] || rank(a[3]) - rank(b[3]));
  // de-duplicate (a cyclic minus hit seen on both laps) and merge palindromes
  for (const h of all) {
    const prev = out[out.length - 1];
    if (prev && prev[1] === h[1] && prev[2] === h[2]) {
      if (prev[3] === 1 && h[3] === -1) prev[3] = 0;
      continue;
    }
    out.push([...h]);
  }
  return out;
}
