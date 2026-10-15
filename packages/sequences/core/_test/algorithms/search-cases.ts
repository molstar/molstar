/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The search matchers' cases, for their heap guard and its quick twin (search.heap.ts, search.quick.ts): a search over a
// text with hundreds of hits, a step being one pass of one matcher over it. `regexMatcher` allocates a match array per
// hit and is not guarded.

import {
  approximateMatcher,
  classMatcher,
  compileMotif,
  exactMatcher,
  type Matcher,
  toCodes,
} from '../../src/algorithms/search/index.js';

/** A text of 2,000 letters in which every pattern below has 250 to 500 hits. */
export const seq = toCodes('ACGTACGAAC'.repeat(200));

let sum = 0;
const visit = (start: number, stop: number, score: number) => {
  sum = (sum + start + stop + score) & 0xfffff;
};

export interface SearchCase {
  name: string;
  matcher: Matcher;
  cyclic: boolean;
}

export const SEARCH_CASES: SearchCase[] = [
  { name: 'exact, linear and across the end of a circle', matcher: exactMatcher(toCodes('AC')), cyclic: true },
  {
    name: 'class (a degenerate motif), linear',
    matcher: classMatcher(compileMotif('ARG', { R: 'CGT' })),
    cyclic: false,
  },
  {
    name: 'class across the end of a circle, 40 positions (two words)',
    matcher: classMatcher(compileMotif('ACGTACGAAC'.repeat(4), {})),
    cyclic: true,
  },
  { name: 'approximate, one word, 1 edit', matcher: approximateMatcher(toCodes('ACGAA'), 1), cyclic: false },
  {
    name: 'approximate, across the end of a circle',
    matcher: approximateMatcher(toCodes('GTACGAAC'), 2),
    cyclic: true,
  },
  {
    name: 'approximate, several words, 2 edits',
    matcher: approximateMatcher(toCodes('ACGTACGAAC'.repeat(5)), 2),
    cyclic: false,
  },
  {
    name: 'approximate, every end (ends: all), 2 edits',
    matcher: approximateMatcher(toCodes('ACGAAC'), 2, { ends: 'all' }),
    cyclic: true,
  },
];

/** A case's step (one search over the text) and what it has visited, so a guard over a step that finds nothing fails. */
export function searchStep({ matcher, cyclic }: SearchCase) {
  let hits = 0;
  let calls = 0;
  const step = () => {
    calls++;
    hits += matcher.find(seq, 0, seq.length, cyclic, visit);
  };
  return { step, hits: () => hits, calls: () => calls };
}
