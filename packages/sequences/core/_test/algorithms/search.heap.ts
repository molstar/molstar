/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the search matchers (SR1; seqspec:libraries#search "Tests"): a search over a text with hundreds of
// hits, repeated until a million hits have been visited, on RB0b's stepped guard: runs of `n` and `4n` steps, the
// longer growing at most 64 KB beyond the shorter, and no GC in either. Untimed, as every core guard: the steps run in
// `heapScaling`'s own warm-up, by count. `regexMatcher` allocates a match array per hit and is not guarded. The cases
// are search-cases.ts's, shared with the quick twin (search.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../model/heap.js';
import { SEARCH_CASES, searchStep } from './search-cases.js';

/** Steps in the shorter run; the longer runs four times as many. */
const N = 2000;

describe('search heap guards', () => {
  for (const c of SEARCH_CASES) {
    it(c.name, () => {
      const { step, hits } = searchStep(c);
      const r = heapScaling(N, step);
      expect(hits(), 'the guard visits hits').toBeGreaterThan(N * 200);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      // The template's control (CB-05): the same step with a leak is caught.
      expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
    }, 60_000);
  }
});
