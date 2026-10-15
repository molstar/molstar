/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for SeqQL's runtime (seqspec:seqql#performance, QR1): a query allocates its typed result arrays (off the JS
// heap), its scratch and one closure tree at compile, never an object per run or per position. Each guard runs a query
// over 100,000 runs or groups (a million positions scanned) and allows 1 MB of heap growth: 10 bytes a run. Untimed: the
// runs are warmed by count. The cases are query-cases.ts's, shared with the quick twin (query.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, heapGrowth, LIMIT } from './heap.js';
import { prepare, QUERIES } from './query-cases.js';

/** The least of three measured runs after ten warm-ups, and the groups selected. */
function measure(text: string): { groups: number; growth: number } {
  const { run, groups } = prepare(text);
  for (let k = 0; k < 10; k++) run();
  // The least of three measured runs: optimized code that arrives late from a background compile (on a loaded
  // machine) can make one run box a few values; an allocation per run shows in every run.
  const growth = Math.min(heapGrowth(run), heapGrowth(run), heapGrowth(run));
  return { groups: groups(), growth };
}

describe("SeqQL's runtime allocates nothing per run", () => {
  for (const text of QUERIES) {
    it(text, () => {
      const { groups, growth } = measure(text);
      expect(groups).toBeGreaterThan(text.startsWith('(not') ? 0 : 40_000);
      expect(growth).toBeLessThan(LIMIT);
    });
  }

  // CM-05: the template's control, an object a run through `(leak Q)`, caught by the same measurement.
  it('an object a run is caught (the positive control)', () => {
    for (const text of [QUERIES[0], QUERIES[1], QUERIES[7]])
      expect([text, controlCaught(() => measure(`(leak ${text})`).growth)]).toEqual([text, true]);
  });
});
