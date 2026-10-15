/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for QR2a's runtimes (seqspec:seqql#performance): `onto` through a row's binding and a codon mapping,
// `partners` over 100,000 pairs, `expand`, `group`, `pick`, `first`, `ids` and an item test each run over 100,000 runs
// or 50,000 groups and allow 1 MB of heap growth (10 bytes a run). The positive control is the template's: an object
// a run, through `(leak Q)` (query-leak.ts), caught by the same measurement. `each` is not guarded: it builds a
// selection per group by definition (its `(current)` and E's result). Untimed: the runs are warmed by count. The cases
// are query-qr2a-cases.ts's, shared with the quick twin (query-qr2a.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, heapGrowth, LIMIT } from './heap.js';
import { prepare, QUERIES } from './query-qr2a-cases.js';

/** The least of three measured runs after ten warm-ups (a loaded machine's late compiles can box in one run). */
function growth(text: string): { groups: number; growth: number } {
  const { run, groups } = prepare(text);
  for (let k = 0; k < 10; k++) run();
  const least = Math.min(heapGrowth(run), heapGrowth(run), heapGrowth(run));
  return { groups: groups(), growth: least };
}

describe("QR2a's runtimes allocate nothing per run", () => {
  for (const [text, groups] of QUERIES) {
    it(text, () => {
      const r = growth(text);
      expect(r.groups).toBeGreaterThanOrEqual(groups);
      expect(r.growth).toBeLessThan(LIMIT);
    });
  }

  it("an object a run is caught by the same measurement (the template's control)", () => {
    for (const text of ['(onto (values :of v :above 0.5) aln)', '(partners y :via pp)']) {
      expect([text, controlCaught(() => growth(`(leak ${text})`).growth)]).toEqual([text, true]);
    }
  });
});
