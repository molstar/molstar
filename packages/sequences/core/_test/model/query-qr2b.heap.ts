/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for QR2b's runtimes (seqspec:seqql#performance): `letters`, `positions` with a test per position,
// `search` with each matcher but the regex one (which allocates its match per hit, seqspec:libraries#search), the group
// tests over a selection unioned once, and system tests over a 100,000-row alignment. Each runs over 100,000 runs,
// groups or rows (a million positions scanned or tested) and allows 1 MB of heap growth (10 bytes a run). Label tests
// (`=`, `!=`, `in` on `(label)`, ordinal and through a labels and an offset numbering) compare as the numbering holds
// labels, no string per position (CM-04). The positive control is the template's: an object a run, through `(leak Q)`
// (query-leak.ts), caught by the same measurement. Untimed: the runs are warmed by count. The cases are
// query-qr2b-cases.ts's, shared with the quick twin (query-qr2b.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, heapGrowth, LIMIT } from './heap.js';
import { prepare, QUERIES } from './query-qr2b-cases.js';

/** The least of three measured runs after ten warm-ups (a loaded machine's late compiles can box in one run). */
function growth(text: string): { groups: number; growth: number } {
  const { run, groups } = prepare(text);
  for (let k = 0; k < 10; k++) run();
  const least = Math.min(heapGrowth(run), heapGrowth(run), heapGrowth(run));
  return { groups: groups(), growth: least };
}

describe("QR2b's runtimes allocate nothing per run", () => {
  for (const [text, groups] of QUERIES) {
    it(text, () => {
      const r = growth(text);
      expect(r.groups).toBeGreaterThanOrEqual(groups);
      expect(r.growth).toBeLessThan(LIMIT);
    });
  }

  it("an object a run is caught by the same measurement (the template's control)", () => {
    for (const text of ['(letters A :in big)', '(search ACG :in big)', '(positions :in big :where (= (letter) T))']) {
      expect([text, controlCaught(() => growth(`(leak ${text})`).growth)]).toEqual([text, true]);
    }
  });
});
