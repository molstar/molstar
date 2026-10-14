/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for SeqQL's runtime (seqspec:seqql#performance, QR1): a query allocates its typed result arrays (off the
// JS heap), its scratch and one closure tree at compile, never an object per run or per position. Each guard runs a
// query over 100,000 runs or groups (a million positions scanned) and allows 1 MB of heap growth: 10 bytes a run.

import { checkText } from "@molstar/seq-core/query";
import { describe, expect, it } from "vitest";
import { SeqSets } from "../../src/model/index.ts";
import { compileQuery, createQueryEnv, SeqSelection } from "../../src/model/query/index.ts";
import { controlCaught, heapGrowth, LIMIT } from "./heap.ts";
import { leakExtension } from "./query-leak.ts";

const N = 1_000_000;
const { set } = SeqSets.build({
  molecules: [{ key: "big", residues: "ACGT".repeat(N / 4), ownUnit: true, alphabet: "dna" }],
  annotations: [
    // passes in runs of 5 every 10: 100,000 runs
    {
      key: "v",
      geometry: "values",
      target: "big",
      columns: { value: Float64Array.from({ length: N }, (_, i) => (i % 10 < 5 ? 1 : 0)) },
    },
    // 50,000 items of 4 positions, every 20
    {
      key: "reps",
      geometry: "ranges",
      target: "big",
      columns: {
        start: Int32Array.from({ length: N / 20 }, (_, i) => 20 * i + 3),
        stop: Int32Array.from({ length: N / 20 }, (_, i) => 20 * i + 7),
      },
    },
  ],
});
const env = createQueryEnv([leakExtension()]);
const SELECTION = SeqSelection.fromExtents(set, {
  count: 1,
  kind: Uint8Array.of(1),
  id: Int32Array.of(0),
  row: Int32Array.of(-1),
  start: Int32Array.of(1000),
  stop: Int32Array.of(900_000),
});

const QUERIES = [
  "(values :of v :above 0.5)",
  "(and (values :of v :above 0.5) (not reps))",
  "(and (values :of v :above 0.5) reps (input selection))",
  "(range (values :of v :above 0.5) 1-500000 600000-900000)",
  "(slice (values :of v :above 0.5) 0 500000 600000 900000)",
  "(or reps (values :of v :above 0.5) reps)",
  "(not reps :in big)",
  "(items :of reps)",
];

/** The least of three measured runs after ten warm-ups, and the groups selected. */
function measure(text: string): { groups: number; growth: number } {
  const q = compileQuery(checkText(text, env.language).expr as never, env);
  let groups = 0;
  const run = () => {
    groups = q.run({ set, inputs: { selection: SELECTION } }).selection.groupCount;
  };
  for (let k = 0; k < 10; k++) run();
  // The least of three measured runs: optimized code that arrives late from a background compile (on a loaded
  // machine) can make one run box a few values; an allocation per run shows in every run.
  return { groups, growth: Math.min(heapGrowth(run), heapGrowth(run), heapGrowth(run)) };
}

describe("SeqQL's runtime allocates nothing per run", () => {
  for (const text of QUERIES) {
    it(text, () => {
      const { groups, growth } = measure(text);
      expect(groups).toBeGreaterThan(text.startsWith("(not") ? 0 : 40_000);
      expect(growth).toBeLessThan(LIMIT);
    });
  }

  // CM-05: the template's control, an object a run through `(leak Q)`, caught by the same measurement.
  it("an object a run is caught (the positive control)", () => {
    for (const text of [QUERIES[0], QUERIES[1], QUERIES[7]])
      expect([text, controlCaught(() => measure(`(leak ${text})`).growth)]).toEqual([text, true]);
  });
});
