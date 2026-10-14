/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL's 3M-position bench (seqspec:seqql#performance, "The two scale cases"; QR1): a chromosome of 3,000,000
// nucleotides with a dense GC track (isochores of 20-200 kb, noisy per position) and 30,000 repeats, queried with QR1's
// symbols. Targets: a label range under 1 ms; `(values :of gc :above 0.6)` under 15 ms; `and` with `not` linear in
// runs. `expand` (QR2a) on the GC runs; `letters` (under 20 ms) and `search` (under 15 ms), QR2b's. Run with
// `pnpm --filter @molstar/seq-core/model bench`.

import { checkText } from "@molstar/seq-core/query";
import { bench, describe } from "vitest";
import { prng } from "../../_test/model/support/generators.ts";
import { Extents, SeqSets } from "../../src/model/index.ts";
import { compileQuery, createQueryEnv, SeqSelection } from "../../src/model/query/index.ts";

const N = 3_000_000;
const rand = prng(20261012);
// Nucleotides with GC drifting over isochores of 20-200 kb (0.35-0.65), as the genome story's data; GC per position
// is its 1 kb window's fraction, as the story's track (genome-fixture.ts), and, for the worst case, noisy per position.
const code = new Uint8Array(N);
const noisy = new Float64Array(N);
for (let i = 0, base = 0.45, next = 0; i < N; i++) {
  if (i === next) {
    base = 0.35 + rand() * 0.3;
    next = i + 20_000 + Math.floor(rand() * 180_000);
  }
  const r = rand();
  code[i] = r < base ? (r < base / 2 ? 67 : 71) : r < (1 + base) / 2 ? 65 : 84;
  noisy[i] = Math.min(1, Math.max(0, base + (rand() - 0.5) * 0.3));
}
const gc = new Float64Array(N);
for (let w = 0; w < N; w += 1000) {
  let n = 0;
  for (let i = w; i < w + 1000; i++) n += +(code[i] === 67 || code[i] === 71);
  gc.fill(n / 1000, w, w + 1000);
}
const REPEATS = 30_000;
const at = Int32Array.from({ length: REPEATS }, (_, k) => Math.floor((k * N) / REPEATS) + Math.floor(rand() * 50));
const { set } = SeqSets.build({
  molecules: [{ key: "chr1", residues: new TextDecoder().decode(code), alphabet: "dna", ownUnit: true }],
  annotations: [
    { key: "gc", geometry: "values", target: "chr1", columns: { value: gc } },
    { key: "noisy", geometry: "values", target: "chr1", columns: { value: noisy } },
    {
      key: "repeats",
      geometry: "ranges",
      target: "chr1",
      columns: { start: at, stop: Int32Array.from(at, (s) => s + 20 + Math.floor(rand() * 300)) },
    },
  ],
});
const env = createQueryEnv([{ namespace: "app", symbols: [], classes: { iupac: { N: "ACGT", R: "AG", Y: "CT" } } }]);
const query = (text: string) => {
  const q = compileQuery(checkText(text, env.language).expr as never, env);
  const r = q.run({ set });
  const s = r.selection;
  console.log(`${text}: ${s.groupCount} groups, ${s.runCount} runs; ${SeqSelection.summary(s, { maxParts: 1 })}`);
  return () => {
    q.run({ set });
  };
};

describe("SeqQL on 3,000,000 positions", () => {
  const cases = [
    '(sel "chr1[1000000-1200000]")',
    '(sel "chr1@1000000:1200000")',
    "(values :of gc :above 0.6)",
    "(values :of noisy :above 0.6)",
    "(and (values :of gc :above 0.6) (not repeats))",
    "(and chr1 repeats)",
    "(items :of repeats)",
    "(range (values :of gc :above 0.6) 1-1500000)",
    "(or (values :of gc :above 0.6) repeats)",
    "(not repeats :in chr1)",
    // QR2a: the "linear in runs" case
    "(expand (values :of gc :above 0.6) 50)",
    "(and (expand (values :of gc :above 0.6) 50) (not repeats))",
    // QR2b: the scans
    "(letters G)",
    "(letters GC)",
    "(and (values :of gc :above 0.6) (letters G))",
    "(search GAATTC)",
    "(search GAATTC :edits 1)",
    "(search GANTTC :classes iupac)",
  ];
  for (const text of cases) bench(text, query(text), { time: 300 });
  const wide = compileQuery(checkText("(values :of noisy :above 0.6)").expr as never, env).run({ set }).selection;
  bench("toLoci of the GC runs", () => void SeqSelection.toLoci(wide), { time: 300 });
  const ext = Extents.create(1 << 20);
  const toExtents = () => {
    ext.count = 0;
    SeqSelection.toExtents(wide, ext);
  };
  bench("toExtents of the GC runs", toExtents, { time: 300 });
  bench("summary of the GC runs", () => void SeqSelection.summary(wide, { indices: true }), { time: 300 });
});
