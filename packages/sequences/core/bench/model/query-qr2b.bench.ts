/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QR2b's bench (seqspec:seqql#performance, "The two scale cases"): a 100,000-row alignment of 600 columns with a rows
// key table of identities. `(rows :of msa :where (>= (field identity :of stats) 0.9))` evaluates 100,000 system tests
// (target: under 10 ms); `(and (rows :of msa) (letters C))` would visit about 5·10^7 cells and stops at the positions
// limit with its message; `positions` tests a million positions. Run with `pnpm --filter @molstar/seq-core/model bench`.

import { checkText, renderIssue } from "@molstar/seq-core/query";
import { bench, describe } from "vitest";
import { prng } from "../../_test/model/support/generators.js";
import { SeqSets } from "../../src/model/index.js";
import { compileQuery, createQueryEnv } from "../../src/model/query/index.js";

const rand = prng(20261014);
const ROWS = 100_000;
const COLS = 600;
const AA = "ACDEFGHIKLMNPQRSTVWY";
const gapped = Array.from({ length: ROWS }, () => {
  let s = "";
  for (let c = 0; c < COLS; c++) s += rand() < 0.125 ? "-" : AA[Math.floor(rand() * 20)];
  return s;
});
const { set } = SeqSets.build({
  molecules: [],
  alignments: [{ key: "msa", gapped, flavor: "aligned" }],
  annotations: [
    {
      key: "stats",
      geometry: "rows",
      target: "msa",
      columns: {
        row: Int32Array.from({ length: ROWS }, (_, r) => r),
        identity: Float64Array.from({ length: ROWS }, () => rand()),
      },
    },
  ],
});
const env = createQueryEnv();
const query = (text: string) => {
  const q = compileQuery(checkText(text).expr as never, env);
  const t0 = performance.now();
  const r = q.run({ set });
  const first = performance.now() - t0;
  const s = r.selection;
  const said = r.issues.map((i) => renderIssue(i).split("\n")[0]).join(" | ");
  console.log(`${text}: ${s.groupCount} groups, ${s.runCount} runs, first run ${first.toFixed(1)} ms; ${said}`);
  return () => {
    q.run({ set });
  };
};

describe("system tests and scans over 100,000 rows", () => {
  const cases = [
    "(rows :of msa :where (>= (field identity :of stats) 0.9))",
    "(rows :of msa :where (> (length) 530))",
    "(and (rows :of msa) (letters C))",
    "(and (rows :of msa :where (>= (field identity :of stats) 0.99)) (letters C))",
    "(positions :in (rows :of msa :where (>= (field identity :of stats) 0.98)) :where (= (letter) W))",
    "(search WW :in (rows :of msa :where (>= (field identity :of stats) 0.9)))",
  ];
  for (const text of cases) bench(text, query(text), { time: 300 });
});
