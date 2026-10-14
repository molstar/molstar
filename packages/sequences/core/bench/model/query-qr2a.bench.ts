/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QR2a's bench (seqspec:seqql#performance, "The two scale cases"): `onto` over a 100,000-row alignment (the person's
// selection onto it costs the positions it carries; a column range onto every row), across a codon mapping of a
// 1,000,000-residue protein over 3,000,000 nucleotides, and `partners` over a million pairs; `expand` on the genome's
// GC runs. Run with `pnpm --filter @molstar/seq-core/model bench`.

import { checkText } from "@molstar/seq-core/query";
import { bench, describe } from "vitest";
import { prng } from "../../_test/model/support/generators.ts";
import { SeqSets } from "../../src/model/index.ts";
import { compileQuery, createQueryEnv, SeqSelection } from "../../src/model/query/index.ts";

const rand = prng(20261013);
const ROWS = 100_000;
const COLS = 300;
const AA = "ACDEFGHIKLMNPQRSTVWY";
/** Rows of 300 columns, about one gap in eight; row 0 is bound to the unit `ref`. */
const gapped = Array.from({ length: ROWS }, () => {
  let s = "";
  for (let c = 0; c < COLS; c++) s += rand() < 0.125 ? "-" : AA[Math.floor(rand() * 20)];
  return s;
});
const ref = gapped[0].replaceAll("-", "");
const P = 1_000_000;
const PAIRS = 1_000_000;
const U = 100_000;
const { set } = SeqSets.build({
  molecules: [
    { key: "ref", residues: ref, ownUnit: true },
    { key: "prot", residues: "MA".repeat(P / 2), ownUnit: true },
    { key: "gene", residues: "ATG".repeat(P), ownUnit: true, alphabet: "dna" },
    { key: "x", residues: "A".repeat(U), ownUnit: true },
    { key: "y", residues: "A".repeat(U), ownUnit: true },
  ],
  alignments: [{ key: "msa", gapped, flavor: "aligned" }],
  mappings: [
    { key: "bind", from: "msa/0", to: "ref", method: "identity" },
    { key: "codon", from: "prot", to: "gene", method: "codon", start: 0 },
  ],
  annotations: [
    {
      key: "pv",
      geometry: "values",
      target: "prot",
      columns: { value: Float64Array.from({ length: P }, (_, i) => (i % 10 < 5 ? 1 : 0)) },
    },
    // a million contacts between x and y: ten per residue of y, to residues of x near it
    {
      key: "pp",
      geometry: "pairs",
      columns: {
        target: Array(PAIRS).fill("x"),
        start: Int32Array.from({ length: PAIRS }, (_, i) => Math.min(U - 1, Math.floor(i / 10) + (i % 7))),
        b_target: Array(PAIRS).fill("y"),
        b_start: Int32Array.from({ length: PAIRS }, (_, i) => Math.floor(i / 10)),
      },
    },
  ],
});
const env = createQueryEnv();
const selection = SeqSelection.fromExtents(set, {
  count: 1,
  kind: Uint8Array.of(1),
  id: Int32Array.of(0),
  row: Int32Array.of(-1),
  start: Int32Array.of(0),
  stop: Int32Array.of(ref.length),
});
const query = (text: string) => {
  const q = compileQuery(checkText(text).expr as never, env);
  const t0 = performance.now();
  const s = q.run({ set, inputs: { selection } }).selection;
  const first = performance.now() - t0;
  const say = `${s.groupCount} groups, ${s.runCount} runs, first run ${first.toFixed(1)} ms`;
  console.log(`${text}: ${say}; ${SeqSelection.summary(s, { maxParts: 1 })}`);
  return () => {
    q.run({ set, inputs: { selection } });
  };
};

describe("onto and partners at scale", () => {
  const cases = [
    "(onto (input selection) msa)",
    '(onto (sel "ref[1-50]") msa :kind row)',
    '(onto (sel "msa[100-120]") msa/0)',
    "(onto msa ref)",
    '(onto (sel "prot[200000-700000]") gene)',
    "(onto (values :of pv :above 0.5) gene)",
    '(onto (sel "gene[1000000-1500000]") prot)',
    "(partners y :via pp)",
    '(partners (sel "y[1-10000]") :via pp)',
    "(and x (partners y :via pp :side a))",
    "(expand (values :of pv :above 0.5) 3)",
  ];
  for (const text of cases) bench(text, query(text), { time: 300 });
});
