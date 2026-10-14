/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The SeqQL table's examples, evaluated (QR2b): every example of every base symbol runs on its `fixture` (a committed
// case, compiled here by the format on the bundled fixtures) to the result recorded below; an example without one uses
// the spec's illustrative refs (domains, variants, msa, fam, hits, sites), which a small set built here holds, so it
// runs too. A new example is added here with its result, or the test fails.

import { type SeqSet, SeqSets, type SetInput } from "@molstar/seq-core/model";
import { compileQuery, createQueryEnv, SeqSelection } from "@molstar/seq-core/model/query";
import { BASE_SYMBOLS, checkText } from "@molstar/seq-core/query";
import { describe, expect, it } from "vitest";
import { compile } from "../../src/format/index.ts";
import { CASES } from "./cases.ts";
import { fixtureLoaders, registries } from "./load-fixtures.ts";

/** The illustrative refs as data: a synthetic chain `1brs/A` with cysteines, insertion codes and a GXGXXG motif. */
function illustrative(): SetInput {
  const A = "MKGCAGCQGCSCGXGXXGLLKTAVDEWHRAGSGAAGTPLKCQKVLNEFGHAPDRSTVWAKLQNGEGAGTTGALDRSLQYVKPHCEEFTWLRNDKSGAEHMTCAA";
  const labels = Array.from({ length: A.length }, (_, i) => (i === 30 ? "30A" : i === 31 ? "30B" : i < 30 ? i + 1 : i));
  return {
    molecules: [
      { key: "1brs/entity:1", residues: A, numberings: [{ name: "auth", method: "labels", labels, isDefault: true }] },
      {
        key: "1brs/entity:2",
        residues: "KKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEG",
      },
    ],
    units: [
      { key: "1brs/A", molecule: "1brs/entity:1", score: 0.9 },
      { key: "1brs/D", molecule: "1brs/entity:2", score: 0.3 },
    ],
    alignments: [
      { key: "msa", gapped: [A.slice(0, 40), `${A.slice(0, 20)}--${A.slice(22, 40)}`], flavor: "aligned" },
      {
        key: "fam",
        gapped: ["MKGCAG", "MKGCSG", "MRGCAG"],
        labels: ["P53_HUMAN", "P53_MOUSE", "p53_human"],
        flavor: "aligned",
      },
    ],
    mappings: [{ key: "bind", from: "msa/0", to: "1brs/A", method: "identity" }],
    annotations: [
      {
        key: "domains",
        geometry: "ranges",
        target: "1brs/A",
        semanticType: "uniprot:Domain",
        columns: {
          start: [0, 40, 60, 80],
          stop: [20, 55, 70, 95],
          id: ["PF00069", "PF07714", "PF07714", "PF00017"],
          label: ["Protein kinase", "Tyrosine kinase", "Tyrosine kinase", "SH2"],
          value: [0.95, 0.7, 0.7, 0.85],
        },
      },
      {
        key: "variants",
        geometry: "ranges",
        target: "1brs/A",
        columns: { start: [5, 30, 50], stop: [6, 31, 51], category: ["pathogenic", "benign", "pathogenic"] },
      },
      { key: "sites", geometry: "ranges", target: "1brs/A", columns: { start: [30, 92, 100], stop: [33, 94, 102] } },
      {
        key: "contacts",
        geometry: "pairs",
        target: "1brs/A",
        columns: { start: [10, 20], b_target: ["1brs/D", "1brs/D"], b_start: [5, 9], category: ["hbond", "contact"] },
      },
      { key: "traits", geometry: "rows", target: "msa", columns: { row: [0, 1], identity: [1, 0.85] } },
      {
        key: "hits",
        geometry: "ranges",
        target: "1brs/D",
        columns: { start: [0, 30], stop: [20, 50], evalue: [1e-5, 0.2] },
      },
    ],
  };
}

/** Each example's result: the summary of what it selects. */
const RESULTS: Record<string, string> = {
  "(items :of regions :where (and (> (length) 30) (matches (label) binding)))":
    "1 group, 191 positions: p53/unit 102-292",
  "(positions :in p53 :where (or (= (letter) C) (= (letter) H)))":
    "20 groups, 22 positions: p53/unit 115, 124, 135, 141 (and 16 more)",
  "(items :of regions :where (not (matches (label) binding)))": "2 groups, 76 positions: p53/unit 1-44, 325-356",
  "(positions :in p53 :where (< (index) 10))": "1 group, 10 positions: p53/unit 1-10",
  "(items :of regions :where (<= (length) 40))": "1 group, 32 positions: p53/unit 325-356",
  "(systems :where (> (length) 100))": "1 group, 393 positions: p53/unit 1-393",
  "(pick (letters C) :where (>= (length) 2))":
    "9 groups, 20 positions: gene/unit 8-9, 11-12, 20-21, 32-35 (and 5 more)",
  "(positions :in p53 :where (in (letter) C H))": "20 groups, 22 positions: p53/unit 115, 124, 135, 141 (and 16 more)",
  '(or (sel "1brs/A[10-20]") (sel "1brs/B[30-35]"))': "2 groups, 17 positions: 1brs/A 10-20; 1brs/B 30-35",
  "(and (all) (letters W))": "1 group, 4 positions: p53/unit 23, 53, 91, 146",
  '(sel "1brs/A[10-20]")': "1 group, 11 positions: 1brs/A 10-20",
  '(sel "p53@0:50")': "1 group, 50 positions: p53/unit 1-50",
  "(systems :of 1brs)": "6 groups, 597 positions: 1brs/A 1-110; 1brs/B 1-110; 1brs/C 1-110; 1brs/D 1-89 (and 2 more)",
  "(items :of regions)": "3 groups, 267 positions: p53/unit 1-44, 102-292, 325-356",
  "(values :of cons :above 0.8)": "12 groups, 15 positions: sh2 columns 1, 4, 6, 8 (and 8 more)",
  "(letters CHW :in p53)": "24 groups, 26 positions: p53/unit 23, 53, 91, 115 (and 20 more)",
  "(each (items :of regions) (expand (current) 5))": "3 groups, 292 positions: p53/unit 1-49, 97-297, 320-361",
  "(range 1brs/A 10-20 45)": "1 group, 12 positions: 1brs/A 10-20, 45",
  "(range p53 100-120)": "1 group, 21 positions: p53/unit 100-120",
  "(slice p53 0 50)": "1 group, 50 positions: p53/unit 1-50",
  '(and (expand (sel "1brs/A[10-20]") 5) (not (sel "1brs/A[10-20]")))': "1 group, 10 positions: 1brs/A 5-9, 21-25",
  '(onto (sel "protein[5-10]") gene)': "1 group, 18 positions: gene/unit 25-42",
  "(group (letters C :in p53) :by none)": "1 group, 10 positions: p53/unit 124, 135, 141, 176 (and 6 more)",
  "(and 1brs/A (partners 1brs/D :via contacts))": "1 group, 19 positions: 1brs/A 27, 35, 37-38, 56 (and 6 more)",
  "(each (items :of regions) (and (current) (letters C)))":
    "1 group, 10 positions: p53/unit 124, 135, 141, 176 (and 6 more)",
  "(first (items :of regions) 2)": "2 groups, 235 positions: p53/unit 1-44, 102-292",
  "(positions :in 1brs/A :where (< (index) 10))": "1 group, 10 positions: 1brs/A 1-10",
  "(positions :in p53 :where (= (letter) W))": "4 groups, 4 positions: p53/unit 23, 53, 91, 146",
  "(positions :in model :where (> (value-of plddt) 70))": "2 groups, 198 positions: model/A 47-156; model/B 2-89",
  "(systems :where (matches (ref) ^1brs/))":
    "6 groups, 597 positions: 1brs/A 1-110; 1brs/B 1-110; 1brs/C 1-110; 1brs/D 1-89 (and 2 more)",
  "(pick (all) :where (= (kind) row))":
    "4 groups, 120 positions: sh2/0 151-180; sh2/1 6-36; sh2/2 60-90; sh2/3 333-360",
  "(systems :where (= (origin) sequence))": "1 group, 393 positions: p53/unit 1-393",
  "(systems :where (= (alphabet) dna))": "1 group, 120 positions: gene/unit 1-120",
  "(systems :where (= (field auth_asym_id :ns mvs) A))": "1 group, 110 positions: 1brs/A 1-110",
  "(items :of regions :where (matches (label) binding))": "1 group, 191 positions: p53/unit 102-292",
  "(items :of regions :where (> (length) 50))": "1 group, 191 positions: p53/unit 102-292",
  // The illustrative set's.
  "(items :where (= (type) uniprot:Domain))": "3 groups, 60 positions: 1brs/A 1-20, 40-54, 60-69, 80-94",
  "(items :of variants :where (!= (category) benign))": "2 groups, 2 positions: 1brs/A 6, 50",
  "(rows :of fam :where (matches (label) HUMAN :case true))": "1 group, 6 positions: fam/0 1-6",
  "(and 1brs/A domains)": "1 group, 60 positions: 1brs/A 1-20, 40-54, 60-69, 80-94",
  "(not domains :in 1brs/A)": "1 group, 44 positions: 1brs/A 21-39, 55-59, 70-79, 95-103",
  "(empty)": "selects nothing",
  "(rows :of msa :where (>= (field identity :of traits) 0.9))": "1 group, 40 positions: msa/0 1-40",
  "(and 1brs/A (letters C))": "1 group, 7 positions: 1brs/A 4, 7, 10, 12 (and 3 more)",
  "(search GXGXXG :in 1brs/A :classes protein)": "2 groups, 12 positions: 1brs/A 30A-35, 65-70",
  '(search "[KR]{2}.{10,12}[KR]{3,}" :in p53 :regex true)': "1 group, 17 positions: p53/unit 305-321",
  "(input selection)": "1 group, 20 positions: 1brs/A 11-30",
  "(onto (input selection) msa)": "1 group, 20 positions: msa columns 11-30",
  "(pick (items :of domains) :where (matches (label) kinase))": "2 groups, 45 positions: 1brs/A 1-20, 40-54, 60-69",
  "(ids (items :of domains) PF00069 PF07714)": "2 groups, 45 positions: 1brs/A 1-20, 40-54, 60-69",
  '(positions :in 1brs/A :where (matches (label) "[A-Z]$"))': "1 group, 2 positions: 1brs/A 30A-30B",
  "(rows :of fam :where (matches (label) HUMAN))": "2 groups, 12 positions: fam/0 1-6; fam/2 1-6",
  "(systems :where (>= (score) 0.5))": "1 group, 104 positions: 1brs/A 1-103",
  "(pick (items :of domains) :where (in (id) PF00069 PF07714))": "2 groups, 45 positions: 1brs/A 1-20, 40-54, 60-69",
  "(partners 1brs/D :via contacts :where (= (category) hbond))": "1 group, 1 position: 1brs/A 11",
  "(pick (items :of domains) :where (> (value) 0.8))": "2 groups, 35 positions: 1brs/A 1-20, 80-94",
  "(items :where (= (annotation) variants))": "3 groups, 3 positions: 1brs/A 6, 30A, 50",
  "(items :of hits :where (< (field evalue) 0.001))": "1 group, 20 positions: 1brs/D 1-20",
  "(pick (items :of domains) :where (> (runs) 1))": "1 group, 25 positions: 1brs/A 40-54, 60-69",
  "(pick (search GXGXXG :edits 1) :where (= (score) 0))": "1 group, 6 positions: 1brs/A 13-18",
  '(pick (items :of variants) :where (overlaps (sel "domains[#PF00069]")))': "1 group, 1 position: 1brs/A 6",
  "(pick sites :where (overlaps (pick domains :where (matches (label) SH2))))": "1 group, 2 positions: 1brs/A 92-93",
  '(pick (items :of sites) :where (within (sel "1brs/A[40-90]") 5))': "1 group, 2 positions: 1brs/A 92-93",
  "(pick (items :of domains) :where (>= (count (letters C)) 4))": "1 group, 20 positions: 1brs/A 1-20",
};

const env = createQueryEnv([{ namespace: "app", symbols: [], classes: { protein: { X: "ACDEFGHIKLMNPQRSTVWY" } } }]);
const sets = new Map<string, Promise<SeqSet>>();
const setOf = (fixture: string) => {
  let s = sets.get(fixture);
  if (!s) {
    s =
      fixture === "illustrative"
        ? Promise.resolve(SeqSets.build(illustrative()).set)
        : compile(CASES[fixture](), registries(), fixtureLoaders()).then((r) => r.set);
    sets.set(fixture, s);
  }
  return s;
};

describe("the table's examples, evaluated", () => {
  const examples = BASE_SYMBOLS.flatMap((s) => s.examples.map((e) => ({ symbol: s.name, ...e })));
  for (const ex of examples) {
    it(`${ex.symbol}: ${ex.text} (${ex.fixture ?? "illustrative"})`, async () => {
      const set = await setOf(ex.fixture ?? "illustrative");
      const checked = checkText(ex.text, env.language);
      expect(checked.issues).toEqual([]);
      const selection = SeqSelection.fromExtents(set, {
        count: 1,
        kind: Uint8Array.of(1),
        id: Int32Array.of(0),
        row: Int32Array.of(-1),
        start: Int32Array.of(10),
        stop: Int32Array.of(30),
      });
      const r = compileQuery(checked.expr as never, env).run({ set, inputs: { selection } });
      expect(r.issues.filter((i) => i.severity === "error")).toEqual([]);
      expect([ex.text, SeqSelection.summary(r.selection, { maxParts: 4 })]).toEqual([ex.text, RESULTS[ex.text]]);
    });
  }
});
