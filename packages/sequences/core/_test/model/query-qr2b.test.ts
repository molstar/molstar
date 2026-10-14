/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL's QR2b runtimes (seqspec:seqql#empty-results, seqspec:seqql#performance): `letters`, `positions` and `search`,
// the position, system and group slots' properties, scans in their `and`'s context, and the diagnosis of an empty
// result. Every symbol against the brute-force oracle over random queries, then the shapes and messages the spec names.

import { bindQuery, checkText, defineSymbol, readQuery, renderIssue } from "@molstar/seq-core/query";
import { describe, expect, it } from "vitest";
import { SeqSets } from "../../src/model/index.ts";
import {
  compileQuery,
  createQueryEnv,
  type QueryEnv,
  type QueryFn,
  type QueryOptions,
  SeqSelection as Selections,
  type SeqSelection,
} from "../../src/model/query/index.ts";
import { queryFixture } from "./query-fixtures.ts";
import { groupsOf, Oracle, OracleError } from "./query-oracle.ts";
import { prng } from "./support/generators.ts";

/** A letter-class table as an app registers one: X any amino acid, B and J the ambiguous pairs. */
const CLASSES = { grp: { X: "ACDEFGHIKLMNPQRSTVWY", B: "DN", J: "IL" } };
/** Counts its runs: the group tests' selections are evaluated once per `pick` (seqspec:seqql#elements-and-groups). */
const counted = defineSymbol({
  ...{ name: "test.counted", alias: "counted", kind: "modifier", result: "query", signature: "(counted Q)" },
  ...{ args: [{ name: "Q", type: "query", positional: true, doc: "the query" }], doc: "counts", examples: [] },
  hints: [],
});
const calls = { n: 0 };
const ENV = createQueryEnv([
  {
    namespace: "test",
    symbols: [counted],
    classes: CLASSES,
    runtimes: [
      {
        symbol: counted,
        compile: (c) => {
          const q = c.query("0") as QueryFn<SeqSelection>;
          return (ctx) => {
            calls.n++;
            return q(ctx);
          };
        },
      },
    ],
  },
]);
const input = queryFixture();
input.units = (input.units ?? []).map((u) =>
  u.key === "1brs/A" ? { ...u, origin: "experimental", score: 0.8, custom: { mvs: { auth_asym_id: "A" } } } : u,
);
input.alignments?.push({ key: "ins", gapped: ["ACdeFG-H", "AC-FGkH"], flavor: "a3m" });
input.molecules.push({ key: "far", alphabet: "protein", lazy: { length: 50 }, ownUnit: true });
input.annotations?.push(
  { key: "stats", geometry: "rows", target: "msa", columns: { row: [0, 1, 2], identity: [1, 0.5, 0.9] } },
  // over barstar's molecule: every unit over it reads it
  {
    key: "hyd",
    geometry: "values",
    target: "1brs/entity:2",
    columns: { value: Array.from({ length: 89 }, (_, i) => (i % 10) / 10) },
  },
);
const { set } = SeqSets.build(input);
const SELECTION = Selections.fromExtents(set, {
  count: 1,
  kind: Uint8Array.of(1),
  id: Int32Array.of(Array.from(set.units.key).indexOf(set.strings.find("1brs/A") as number)),
  row: Int32Array.of(-1),
  start: Int32Array.of(40),
  stop: Int32Array.of(60),
});

function run(text: string, o: Partial<QueryOptions> = {}, env: QueryEnv = ENV) {
  const read = readQuery(text);
  const bound = bindQuery(read.node as never, env.language);
  expect([...read.issues, ...bound.issues]).toEqual([]);
  const r = compileQuery(bound.expr as never, env, bound.source).run({ set, inputs: { selection: SELECTION }, ...o });
  return { ...r, said: r.issues.map((i) => renderIssue(i).split("\n")[0]) };
}
const summary = (text: string, o: Partial<QueryOptions> = {}) =>
  Selections.summary(run(text, o).selection, { maxParts: 12 });

describe("QR2b's runtimes against the oracle", () => {
  const SEED = 20261014;
  let rand = prng(SEED);
  const pick = <T>(xs: readonly T[]) => xs[Math.floor(rand() * xs.length)];
  const int = (lo: number, hi: number) => lo + Math.floor(rand() * (hi - lo));
  const LEAVES = [
    "1brs/A",
    "1brs/D",
    "1brs",
    "ab/H",
    "ab",
    "lonely",
    "msa",
    "msa/1",
    "ins/0",
    "(all)",
    "(input selection)",
    "(items :of sites)",
    "dom",
    '(sel "1brs/A[10-40]")',
    '(sel "ab/H[3-12]")',
    '(sel "msa/2@1:5")',
  ];
  const POSITION = [
    "(< (index) 20)",
    "(= (letter) G)",
    "(in (letter) A C D K)",
    "(matches (label) A$)",
    '(matches (label :numbering seq) "^-")',
    "(> (value-of plddt) 70)",
    "(< (value-of sv) 0.8)",
    "(> (value-of cons) 0.8)",
    "(and (> (index) 3) (< (index) 9))",
    "(not (= (letter) Q))",
    "(= (ref) ab/H)",
    "(= (kind) row)",
  ];
  const SYSTEM = [
    "(> (length) 60)",
    "(= (kind) unit)",
    '(matches (ref) "^1brs/[A-C]")',
    "(matches (label) arn)",
    "(= (label) y)",
    "(>= (field identity :of stats) 0.9)",
    "(= (field species :of species) human)",
    "(= (field molecule) ab)",
    "(= (field complex) 1brs)",
    "(= (field row) 1)",
    "(= (field auth_asym_id :ns mvs) A)",
    "(= (alphabet) protein)",
    "(= (origin) experimental)",
    "(> (score) 0.5)",
    "(not (= (alphabet) dna))",
  ];
  const PATTERNS = ["GG", "vq", "AC", "L", "QVQ", "KK", "XG :classes grp", "BJ :classes grp", "G+ :regex true"];
  function gen(depth: number): string {
    if (depth === 0 || rand() < 0.2) return pick(LEAVES);
    const q = () => gen(depth - 1);
    const within = () => (rand() < 0.6 ? ` :in ${q()}` : "");
    const GROUP = [
      "(> (length) 3)",
      "(= (runs) 1)",
      "(> (runs) 1)",
      "(= (score) 0)",
      '(matches (label) "^[ab]$")',
      ...SYSTEM.slice(0, 6),
      `(overlaps ${q()})`,
      `(overlaps ${pick(LEAVES)})`,
      `(within ${q()} ${int(0, 6)})`,
      `(within ${pick(LEAVES)} ${int(0, 30)})`,
      `(within ${pick(LEAVES)} ${int(0, 9)})`,
      `(>= (count ${q()}) ${int(1, 5)})`,
      `(>= (count ${pick(LEAVES)}) ${int(1, 5)})`,
      `(< (count (letters ${pick(["G", "A", "KR"])})) ${int(1, 4)})`,
      `(and (overlaps ${q()}) (> (length) 1))`,
    ];
    const NEAR = GROUP.filter((x) => /\((overlaps|within|count) /.test(x));
    switch (int(0, 12)) {
      case 0:
        return rand() < 0.5
          ? `(letters ${pick(["G", "chw", "QV", "A", "-"])}${within()})`
          : `(and ${q()} (letters ${pick(["G", "KR", "L"])}))`;
      case 1:
        return rand() < 0.6
          ? `(positions${within()} :where ${pick(POSITION)})`
          : `(and ${q()} (positions :where ${pick(POSITION)}))`;
      case 2:
        return rand() < 0.6 ? `(search ${pick(PATTERNS)}${within()})` : `(and ${q()} (search ${pick(PATTERNS)}))`;
      case 3: {
        const kind = pick(["unit", "molecule", "row", "columns"]);
        return rand() < 0.5
          ? `(systems :kind ${kind} :where ${pick(SYSTEM)})`
          : `(rows${rand() < 0.5 ? " :of msa" : ""} :where ${pick([...SYSTEM.slice(3, 7), "(= (field row) 1)"])})`;
      }
      case 4:
      case 5:
      case 6:
      case 11:
        return `(pick ${rand() < 0.3 ? `(search ${pick(PATTERNS)})` : pick([q(), pick(LEAVES), "(items :of sites)"])} :where ${pick(rand() < 0.5 ? NEAR : GROUP)})`;
      case 7:
        return `(expand ${q()} ${int(0, 4)})`;
      case 8:
        return `(range ${q()} 5-30)`;
      default:
        return `(${pick(["and", "or", "or"])} ${q()} ${q()})`;
    }
  }

  const heads = ["letters", "positions", "search", "systems", "rows", "pick", "overlaps", "within", "count"];
  /** `count` random queries from the seed, each against the oracle; how many of each outcome. */
  function againstOracle(count: number) {
    rand = prng(SEED);
    const oracle = new Oracle(set, { selection: groupsOf(SELECTION) });
    oracle.classes = CLASSES;
    const tally = { runs: 0, errors: 0, empty: 0, ...Object.fromEntries(heads.map((h) => [h, 0])) } as Record<
      string,
      number
    >;
    for (let i = 0; i < count; i++) {
      const text = gen(3);
      const checked = checkText(text, ENV.language);
      if (checked.expanded === null) continue;
      const r = run(text);
      const error = r.issues.find((x) => x.severity === "error");
      let want: string;
      try {
        want = JSON.stringify(oracle.eval(checked.expanded));
      } catch (e) {
        if (!(e instanceof OracleError)) throw e;
        want = e.code;
      }
      const got = error ? error.code : JSON.stringify(groupsOf(r.selection));
      expect([text, got]).toEqual([text, want]);
      tally.runs++;
      if (error) tally.errors++;
      else if (r.selection.groupCount === 0) tally.empty++;
      else for (const h of heads) tally[h] += +text.includes(`(${h} `);
    }
    return tally;
  }

  it("gives the oracle's groups, in its order, or its error, for 300 random queries", () => {
    const tally = againstOracle(300);
    expect(tally.runs).toBeGreaterThan(280);
    expect(tally.runs - tally.errors - tally.empty).toBeGreaterThan(120);
    for (const h of heads) expect([h, tally[h] > 2]).toEqual([h, true]);
  });

  it("gives the oracle's groups, in its order, or its error, for 1,500 random queries", {
    tags: ["slow"],
    timeout: 60_000,
  }, () => {
    const tally = againstOracle(1500);
    expect(tally.runs).toBeGreaterThan(1400);
    expect(tally.runs - tally.errors - tally.empty).toBeGreaterThan(600);
    for (const h of heads) expect([h, tally[h] > 40]).toEqual([h, true]);
  });
});

describe("letters, positions and search", () => {
  it("read a unit's letters through its molecule, a row's own residues with its inserts, case-insensitively", () => {
    expect(summary("(letters chw :in 1brs/A)")).toBe("5 groups, 5 positions: 1brs/A 18, 35, 71, 94, 102");
    // ins/0 is ACdeFG-H: d and e are inserts, the row's positions 3 and 4
    expect(summary("(letters DE :in ins/0)")).toBe("1 group, 2 positions: ins/0 3-4");
    expect(summary("(and msa/1 (letters D))")).toBe("1 group, 1 position: msa/1 2");
  });

  it("skip the columns and letters not loaded, with a warning each", () => {
    expect(run("(letters A :in (or msa/0 msa))").said).toEqual([
      "1:1 seqql.letters.columns: alignment columns have no letters and were skipped; name the rows: (and (rows :of A) (letters C))",
    ]);
    expect(run("(letters A :in far)").said).toEqual([
      "1:1 seqql.lazy: the letters of 1 sequence are not loaded and were not read",
    ]);
  });

  it("read only an and's selection when they are its operand with no :in, a molecule without a unit too", () => {
    const limits = { positions: 120 };
    expect(summary("(and 1brs/A (letters W))", { limits })).toBe("1 group, 3 positions: 1brs/A 35, 71, 94");
    expect(run("(letters W)", { limits }).said.at(-1)).toBe(
      "1:1 seqql.limit.positions: this visits 664 positions in 8 sequences and 5 rows (limit 120); narrow first; scan less: (and 1brs/A (letters C)) reads 1brs/A only",
    );
    expect(run("(and (rows :of msa) (letters C))", { limits: { positions: 10 } }).said).toEqual([
      "1:21 seqql.limit.positions: this visits 18 positions in 3 rows (limit 10); narrow first; scan less: (and 1brs/A (letters C)) reads 1brs/A only",
    ]);
    // `(letters G)` alone reads every unit and row; lonely has no unit, but the and names it
    expect(summary("(and lonely (letters G))")).toBe("1 group, 1 position: lonely 6");
    expect(summary("(and ins/1 (search FGK))")).toBe("1 group, 3 positions: ins/1 3-5");
  });

  it("test each position with its index, letter, label in any numbering, and a values annotation's value", () => {
    // ab/H holds indices 2 to 18 of ab, Kabat-labeled with insertion codes; its `seq` numbering runs -5, -4, …
    expect(summary('(positions :in ab/H :where (matches (label) "[A-Z]$"))')).toBe("1 group, 2 positions: ab/H 3A-3B");
    expect(summary('(positions :in ab/H :where (matches (label :numbering seq) "^-"))')).toBe(
      "1 group, 3 positions: ab/H 3-3B",
    );
    expect(summary("(positions :in 1brs/E :where (< (index) 3))")).toBe("1 group, 3 positions: 1brs/E 1-3");
    // hyd is over barstar's molecule: a unit reads it; plddt is chain A's own
    expect(summary("(positions :in 1brs/E :where (>= (value-of hyd) 0.8))")).toMatch(
      /^9 groups, 17 positions: 1brs\/E 9-10, 19-20/,
    );
    expect(summary("(positions :in 1brs/D :where (> (value-of sv) 0.8))")).toBe(
      "2 groups, 13 positions: 1brs/D 1-5, 13-20",
    );
    expect(summary("(positions :in msa :where (> (value-of cons) 0.85))")).toBe(
      "3 groups, 3 positions: msa columns 1, 4, 8",
    );
    expect(run("(positions :in 1brs/B :where (> (value-of dom) 1))").said).toEqual([
      '1:43 seqql.ref.kind: "dom" is a ranges annotation; value-of takes a values annotation (cons, plddt)',
    ]);
  });

  it("search with each matcher: exact, a motif of classes, a regular expression, edits with their ends", () => {
    // ab/H holds indices 2 to 18: QVQ at 0 is ab/L's only
    expect(summary("(search QVQ :in ab)")).toBe("1 group, 3 positions: ab/L 1-3");
    expect(summary("(search XG :classes grp :in ab/H)")).toBe("5 groups, 7 positions: ab/H @6:10, 10-12");
    expect(summary('(search "G+" :regex true :in ab)')).toBe(
      "10 groups, 10 positions: ab/H @7:10, 11-12; ab/L @7:10, 11-12",
    );
    expect(summary("(search QVQLVQ :edits 1 :in ab)")).toBe("1 group, 6 positions: ab/L 1-4");
    const all = run("(search QVQLVQ :edits 1 :ends all :in ab)").selection;
    expect([all.groupCount, Array.from(all.groupScore)]).toEqual([2, [1, 1]]);
    expect(summary("(pick (search QVQLVQ :edits 1 :ends all :in ab) :where (= (score) 1))")).toBe(
      "2 groups, 6 positions: ab/L 1-4",
    );
  });

  it("search a cyclic molecule read whole as a circle: a hit across its origin is one group of two runs", () => {
    const { set: ring } = SeqSets.build({
      molecules: [{ key: "pla", residues: "TTCCAAAAGAA", alphabet: "dna", cyclic: true, ownUnit: true }],
    });
    const env = createQueryEnv();
    const s = compileQuery(checkText("(search AATT)").expr as never, env).run({ set: ring }).selection;
    expect(Selections.summary(s, { indices: true })).toBe("1 group, 4 positions: pla/unit 1-2 @0:2, 10-11 @9:11");
    const part = compileQuery(checkText('(search AATT :in (sel "pla@0:11"))').expr as never, env).run({ set: ring });
    expect(part.selection.groupCount).toBe(1);
    expect(
      compileQuery(checkText('(search AATT :in (sel "pla@1:11"))').expr as never, env).run({ set: ring }).selection
        .groupCount,
    ).toBe(0);
  });

  it("say what is wrong with a pattern or a class table", () => {
    expect(run("(search QV :classes nope)").said).toEqual([
      '1:21 seqql.classes.unknown: no letter-class table "nope" here; its tables: grp',
    ]);
    expect(run('(search "(" :regex true)').said).toEqual([
      '1:9 seqql.type: search cannot use the pattern "(": Invalid regular expression: /(/i: Unterminated group',
    ]);
    expect(run("(search QV :edits 2)").said).toEqual([
      '1:9 seqql.type: search cannot use the pattern "QV": 2 edits in a pattern of 2 letters would match anything; use at most 1',
    ]);
    expect(run(`(search ${"X".repeat(65)} :classes grp)`).said[0]).toMatch(/65 positions; the limit is 64/);
  });
});

describe("the system and group slots", () => {
  it("read a system's ref, kind, label, length, origin, score, alphabet and fields", () => {
    expect(summary("(rows :of msa :where (>= (field identity :of stats) 0.9))")).toBe(
      "2 groups, 12 positions: msa/0 11-16; msa/2 1-6",
    );
    expect(summary("(rows :where (= (field species :of species) human))")).toBe("1 group, 6 positions: msa/0 11-16");
    expect(summary("(systems :where (= (field auth_asym_id :ns mvs) A))")).toBe("1 group, 110 positions: 1brs/A 1-110");
    expect(summary("(systems :where (and (= (origin) experimental) (> (score) 0.5)))")).toBe(
      "1 group, 110 positions: 1brs/A 1-110",
    );
    expect(summary("(systems :kind row :where (matches (label) y))")).toBe("1 group, 6 positions: msa/1 1-6");
    expect(summary("(systems :kind columns :where (= (ref) msa))")).toBe("1 group, 8 positions: msa columns 1-8");
    expect(summary("(rows :where (= (field alignment) ins))")).toBe("2 groups, 13 positions: ins/0 1-7; ins/1 1-6");
    expect(summary("(systems :where (> (length) 100))")).toBe(
      "3 groups, 330 positions: 1brs/A 1-110; 1brs/B 1-110; 1brs/C 1-110",
    );
    expect(run("(rows :of msa :where (= (field nope :of species) 1))").said).toEqual([
      '1:25 seqql.field.unknown: species has no field "nope"; its fields: species; a system\'s own fields: molecule, complex, alignment, row; others need :of or :ns',
    ]);
  });

  it("test a group by its length, runs and score, and by a selection evaluated once per pick", () => {
    // sites: 1brs/A [5, 10), [30, 40), [60, 70), [8, 12)
    expect(summary('(pick (items :of sites) :where (overlaps (sel "1brs/A[30-35]")))')).toBe(
      "1 group, 10 positions: 1brs/A 31-40",
    );
    // [30, 40) ends 5 before index 44, [60, 70) starts 11 after index 49
    expect(summary('(pick (items :of sites) :where (within (sel "1brs/A[45-50]") 10))')).toBe(
      "1 group, 10 positions: 1brs/A 31-40",
    );
    expect(summary('(pick (items :of sites) :where (within (sel "1brs/A[45-50]") 11))')).toBe(
      "2 groups, 20 positions: 1brs/A 31-40, 61-70",
    );
    expect(summary("(pick (items :of sites) :where (>= (count (letters G :in 1brs/A)) 2))")).toBe(
      "2 groups, 20 positions: 1brs/A 31-40, 61-70",
    );
    calls.n = 0;
    expect(summary("(pick (items :of sites) :where (overlaps (counted (letters G :in 1brs/A))))")).toMatch(/^4 groups/);
    expect(calls.n).toBe(1);
    // once per evaluation of the pick: per group of each (dom's two items)
    calls.n = 0;
    expect(summary("(each (items :of dom) (pick (items :of sites) :where (overlaps (counted (current)))))")).toBe(
      "1 group, 10 positions: 1brs/A 31-40",
    );
    expect(calls.n).toBe(2);
  });
});

describe("an empty result says why", () => {
  it("names the first node that gave nothing from something, with what it knows", () => {
    const why = (text: string) => run(text).said.find((x) => x.includes("seqql.empty"));
    expect(why("(values :of plddt :above 99)")).toBe(
      "1:1 seqql.empty: (values :of plddt :above 99) selected nothing: plddt runs 5 to 95",
    );
    expect(why("(and 1brs/A (letters c))")).toBe(
      "1:13 seqql.empty: (letters c) selected nothing: no C in 110 positions",
    );
    expect(why("(search WWW :in 1brs/A)")).toBe(
      "1:1 seqql.empty: (search WWW :in 1brs/A) selected nothing: no hit in 110 positions; allow edits with :edits 1",
    );
    expect(why("(positions :in 1brs/B :where (> (value-of plddt) 1))")).toBe(
      "1:1 seqql.empty: (positions :in 1brs/B :where (> (value-of plddt) 1)) selected nothing: none of 110 positions passed :where; (value-of plddt) was missing 110 times",
    );
    expect(why("(pick (items :of sites) :where (> (value) 1))")).toBe(
      "1:1 seqql.empty: (pick (items :of sites) :where (> (value) 1)) selected nothing: none of 4 groups passed :where; (value) was missing 4 times",
    );
    expect(why("(rows :of msa :where (= (field species :of species) cat))")).toBe(
      "1:1 seqql.empty: (rows :of msa :where (= (field species :of species) cat)) selected nothing: none of 3 rows passed :where; (field species :of species) was missing once",
    );
    expect(why("(items :of sites :where (= (label) zz))")).toMatch(/none of 4 items passed :where$/);
    expect(why("(partners 1brs/B :via contacts)")).toMatch(/selected nothing: no pair of contacts touches it$/);
    expect(why("(partners 1brs/E :via contacts :where (= (category) x))")).toMatch(
      /1 pair end of contacts touches it; none passed :where or :side$/,
    );
    expect(why("(and 1brs/A (not 1brs))")).toBe(
      "1:1 seqql.empty: (and 1brs/A (not 1brs)) selected nothing: 1brs covers all of what it is taken from (1 group, 110 positions: 1brs/A 1-110)",
    );
    // a modifier without more to say: what it had
    expect(why("(range (letters W :in 1brs/A) 40-60)")).toBe(
      "1:1 seqql.empty: (range (letters W :in 1brs/A) 40-60) selected nothing from 3 groups, 3 positions: 1brs/A 35 @34:35, 71 @70:71, 94 @93:94",
    );
  });

  it("stays quiet when the result is not empty, or a warning of the node named already says why", () => {
    expect(run("(or (values :of plddt :above 99) 1brs/A)").said).toEqual([]);
    expect(run("(range 1brs/A 300)").said.map((x) => x.split(" ")[1])).toEqual(["seqql.label.unknown:"]);
    expect(run("(onto msa/1 1brs/A)").said.map((x) => x.split(" ")[1])).toEqual(["seqql.onto.unmapped:"]);
    // an input the host passed empty is named as such
    expect(run("(and (input none) 1brs/A)", { inputs: { none: Selections.empty(set) } }).said).toEqual([
      "1:6 seqql.empty: (input none) selected nothing: the app passed it empty",
    ]);
  });
});
