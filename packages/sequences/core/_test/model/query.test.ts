/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL's runtime core (seqspec:seqql#the-selection, seqspec:seqql#extension, seqspec:seqql#errors, QR1): every symbol
// against a brute-force oracle over random queries, the shapes the spec names (labels with insertion codes and gaps,
// slices at the edges, molecule targets on units), the runtime's issues, the limits, `toLoci` against the model's loci,
// and extension runtimes.

import {
  BASE_SYMBOLS,
  bindQuery,
  checkText,
  defineSymbol,
  type QueryIssue,
  readQuery,
  renderIssue,
} from "@molstar/seq-core/query";
import { describe, expect, it, vi } from "vitest";
import { Extents, Loci, OrderedSet, type SeqIndex, type SeqLoci, SeqSets, type UnitId } from "../../src/model/index.ts";
import {
  compileQuery,
  createQueryEnv,
  type QueryEnv,
  type QueryOptions,
  SelectionBuilder,
  SeqSelection,
} from "../../src/model/query/index.ts";
import { queryFixture, REFS } from "./query-fixtures.ts";
import { groupsOf, Oracle, OracleError } from "./query-oracle.ts";
import { prng } from "./support/generators.ts";

const { set } = SeqSets.build(queryFixture());
const ENV = createQueryEnv();
const unit = (key: string) => set.strings.find(key) && Array.from(set.units.key).indexOf(set.strings.find(key));
const seqLoci = (elements: [string, number, number][]): SeqLoci => ({
  kind: "seq-loci",
  set,
  elements: elements.map(([k, s, t]) => ({
    unit: unit(k) as UnitId,
    indices: OrderedSet.ofRange(s as SeqIndex, t as SeqIndex),
  })),
});
const SELECTION = SeqSelection.fromLoci(set, [seqLoci([["1brs/A", 40, 45]]), seqLoci([["ab/H", 5, 9]])]);

/** Reads, binds and runs a text; the issues rendered (`line:column code: message`). */
function run(text: string, o: Partial<QueryOptions> = {}, env: QueryEnv = ENV) {
  const read = readQuery(text);
  const bound = bindQuery(read.node as never, env.language);
  expect([...read.issues, ...bound.issues]).toEqual([]);
  const r = compileQuery(bound.expr as never, env, bound.source).run({ set, inputs: { selection: SELECTION }, ...o });
  return { ...r, said: r.issues.map((i) => renderIssue(i).split("\n")[0]) };
}
const summary = (text: string) => SeqSelection.summary(run(text).selection, { indices: true });
const codes = (issues: readonly QueryIssue[]) => issues.map((i) => `${i.severity} ${i.code}`);

describe("the runtime against the oracle", () => {
  const SEED = 20261012;
  let rand = prng(SEED);
  const pick = <T>(xs: readonly T[]) => xs[Math.floor(rand() * xs.length)];
  const int = (lo: number, hi: number) => lo + Math.floor(rand() * (hi - lo));
  const LABELS = [
    "1",
    "3",
    "3A",
    "3B",
    "4",
    "5",
    "6",
    "7",
    "9",
    "9.1",
    "9.2",
    "10",
    "13",
    "16",
    "-3",
    "-1",
    "2",
    "45",
    "108",
    "111",
  ];
  const label = () => (rand() < 0.5 ? pick(LABELS) : String(int(-4, 115)));
  const part = () => (rand() < 0.6 ? `${label()}-${label()}` : label());
  const NUMBERINGS = [
    ...Array(8).fill(""),
    " :numbering kabat",
    " :numbering seq",
    " :numbering auth",
    " :numbering label",
    " :numbering row",
    " :numbering nope",
  ];
  const LEAVES = [
    ...REFS,
    "(all)",
    "(empty)",
    "(input selection)",
    "(items)",
    "(items :of sites)",
    "(items :of multi)",
    "(items :of contacts)",
    "(values :of cons :above 0.8)",
    "(values :of plddt :min 60 :below 90)",
    "(values :of sv :above 0.5)",
    "(values :of sites :above 0.5)",
    "(items :of cons)",
    "(systems)",
    "(systems :of 1brs)",
    "(systems :kind molecule)",
    "(systems :of ab :kind unit)",
    "(systems :kind columns)",
    "(rows :of msa)",
    "(rows)",
  ];
  function gen(depth: number): string {
    if (depth === 0 || rand() < 0.25) return pick(LEAVES);
    const q = () => gen(depth - 1);
    switch (int(0, 7)) {
      case 0:
        return `(range ${q()} ${Array.from({ length: int(1, 4) }, part).join(" ")}${pick(NUMBERINGS)})`;
      case 1: {
        const pairs = Array.from({ length: int(1, 3) }, () => {
          const a = rand() < 0.5 ? pick([0, 1, 2, 3, 16, 17, 18, 19, 20, 88, 89, 109, 110, 111]) : int(0, 112);
          return `${a} ${a + (rand() < 0.5 ? 1 : int(1, 40))}`;
        });
        return `(slice ${q()} ${pairs.join(" ")})`;
      }
      case 2:
        return `(and ${q()} ${rand() < 0.4 ? `(not ${q()})` : q()}${rand() < 0.3 ? ` ${q()}` : ""})`;
      case 3:
        return `(or ${q()} ${q()}${rand() < 0.3 ? ` ${q()}` : ""})`;
      case 4:
        return rand() < 0.5 ? `(not ${q()})` : `(not ${q()} :in ${q()})`;
      default:
        return `(and ${pick(["1brs/A", "ab", "msa", "1brs", "(all)"])} ${q()})`;
    }
  }

  /** `count` random queries from the seed, each against the oracle; how many of each outcome. */
  function againstOracle(count: number) {
    rand = prng(SEED);
    const oracle = new Oracle(set, { selection: groupsOf(SELECTION) });
    const tally = { runs: 0, errors: 0, empty: 0, range: 0, slice: 0, warned: 0 };
    for (let i = 0; i < count; i++) {
      const text = gen(3);
      const checked = checkText(text);
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
      else {
        tally.range += +text.includes("(range");
        tally.slice += +text.includes("(slice");
      }
      tally.warned += +r.issues.some((x) => x.severity === "warning");
    }
    return tally;
  }

  it("gives the oracle's groups, in its order, or its error, for 300 random queries", () => {
    const tally = againstOracle(300);
    // Enough of each outcome that the comparison means something.
    expect(tally.runs).toBeGreaterThan(280);
    expect(tally.errors).toBeGreaterThan(6);
    expect(tally.runs - tally.errors - tally.empty).toBeGreaterThan(100);
    expect([tally.range > 16, tally.slice > 16, tally.warned > 30]).toEqual([true, true, true]);
  });

  it("gives the oracle's groups, in its order, or its error, for 1,500 random queries", {
    tags: ["slow"],
    timeout: 60_000,
  }, () => {
    const tally = againstOracle(1500);
    // Enough of each outcome that the comparison means something.
    expect(tally.runs).toBeGreaterThan(1400);
    expect(tally.errors).toBeGreaterThan(30);
    expect(tally.runs - tally.errors - tally.empty).toBeGreaterThan(500);
    // Seeded: 1,500 queries give about 700 non-empty results, 100 each through range and slice, 400 errors.
    expect([tally.range > 80, tally.slice > 80, tally.warned > 150]).toEqual([true, true, true]);
    // A second under a normal load; a shared machine's load can make it many.
  });
});

describe("generators and refs", () => {
  it("name systems by ref, one group each: a complex its units, a molecule its alive units, else itself", () => {
    expect(summary("1brs/A")).toBe("1 group, 110 positions: 1brs/A 1-110 @0:110");
    expect(summary("ab")).toBe("2 groups, 36 positions: ab/H 3-14 @2:18; ab/L 1-16 @0:20");
    expect(summary("lonely")).toBe("1 group, 10 positions: lonely 1-10 @0:10");
    expect(summary("msa/1")).toBe("1 group, 6 positions: msa/1 1-6 @0:6");
    expect(summary("msa")).toBe("1 group, 8 positions: msa columns 1-8 @0:8");
    expect(run("1brs").selection.groupCount).toBe(6);
    expect(run("(all)").selection.groupCount).toBe(6 + 2 + 3);
    expect(run("(empty)").selection.groupCount).toBe(0);
  });

  it("give items with their extents, tagged: a molecule's target on every unit over it, a pair's both ends", () => {
    const dom = run("dom").selection;
    expect([dom.groupCount, Array.from(dom.groupItem)]).toEqual([2, [0, 1]]);
    expect(SeqSelection.summary(dom, { maxParts: 3 })).toBe(
      "2 groups, 135 positions: 1brs/A 21-50, 91-105; 1brs/B 21-50 (and 3 more)",
    );
    expect(summary("(items :of multi)")).toBe("2 groups, 9 positions: 1brs/D 1-4 @0:4; ab/H 3 @2:3; ab/L 12-15 @15:19");
    expect(summary("(items :of contacts)")).toBe(
      "3 groups, 95 positions: 1brs/A 34 @33:34, 59 @58:59, 102 @101:102; 1brs/D 29-30 @28:30, 38 @37:38; 1brs/E 1-89 @0:89",
    );
    expect(summary("species")).toBe("2 groups, 12 positions: msa/0 11-16 @0:6; msa/2 1-6 @0:6");
    expect(summary("lone")).toBe("1 group, 4 positions: lonely 3-6 @2:6");
  });

  it("select by value, one group per run; NaN never passes", () => {
    expect(summary("(values :of cons :above 0.8)")).toBe(
      "3 groups, 5 positions: msa columns 1 @0:1, 3-4 @2:4, 7-8 @6:8",
    );
    expect(summary("(values :of cons :min 0.5 :max 0.9)")).toBe(
      "3 groups, 4 positions: msa columns 1 @0:1, 3 @2:3, 6-7 @5:7",
    );
    expect(summary("(values :of sv :above 0.8)")).toBe("2 groups, 13 positions: 1brs/D 1-5 @0:5, 13-20 @12:20");
  });

  it("give systems by kind within :of", () => {
    expect(summary("(systems :of 1brs :kind molecule)")).toBe(
      "2 groups, 199 positions: 1brs/entity:1 1-110 @0:110; 1brs/entity:2 1-89 @0:89",
    );
    expect(summary("(systems :of 1brs/A :kind molecule)")).toBe("1 group, 110 positions: 1brs/entity:1 1-110 @0:110");
    expect(summary("(rows :of msa/2)")).toBe("1 group, 6 positions: msa/2 1-6 @0:6");
    expect(summary("(systems :kind columns)")).toBe("1 group, 8 positions: msa columns 1-8 @0:8");
  });

  it("give a complex's molecules in canonical order, whatever its units' order (CM-12)", () => {
    const { set: s } = SeqSets.build({
      molecules: [
        { key: "m0", residues: "MKV" },
        { key: "m1", residues: "ACDEF" },
      ],
      units: [
        { key: "u0", molecule: "m1" },
        { key: "u1", molecule: "m0" },
      ],
      complexes: [{ key: "c", units: ["u0", "u1"] }],
    } as never);
    const bound = bindQuery(readQuery("(systems :of c :kind molecule)").node as never, ENV.language);
    const r = compileQuery(bound.expr as never, ENV).run({ set: s });
    // The groups themselves in order (a summary sorts what it prints).
    expect(Array.from(r.selection.groupStart.subarray(0, 2), (k) => r.selection.id[k])).toEqual([0, 1]);
  });
});

describe("range and slice", () => {
  it("read labels in the system's default numbering, with insertion codes, decimals and gaps", () => {
    expect(summary('(sel "ab/L[3A-4]")')).toBe("1 group, 3 positions: ab/L 3A-4 @3:6");
    expect(summary('(sel "ab/L[9.1-10]")')).toBe("1 group, 3 positions: ab/L 9.1-10 @11:14");
    // absent endpoints placed by comparison: 5 lands on 7 (after the gap), 13 on 14
    expect(summary('(sel "ab/L[5-8, 13-14]")')).toBe("1 group, 3 positions: ab/L 7-8 @8:10, 14 @17:18");
    expect(summary('(sel "ab/L[3B-9.1]")')).toBe("1 group, 8 positions: ab/L 3B-9.1 @4:12");
    // a ranged unit: its labels are its molecule's, its indices too
    expect(summary('(sel "ab/H[1-3A]")')).toBe("1 group, 2 positions: ab/H 3-3A @2:4");
    expect(summary("(range ab/L -1-1 :numbering seq)")).toBe("1 group, 2 positions: ab/L 3B-4 @4:6");
    expect(summary("(range msa 2-3 7)")).toBe("1 group, 3 positions: msa columns 2-3 @1:3, 7 @6:7");
    // a run from an unlabeled position is shown by its slice
    expect(SeqSelection.summary(run('(sel "ab/L@6:9")').selection)).toBe("1 group, 3 positions: ab/L @6:9");
    expect(summary("(range msa/0 12-13)")).toBe("1 group, 2 positions: msa/0 12-13 @1:3");
  });

  it("warn of a label range covering nothing, and fail on a numbering no system has", () => {
    const r = run('(sel "1brs/A[300-310]")');
    expect([r.selection.groupCount, r.said]).toEqual([
      0,
      [
        "1:1 seqql.label.unknown: 1brs/A has no labels 300-310; its labels run 1-110; labels are as the ruler shows them; for indices from 0 use (slice Q start stop)",
      ],
    ]);
    expect(run("(range ab 6)").said[0]).toMatch(/ab\/L has no labels 6; its labels run 1-16/);
    expect(run("(range 1brs/A 10-20 :numbering uniprot)").said).toEqual([
      '1:32 seqql.numbering.unknown: 1brs/A has no numbering "uniprot"; it has auth, ordinal, label',
    ]);
    // systems without it are left out when another has it
    expect(summary("(range (or 1brs/A ab/L) 3A :numbering kabat)")).toBe("1 group, 1 position: ab/L 3A @3:4");
  });

  it("clip slices to each domain at its edges: past it a warning, wholly outside an error", () => {
    expect(summary('(sel "ab/H@2:3")')).toBe("1 group, 1 position: ab/H 3 @2:3");
    expect(summary('(sel "ab/H@17:18")')).toBe("1 group, 1 position: ab/H 14 @17:18");
    const clipped = run('(sel "ab/H@0:3")');
    expect([SeqSelection.summary(clipped.selection, { indices: true }), clipped.said]).toEqual([
      "1 group, 1 position: ab/H 3 @2:3",
      [
        "1:1 seqql.slice.clipped: (slice ab/H 0 3) reaches past it, and was clipped: ab/H holds indices 2 to 18 of its molecule (labels 3-14); indices run from 0 to the length; the summary gives each run's slice",
      ],
    ]);
    expect(run('(sel "ab/H@18:19")').said[0]).toMatch(/^1:1 seqql.slice.range: \(slice ab\/H 18 19\) is outside it/);
    expect(run("(slice 1brs/A 0 1 110 111)").said[0]).toMatch(
      /\(slice 1brs\/A 110 111\) is outside it: 1brs\/A has 110/,
    );
    expect(summary("(slice (or 1brs/A 1brs/D) 100 120)")).toBe("1 group, 10 positions: 1brs/A 101-110 @100:110");
    expect(summary("(slice 1brs/A 30 40 9 20 15 25)")).toBe("1 group, 26 positions: 1brs/A 10-25 @9:25, 31-40 @30:40");
    expect(run("(slice (empty) 5 6)").issues).toEqual([]);
  });

  it("drop the item a group came from when they change its runs", () => {
    const r = run("(range sites 1-200)").selection;
    expect([r.groupCount, Array.from(r.groupAnnotation)]).toEqual([4, [-1, -1, -1, -1]]);
  });
});

describe("combinators", () => {
  it("intersect each group of the first with each other operand in turn; remove a (not B) without its complement", () => {
    expect(summary("(and 1brs dom)")).toBe(
      "3 groups, 135 positions: 1brs/A 21-50 @20:50, 91-105 @90:105; 1brs/B 21-50 @20:50, 91-105 @90:105; 1brs/C 21-50 @20:50, 91-105 @90:105",
    );
    expect(summary("(and 1brs/A dom sites)")).toBe("1 group, 10 positions: 1brs/A 31-40 @30:40");
    expect(summary("(and 1brs/A (not dom))")).toBe(
      "1 group, 65 positions: 1brs/A 1-20 @0:20, 51-90 @50:90, 106-110 @105:110",
    );
    // a molecule with no unit is outside (all), so its complement would be empty: removal keeps it
    expect(summary("(and lonely (not lone))")).toBe("1 group, 6 positions: lonely 1-2 @0:2, 7-10 @6:10");
    expect(summary("(and lonely (not lone :in (all)))")).toBe("selects nothing");
  });

  it("merge in operand order, dropping a group equal to an earlier one (runs and item)", () => {
    const r = run('(or sites 1brs/A sites (items :of sites) (sel "1brs/A[1-110]"))').selection;
    expect([r.groupCount, Array.from(r.groupItem)]).toEqual([5, [0, 1, 2, 3, -1]]);
  });

  it("complement within :in, one group", () => {
    expect(summary("(not sites :in 1brs/A)")).toBe(
      "1 group, 83 positions: 1brs/A 1-5 @0:5, 13-30 @12:30, 41-60 @40:60, 71-110 @70:110",
    );
    expect(run("(not (all))").selection.groupCount).toBe(0);
  });
});

describe("the runtime's issues", () => {
  it("names unknown refs and inputs with did-you-mean, and refs of the wrong kind with the symbol's hint", () => {
    const said = (t: string) => run(t).said[0];
    expect(said("1brs/Q")).toBe('1:1 seqql.ref.unknown: unknown ref "1brs/Q"; did you mean "1brs/A"?');
    expect(said("(input focus)")).toBe(
      '1:8 seqql.input.unknown: no input "focus" here; inputs given: selection; inputs are named by the app: (input selection)',
    );
    expect(said("(values :of sites)")).toBe(
      '1:13 seqql.ref.kind: "sites" is a ranges annotation; :of takes a values annotation (cons, plddt); for items by value: (pick (items :of A) :where (> (value) 0.8))',
    );
    expect(said("(items :of cons)")).toMatch(
      /^1:12 seqql.ref.kind: "cons" is a values annotation; values and matrices/,
    );
    expect(said("(items :of m)")).toMatch(/"m" is a matrix annotation/);
    expect(said("(rows :of 1brs/A)")).toBe(
      '1:11 seqql.ref.kind: "1brs/A" is a unit; :of takes an alignment (msa); a row alone is msa/3',
    );
    expect(said("(systems :of msa)")).toBe(
      '1:14 seqql.ref.kind: "msa" is an alignment; systems of kind unit take :of a molecule or a complex',
    );
    // an unknown ref where a kind is wanted: the refs of that kind, when none is near
    expect(said("(values :of stats)")).toBe(
      '1:13 seqql.ref.unknown: unknown ref "stats"; values annotations: cons, plddt, sv',
    );
    expect(said("(rows :of mas)")).toBe('1:11 seqql.ref.unknown: unknown ref "mas"; did you mean "msa"?');
    // and over sequences that share none: a warning that says to carry one onto the other
    expect(run("(and 1brs/A msa)").said).toEqual([
      "1:1 seqql.empty: and meets nothing: its operands share no sequence (1brs/A, msa columns); carry one onto the other's first: (onto Q R)",
    ]);
    // over the same sequence: the empty-result diagnosis (QR2b) names the operand that met nothing
    expect(run('(and 1brs/A dom (sel "1brs/A[60-70]"))').said).toEqual([
      "1:1 seqql.empty: (and 1brs/A dom (range 1brs/A 60-70)) selected nothing: (range 1brs/A 60-70) meets none of what came before it (1 group, 45 positions: 1brs/A 21-50, 91-105)",
    ]);
    const issue = run("(and 1brs/A\n  (values :of sites))").issues[0];
    expect([issue.at, issue.signature]).toEqual([
      { line: 2, column: 15 },
      "(values :of A :min x :max x :above x :below x)",
    ]);
  });

  it("emits every runtime code a hint of a symbol it runs carries (skill.test.ts' RUNTIME)", () => {
    const EMITS: Record<string, string> = {
      "seq.generator.rows seqql.ref.kind": "(rows :of 1brs/A)",
      "seq.generator.items seqql.ref.kind": "(items :of cons)",
      "seq.generator.values seqql.ref.kind": "(values :of sites)",
      "seq.generator.input seqql.input.unknown": "(input focus)",
      "seq.modifier.range seqql.label.unknown": "(range 1brs/A 300)",
      "seq.modifier.slice seqql.slice.range": "(slice 1brs/A 200 300)",
      "seq.modifier.slice seqql.slice.clipped": "(slice 1brs/A 0 300)",
      "seq.modifier.onto seqql.onto.unmapped": "(onto msa/1 1brs/A)",
      "seq.modifier.partners seqql.ref.kind": "(partners 1brs/A :via sites)",
      "seq.item.field seqql.field.unknown": "(items :of sites :where (= (field nope) 1))",
      "seq.generator.positions seqql.limit.positions": "(positions :in (all) :where (> (index) 5))",
      "seq.generator.letters seqql.limit.positions": "(letters C)",
      "seq.generator.search seqql.limit.positions": "(search GG)",
      "seq.position.value_of seqql.ref.kind": "(positions :in 1brs/A :where (> (value-of sites) 1))",
      "seq.system.field seqql.field.unknown": "(rows :of msa :where (= (field nope :of species) 1))",
    };
    const runtime = /^seqql\.(ref|input|label|slice|onto|limit|field|numbering|lazy|empty|classes|version)\./;
    const reached = BASE_SYMBOLS.flatMap((s) =>
      ENV.runtimes.has(s.name) ? s.hints.filter((h) => runtime.test(h.code)).map((h) => `${s.name} ${h.code}`) : [],
    );
    expect(reached.sort()).toEqual(Object.keys(EMITS).sort());
    for (const [k, text] of Object.entries(EMITS)) {
      const [name, code] = k.split(" ");
      const limits = code === "seqql.limit.positions" ? { positions: 100 } : {};
      const issue = run(text, { limits }).issues.find((i) => i.code === code);
      expect([k, issue?.signature]).toEqual([k, BASE_SYMBOLS.find((s) => s.name === name)?.signature]);
    }
  });

  it("fails to compile a symbol without a runtime, at its place", () => {
    const orphan = defineSymbol({
      ...{ name: "test.orphan", alias: "orphan", kind: "generator", args: [], result: "query", signature: "(orphan)" },
      ...{ doc: "a symbol with no runtime", examples: [], hints: [] },
    });
    const r = run("(and 1brs/A (orphan))", {}, createQueryEnv([{ namespace: "test", symbols: [orphan] }]));
    expect(r.said).toEqual([
      "1:13 seqql.symbol.unavailable: orphan cannot run here: this host has no runtime for test.orphan",
    ]);
    expect(r.selection.groupCount).toBe(0);
  });

  it("stops at the limits: positions scanned, runs and groups built, time, an abort", () => {
    const limit = (text: string, limits: QueryOptions["limits"]) => codes(run(text, { limits }).issues);
    expect(limit("(values :of plddt :above 60)", { positions: 100 })).toEqual(["error seqql.limit.positions"]);
    expect(run("(values :of plddt :above 60)", { limits: { positions: 100 } }).said[0]).toMatch(
      /this visits 110 positions in plddt \(limit 100\); narrow first/,
    );
    expect(limit("(all)", { groups: 10 })).toEqual(["error seqql.limit.groups"]);
    expect(limit("(values :of plddt :above 40)", { groups: 1e6, runs: 1 })).toEqual([]);
    const big = SeqSets.build({
      molecules: [{ key: "big", residues: "A".repeat(5000), ownUnit: true }],
      annotations: [
        {
          key: "v",
          geometry: "values",
          target: "big",
          columns: { value: Float64Array.from({ length: 5000 }, (_, i) => i % 2) },
        },
      ],
    }).set;
    const bigRun = (o: Partial<QueryOptions>) => {
      const q = compileQuery(checkText("(values :of v :above 0.5)").expr as never, ENV);
      return codes(q.run({ set: big, ...o }).issues);
    };
    expect(bigRun({ limits: { runs: 1500 } })).toEqual(["error seqql.limit.runs"]);
    const now = vi.spyOn(Date, "now");
    now.mockReturnValueOnce(0).mockReturnValue(10_000);
    expect(bigRun({ limits: { timeMs: 50 } })).toEqual(["error seqql.limit.time"]);
    now.mockRestore();
    const stop = new AbortController();
    stop.abort();
    expect(bigRun({ signal: stop.signal })).toEqual(["error seqql.aborted"]);
  });
});

describe("the selection", () => {
  it("builds canonical groups: touching runs merged, an out-of-order group sorted, empty groups dropped", () => {
    const b = new SelectionBuilder();
    b.begin();
    b.add(1, 3, -1, 10, 20);
    b.add(1, 3, -1, 20, 25);
    b.add(1, 0, -1, 5, 6);
    b.add(1, 3, -1, 0, 2);
    b.add(1, 3, -1, 1, 4);
    b.end();
    b.begin();
    b.add(1, 0, -1, 4, 4);
    b.end();
    const s = b.build(set);
    expect([s.groupCount, s.runCount, Array.from(s.id), Array.from(s.start), Array.from(s.stop)]).toEqual([
      1,
      3,
      [0, 3, 3],
      [5, 0, 10],
      [6, 4, 25],
    ]);
    expect([SeqSelection.positionCount(s), SeqSelection.positionCount(s, 0)]).toEqual([20, 20]);
  });

  it("goes to the model's loci, one per kind and alignment, and back through Extents", () => {
    const loci = (text: string, o?: Parameters<typeof SeqSelection.toLoci>[1]) =>
      SeqSelection.toLoci(run(text).selection, o);
    const [a] = loci('(sel "1brs/A[10-20, 31-40]")');
    const want = {
      kind: "seq-loci",
      set,
      elements: [
        {
          unit: unit("1brs/A"),
          indices: OrderedSet.ofSortedArray(
            Int32Array.from([...Array(11).keys()].map((i) => i + 9).concat([...Array(10).keys()].map((i) => i + 30))),
          ),
        },
      ],
    };
    expect(Loci.areEqual(a, want as unknown as SeqLoci)).toBe(true);
    expect(Loci.areEqual(loci('(sel "1brs/A[10-20]")')[0], seqLoci([["1brs/A", 9, 20]]))).toBe(true);
    // one SeqLoci for every unit, then rows and columns per alignment, then items
    expect(loci("(or 1brs/D msa/2 msa/0 msa 1brs/A sites)").map((l) => l.kind)).toEqual([
      "seq-loci",
      "row-loci",
      "column-loci",
      "item-loci",
    ]);
    const [rows] = loci("(or msa/2 msa/0)");
    expect(rows.kind === "row-loci" && rows.elements.map((e) => e.row)).toEqual([0, 2]);
    const [items] = loci("(or (items :of sites) (items :of dom))");
    expect(items.kind === "item-loci" && items.elements.map((e) => [e.annotation, OrderedSet.size(e.items)])).toEqual([
      [0, 4],
      [1, 2],
    ]);
    expect(loci("sites", { items: "positions" }).map((l) => l.kind)).toEqual(["seq-loci"]);
    expect(loci("dom", { group: 1 })[0].kind).toBe("item-loci");
    // a molecule run on every alive unit over it; with none, dropped and counted
    const report = { dropped: 0 };
    expect(loci("(systems :kind molecule)", { report }).map((l) => l.kind)).toEqual(["seq-loci"]);
    expect(report.dropped).toBe(10);
    // back: the loci's extents are the selection's runs
    const s = run('(or 1brs/A msa/1 msa (sel "ab/L[3A-4]"))').selection;
    const there = Extents.create();
    for (const l of SeqSelection.toLoci(s)) Extents.addLoci(there, l);
    const back = SeqSelection.fromExtents(set, there);
    expect(groupsOf(back)).toEqual(groupsOf(SeqSelection.union(s)));
    const out = Extents.create(1);
    expect(SeqSelection.toExtents(s, out, 1)).toBe(1);
    expect([out.count, out.kind[0], out.stop[0]]).toEqual([1, 2, 6]);
  });

  it("takes a host's loci as one group each, and summarizes as people and the agent read it", () => {
    expect(groupsOf(SELECTION).map((g) => g.pos.length)).toEqual([5, 4]);
    expect(SeqSelection.summary(SELECTION)).toBe("2 groups, 9 positions: 1brs/A 41-45; ab/H 4-7");
    expect(SeqSelection.summary(run("(all)").selection, { maxParts: 2 })).toBe(
      "11 groups, 651 positions: 1brs/A 1-110; 1brs/B 1-110 (and 9 more)",
    );
    expect(SeqSelection.summary(SeqSelection.empty(set))).toBe("selects nothing");
    expect(summary("(input selection)")).toBe("2 groups, 9 positions: 1brs/A 41-45 @40:45; ab/H 4-7 @5:9");
  });
});

describe("extensions", () => {
  // A property and its runtime, as an extension ships them: the item slot's index, and a system's id.
  const nth = defineSymbol({
    name: "x.item.nth",
    alias: "nth",
    kind: "property",
    reads: "item",
    args: [],
    result: "int",
    signature: "(nth)",
    doc: "the item's index",
    examples: [],
    hints: [],
  });
  const env = createQueryEnv([
    { namespace: "x", symbols: [nth], runtimes: [{ symbol: nth, compile: () => (ctx) => ctx.slot.item }] },
  ]);

  it("run the core's logic and relations on the slots a :where binds", () => {
    const items = (where: string) => Array.from(run(`(items :of sites :where ${where})`, {}, env).selection.groupItem);
    expect(items("(>= (nth) 2)")).toEqual([2, 3]);
    expect(items("(and (> (nth) 0) (< (nth) 3))")).toEqual([1, 2]);
    expect(items("(or (= (nth) 0) (in (nth) 2 3))")).toEqual([0, 2, 3]);
    expect(items("(not (!= (nth) 1))")).toEqual([1]);
    expect(items("(<= (nth) 1)")).toEqual([0, 1]);
    expect(items("(matches s1 ^S :case true)")).toEqual([]);
    expect(items("(matches s1 ^S)")).toEqual([0, 1, 2, 3]);
    // missing values fail every comparison, and their negation holds (seqspec:seqql#labels-and-indices)
    expect(items('(< "x" 1)')).toEqual([]);
    expect(items('(not (= "x" 1))')).toEqual([0, 1, 2, 3]);
  });

  it("run only with their own definition object; others' are unavailable", () => {
    const copy = defineSymbol({ ...nth, name: "x.item.nth" });
    const other = createQueryEnv([
      { namespace: "x", symbols: [nth], runtimes: [{ symbol: copy, compile: () => () => 0 }] },
    ]);
    expect(other.runtimes.has("x.item.nth")).toBe(false);
    expect(run("(items :of sites :where (> (nth) 1))", {}, other).said[0]).toMatch(/seqql.symbol.unavailable/);
  });
});
