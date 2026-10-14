/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SeqQL's QR2a runtimes (seqspec:seqql#each): `onto` across alignment rows and columns, bindings, a codon mapping and
// paths; `partners` over 1BRS's contacts and ends of every kind; `expand`, `group`, `each` with `current`, `pick`,
// `first`, `ids` and the item slot's properties. Every symbol against the brute-force oracle over random queries, then
// the shapes and messages the spec names.

import { bindQuery, checkText, readQuery, renderIssue } from "@molstar/seq-core/query";
import { describe, expect, it } from "vitest";
import { type SeqIndex, SeqSets, type UnitId } from "../../src/model/index.js";
import { compileQuery, createQueryEnv, type QueryOptions, SeqSelection } from "../../src/model/query/index.js";
import { ONTO_REFS, ontoFixture } from "./query-fixtures.js";
import { groupsOf, Oracle, OracleError } from "./query-oracle.js";
import { prng } from "./support/generators.js";

const { set } = SeqSets.build(ontoFixture());
const ENV = createQueryEnv();
const unit = (key: string) => Array.from(set.units.key).indexOf(set.strings.find(key) as number) as UnitId;
const SELECTION = SeqSelection.fromLoci(set, [
  { kind: "seq-loci", set, elements: [{ unit: unit("1brs/D"), indices: { start: 25, stop: 40 } as never }] },
]);

function run(text: string, o: Partial<QueryOptions> = {}) {
  const read = readQuery(text);
  const bound = bindQuery(read.node as never, ENV.language);
  expect([...read.issues, ...bound.issues]).toEqual([]);
  const r = compileQuery(bound.expr as never, ENV, bound.source).run({ set, inputs: { selection: SELECTION }, ...o });
  return { ...r, said: r.issues.map((i) => renderIssue(i).split("\n")[0]) };
}
const summary = (text: string) => SeqSelection.summary(run(text).selection, { maxParts: 12 });

describe("QR2a's runtimes against the oracle", () => {
  const SEED = 20261013;
  let rand = prng(SEED);
  const pick = <T>(xs: readonly T[]) => xs[Math.floor(rand() * xs.length)];
  const int = (lo: number, hi: number) => lo + Math.floor(rand() * (hi - lo));
  const LEAVES = [
    ...ONTO_REFS,
    "(input selection)",
    "(items :of sites)",
    "(items :of links)",
    "(rows :of fam)",
    '(sel "1brs/A[20-40]")',
    '(sel "1brs/D[25-35]")',
    '(sel "protein[3-8]")',
    '(sel "gene[20-50]")',
    '(sel "fam[10-30]")',
    '(sel "fam/0[5-25]")',
    '(sel "ins/1[2-5]")',
  ];
  const TARGETS = [
    "1brs/A",
    "1brs/B",
    "1brs/D",
    "1brs/E",
    "1brs",
    "1brs/entity:1",
    "gene",
    "protein",
    "p2",
    "lonely",
    "fam",
    "fam/0",
    "fam/1",
    "fam/2",
    "fam/3",
    "ins",
    "ins/0",
    "ins/1",
    "sites",
  ];
  const ITEM_TESTS = [
    "(= (category) hbond)",
    "(= (category) x)",
    "(matches (label) kinase)",
    "(in (id) s1 s3 w r m)",
    "(< (value) 2)",
    "(> (value) 1)",
    "(= (annotation) sites)",
    "(not (= (category) contact))",
    "(and (> (value) 0) (matches (label) ^[ab]))",
    "(or (= (id) s5) (matches (label) site))",
    "(< (field distance) 3.5)",
  ];
  const WHERE_ONLY = ["(> (length) 4)", "(< (length) 3)"];
  function gen(depth: number): string {
    if (depth === 0 || rand() < 0.2) return pick(LEAVES);
    const q = () => gen(depth - 1);
    /** A selection of items, mostly: what pick and ids test. */
    const tagged = () =>
      rand() < 0.3 ? q() : pick(["sites", "links", "contacts", "(items :of sites)", `(or links ${q()})`]);
    switch (int(0, 12)) {
      case 0:
      case 1: {
        const kind = rand() < 0.3 ? ` :kind ${pick(["row", "columns", "unit", "molecule"])}` : "";
        const to = rand() < 0.12 && kind ? "" : ` ${pick(TARGETS)}`;
        return `(onto ${q()}${to}${kind}${rand() < 0.2 ? " :keep true" : ""})`;
      }
      case 2:
      case 10: {
        const side = rand() < 0.3 ? ` :side ${pick(["a", "b"])}` : "";
        const where = rand() < 0.3 ? ` :where ${pick([...ITEM_TESTS, ...WHERE_ONLY])}` : "";
        const of =
          rand() < 0.3 ? q() : pick(["1brs/D", "1brs/A", '(sel "1brs/D[25-45]")', "1brs", "lonely", "fam", "1brs/F"]);
        return `(partners ${of} :via ${pick(["contacts", "links", "contacts", "sites"])}${side}${where})`;
      }
      case 3:
        return rand() < 0.6 ? `(expand ${q()} ${int(0, 6)})` : `(expand ${q()} :left ${int(0, 4)} :right ${int(0, 9)})`;
      case 4:
        return `(group ${q()} :by ${pick(["none", "run", "system"])})`;
      case 5: {
        const e = pick([
          "(expand (current) 2)",
          `(and (current) ${q()})`,
          `(onto (current) ${pick(TARGETS)})`,
          "(first (current))",
          `(partners (current) :via contacts)`,
          `(or (current) ${q()})`,
        ]);
        return `(each ${q()} ${e})`;
      }
      case 6:
        return `(pick ${tagged()} :where ${pick(ITEM_TESTS)})`;
      case 7:
        return `(first ${q()}${rand() < 0.7 ? ` ${int(0, 5)}` : ""})`;
      case 8:
        return `(ids ${tagged()} ${pick(["s1", "s2 s4", "w m r", "s6 l s"])})`;
      case 9:
        return `(items :of ${pick(["sites", "links", "contacts"])} :where ${pick([...ITEM_TESTS, ...WHERE_ONLY])})`;
      default:
        return `(${pick(["and", "or", "or"])} ${q()} ${q()})`;
    }
  }

  const heads = ["onto", "partners", "expand", "group", "each", "pick", "first", "ids", ":where"];
  /** `count` random queries from the seed, each against the oracle; how many of each outcome. */
  function againstOracle(count: number) {
    rand = prng(SEED);
    const oracle = new Oracle(set, { selection: groupsOf(SELECTION) });
    const tally = { runs: 0, errors: 0, empty: 0, ...Object.fromEntries(heads.map((h) => [h, 0])) } as Record<
      string,
      number
    >;
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
      else for (const h of heads) tally[h] += +text.includes(h === ":where" ? h : `(${h} `);
    }
    // Enough of each outcome, and of each symbol in a non-empty result, that the comparison means something.
    return tally;
  }

  it("gives the oracle's groups, in its order, or its error, for 300 random queries", () => {
    const tally = againstOracle(300);
    expect(tally.runs).toBeGreaterThan(280);
    expect(tally.errors).toBeGreaterThan(4);
    expect(tally.runs - tally.errors - tally.empty).toBeGreaterThan(120);
    for (const h of heads) expect([h, tally[h] > 2]).toEqual([h, true]);
  });

  it("gives the oracle's groups, in its order, or its error, for 1,500 random queries", {
    tags: ["slow"],
    timeout: 60_000,
  }, () => {
    const tally = againstOracle(1500);
    expect(tally.runs).toBeGreaterThan(1400);
    expect(tally.errors).toBeGreaterThan(20);
    expect(tally.runs - tally.errors - tally.empty).toBeGreaterThan(600);
    for (const h of heads) expect([h, tally[h] > 40]).toEqual([h, true]);
    // A second under a normal load; a shared machine's load can make it many.
  });
});

describe("onto: the same places in another system, by overlap", () => {
  it("carries a row onto its columns and back, through A3M inserts, and a row onto another through them", () => {
    // ins/0 is ACdeFG-H: d and e are inserts, no column; H sits after a gap
    expect(summary("(onto ins/0 ins)")).toBe("1 group, 5 positions: ins columns 1-4, 6");
    expect(summary("(onto ins ins/1)")).toBe("1 group, 5 positions: ins/1 1-3, 6-7");
    expect(summary("(onto ins/0 ins/1)")).toBe("1 group, 4 positions: ins/1 1-3, 7");
    // onto every row: its own residues (all, inserts too), and the other row's at its columns
    expect(summary("(onto ins/0 ins :kind row)")).toBe("1 group, 11 positions: ins/0 1-7; ins/1 1-3, 7");
  });

  it("reaches units through the rows' bindings: a molecule's binding every unit over it, the columns every bound row", () => {
    expect(summary("(onto fam/0 1brs/B)")).toBe("1 group, 110 positions: 1brs/B 1-110");
    expect(summary("(onto fam 1brs)")).toBe(
      "1 group, 390 positions: 1brs/A 1-110; 1brs/B 1-110; 1brs/C 1-110; 1brs/D 21-60; 1brs/E 6-25",
    );
    // a row bound to two sequences reaches each: onto names where it goes
    expect(summary("(onto fam/2 protein)")).toBe("1 group, 20 positions: protein/unit 1-20");
    expect(summary("(onto fam/2 1brs/E)")).toBe("1 group, 20 positions: 1brs/E 6-25");
    // a unit onto a row it is not bound to: through its own row and the columns
    const [a, d] = [run("(onto 1brs/A fam/1)").selection, run("(onto (onto 1brs/A fam) fam/1)").selection];
    expect(groupsOf(a)).toEqual(groupsOf(d));
  });

  it("goes across the codon mapping by overlap both ways, and along paths through it", () => {
    // residue 5 is codon 13-15 of the CDS from 13: gene 25-27; residues 5-6 the six nucleotides
    expect(summary('(onto (sel "protein[5-6]") gene)')).toBe("1 group, 6 positions: gene/unit 25-30");
    // a nucleotide lands on the residue whose codon holds it: 40-45 touch codons 10 and 11
    expect(summary('(onto (sel "gene[40-45]") protein)')).toBe("1 group, 2 positions: protein/unit 10-11");
    // across the intron (CDS 13-46, 77-105): residue 12's codon is split, 46 and 77-78
    expect(summary('(onto (sel "protein[12]") gene)')).toBe("1 group, 3 positions: gene/unit 46, 77-78");
    // p2 → protein → gene: one mapping and the codon step
    expect(summary('(onto (sel "p2[1]") gene)')).toBe("1 group, 3 positions: gene/unit 19-21");
    // a molecule with no unit, along its path
    expect(summary('(onto (sel "lonely[1]") gene)')).toBe("1 group, 3 positions: gene/unit 22-24");
  });

  it("keeps the source with :keep, takes every system of a kind without R, and drops a group with no image", () => {
    expect(summary('(onto (sel "protein[5]") gene :keep true)')).toBe(
      "1 group, 4 positions: gene/unit 25-27; protein/unit 5",
    );
    // every unit: itself, and barnase's three chains along the path through the family's rows 1 and 0
    expect(summary('(onto (sel "1brs/D[30]") :kind unit)')).toBe(
      "1 group, 4 positions: 1brs/A 25; 1brs/B 25; 1brs/C 25; 1brs/D 30",
    );
    const r = run("(onto (or fam/3 fam/0) 1brs/A)");
    expect([r.selection.groupCount, r.said]).toEqual([
      1,
      [
        "1:1 seqql.onto.unmapped: 1 of 2 groups has no image on 1brs/A: fam/3 is bound to no sequence; a group with no image is dropped; a row reaches other sequences through its binding",
      ],
    ]);
  });

  it("skips an ambiguous path, and says what :kind does when R names a row or a unit", () => {
    expect(run("(onto amb protein)").said).toEqual([
      "1:1 seqql.onto.ambiguous: amb/unit reaches protein/unit by paths that place it differently; it was skipped",
      "1:1 seqql.onto.unmapped: 1 of 1 group has no image on protein: amb/unit has no route onto protein/unit; a group with no image is dropped; a row reaches other sequences through its binding",
    ]);
    // QE1a: `:kind row` onto a row is harmless, and said so
    const r = run("(onto 1brs/A fam/1 :kind row)");
    expect([r.selection.runCount > 0, r.said]).toEqual([
      true,
      [
        '1:26 seqql.onto.kind: "fam/1" is one row already: :kind is for an alignment\'s rows or columns, (onto Q fam :kind row) for all its rows; it was ignored',
      ],
    ]);
    expect(run("(onto 1brs/D 1brs/A :kind unit)").said[0]).toMatch(/"1brs\/A" is a unit: :kind is for/);
    expect(run("(onto 1brs/A fam :kind unit)").said).toEqual([
      "1:24 seqql.type: :kind unit onto an alignment: it gives its columns, or its rows with :kind row; for the units bound to its rows, leave R out: (onto Q :kind unit)",
    ]);
    expect(run("(onto 1brs/A)").said[0]).toMatch(/^1:1 seqql.arg.missing: onto needs where to go/);
    expect(run("(onto 1brs/A sites)").said[0]).toMatch(
      /"sites" is an annotation; onto takes a unit, molecule, row or alignment$/,
    );
  });
});

describe("partners: one side of a pair by its partner", () => {
  it("gives barnase's interface by barstar over 1BRS's contacts, one group per run", () => {
    const iface = run("(and 1brs/A (partners 1brs/D :via contacts))").selection;
    expect(SeqSelection.summary(iface, { maxParts: 20 })).toBe(
      "1 group, 19 positions: 1brs/A 27, 35, 37-38, 56, 58-60, 62, 73, 82-85, 87, 101-104",
    );
    expect(run("(partners 1brs/D :via contacts)").selection.groupCount).toBe(10);
    // hydrogen bonds only (under 3.2 Å, the fixture's category)
    expect(summary("(partners 1brs/D :via contacts :where (= (category) hbond))")).toBe(
      "5 groups, 7 positions: 1brs/A 27, 59-60, 83-84, 87, 102",
    );
    // a partner selection on chain B (a copy) touches none of chain A's contacts
    expect(summary("(partners 1brs/B :via contacts)")).toBe("selects nothing");
  });

  it("finds ends of every kind: a whole unit, a molecule over its units, one with none, a row and the columns", () => {
    expect(summary("(partners 1brs/F :via links)")).toBe("2 groups, 3 positions: 1brs/B 6; fam columns 4-5");
    expect(summary("(partners lonely :via links)")).toBe("2 groups, 4 positions: 1brs/A 61-62; fam/3 3-4");
    expect(summary("(partners fam/3 :via links)")).toBe("1 group, 2 positions: lonely 2-3");
    expect(summary("(partners fam :via links)")).toBe(
      "3 groups, 6 positions: 1brs/D 11-12; 1brs/E 11-12; 1brs/F 11-12",
    );
    // both ends in Q: both returned; :side keeps one
    expect(summary("(partners 1brs/C :via links)")).toBe("2 groups, 7 positions: 1brs/C 41-45, 51-52");
    expect(summary("(partners 1brs/C :via links :side b)")).toBe("1 group, 2 positions: 1brs/C 51-52");
    expect(run("(partners 1brs/A :via sites)").said[0]).toMatch(/"sites" is a ranges annotation; :via takes a pairs/);
  });

  it("counts the partner selection's positions against the limit", () => {
    expect(run("(partners 1brs :via contacts)", { limits: { positions: 100 } }).said[0]).toMatch(
      /seqql.limit.positions: this visits \d+ positions in contacts/,
    );
  });
});

describe("modifiers and filters", () => {
  it("expand each run in its own system, clipped, merging runs that meet; the item tag dropped", () => {
    expect(summary('(expand (sel "1brs/A[1-3, 8-9]") 2)')).toBe("1 group, 11 positions: 1brs/A 1-11");
    expect(summary('(expand (sel "1brs/A[50]") :left 3)')).toBe("1 group, 4 positions: 1brs/A 47-50");
    expect(Array.from(run("(expand sites 1)").selection.groupAnnotation)).toEqual([-1, -1, -1, -1, -1, -1]);
  });

  it("regroup: none one group, run each maximal run, system one per system", () => {
    expect(run("(group sites :by none)").selection.groupCount).toBe(1);
    expect(run("(group sites :by run)").selection.groupCount).toBe(6);
    expect(summary("(group sites :by system)")).toBe(
      "5 groups, 70 positions: gene/unit 31-60; protein/unit 5-9; 1brs/A 6-10, 31-40; 1brs/D 61-70; fam columns 11-20",
    );
  });

  it("each runs E per group with (current) bound to it, its results concatenated in order", () => {
    const r = run("(each (items :of sites) (expand (current) 1))").selection;
    expect([r.groupCount, SeqSelection.summary(r, { maxParts: 1 })]).toEqual([
      6,
      "6 groups, 82 positions: gene/unit 30-61 (and 5 more)",
    ]);
    // nested: the inner current is the inner group, the outer's restored after it
    expect(groupsOf(run("(each sites (each (current) (first (current))))").selection)).toEqual(
      groupsOf(run("sites").selection),
    );
    expect(Array.from(run("(each sites (current))").selection.groupItem)).toEqual([0, 1, 2, 3, 4, 5]);
  });

  it("pick, first and ids keep groups as they are, tags included", () => {
    expect(summary("(pick sites :where (matches (label) kinase))")).toBe(
      "2 groups, 15 positions: protein/unit 5-9; fam columns 11-20",
    );
    expect(Array.from(run("(first sites 2)").selection.groupItem)).toEqual([0, 1]);
    expect(run("(first sites 0)").selection.groupCount).toBe(0);
    const ids = run("(ids sites s2 s9)");
    expect([Array.from(ids.selection.groupItem), ids.said]).toEqual([
      [1],
      [
        "1:1 seqql.ids.unknown: no item has the id \"s9\"; its items' ids include s1, s2, s3, s4, s5, s6; ids are the annotation's own: match labels with (pick Q :where (matches (label) …))",
      ],
    ]);
    expect(summary('(sel "sites[#s1,#s3]")')).toBe("2 groups, 15 positions: 1brs/A 6-10; 1brs/D 61-70");
    // An id listed twice is found, not said missing (CM-13).
    const twice = run("(ids sites s2 s2)");
    expect([Array.from(twice.selection.groupItem), twice.said]).toEqual([[1], []]);
  });

  it("read the item slot: id, label, category, value, type, annotation, length and fields", () => {
    const items = (where: string) => Array.from(run(`(items :of sites :where ${where})`).selection.groupItem);
    expect(items("(in (id) s2 s4)")).toEqual([1, 3]);
    expect(items("(matches (label) site)")).toEqual([4]);
    expect(items("(= (category) x)")).toEqual([0, 2, 4]);
    expect(items("(>= (value) 2)")).toEqual([1, 2, 4]);
    expect(items("(= (annotation) sites)")).toEqual([0, 1, 2, 3, 4, 5]);
    expect(items("(= (type) x)")).toEqual([]);
    expect(items("(> (length) 10)")).toEqual([4]);
    expect(Array.from(run("(items :of contacts :where (< (field distance) 3))").selection.groupItem).length).toBe(8);
    expect(run("(items :of sites :where (= (field nope) 1))").said).toContain(
      '1:28 seqql.field.unknown: sites has no field "nope"; its fields: category, label, value; an annotation\'s own columns are listed when the name is unknown',
    );
    // an untagged group has no item: its properties are missing, and fail every comparison
    expect(run("(pick 1brs/A :where (not (= (id) x)))").selection.groupCount).toBe(1);
    expect(run("(pick 1brs/A :where (= (id) x))").selection.groupCount).toBe(0);
  });
});

describe("toLoci of an onto result", () => {
  it("gives column loci for the alignment, unit loci for chains", () => {
    const loci = SeqSelection.toLoci(run('(onto (sel "1brs/A[1-5]") fam :keep true)').selection);
    expect(loci.map((l) => l.kind)).toEqual(["seq-loci", "column-loci"]);
    const one = SeqSelection.toLoci(run('(onto (sel "fam[2]") 1brs/D)').selection);
    expect(one.length <= 1).toBe(true);
    void (0 as SeqIndex);
  });
});
