/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QL2b: the checker (seqspec:seqql#types, seqspec:seqql#errors): types, arity, slice's pairs, slots of full names,
// (current), fields; stored JSON checked as text is; macros checked and expanded; Q.read.
import { describe, expect, it } from "vitest";
import {
  checkQuery,
  checkText,
  createLanguage,
  defineSymbol,
  type Expression,
  printQuery,
  Q,
  type QueryError,
  renderIssue,
} from "../../src/query/index.js";

const rendered = (text: string, language?: Parameters<typeof checkText>[1]) =>
  checkText(text, language).issues.map(renderIssue);
const first = (text: string) => rendered(text)[0]?.split("\n")[0];

describe("types", () => {
  it("checks written values against each argument's type, with what was given and how to fix it", () => {
    expect(
      [
        "(expand 1brs/A -1)",
        "(systems :kind units)",
        "(range 1brs/A 10..20)",
        "(onto 123 msa)",
        "(onto 1brs/A 7)",
        "(pick 1brs/A :where 1brs/B)",
        "(search 123)",
        "(onto 1brs/A msa :keep yes)",
        "(first 1brs/A 1.5)",
        "(partners 1brs/A :via contacts :side D)",
        '(onto (sel "protein[5-10]") gene :kind codon)',
        "(onto 1brs/A gene :kind Codons)",
      ].map(first),
    ).toEqual([
      "1:16 seqql.type: expand's n takes a whole number from 0, not -1",
      '1:16 seqql.type: :kind takes one of unit, molecule, row, columns, not "units"; did you mean "unit"?',
      '1:15 seqql.type: range\'s part takes labels or label ranges: 10-20, 45A, -5-10, not "10..20"',
      '1:7 seqql.type: onto\'s Q takes a selection, not 123; a ref that reads as a number is quoted: "123"',
      '1:14 seqql.type: onto\'s R takes a ref, not 7; a ref that reads as a number is quoted: "7"',
      '1:21 seqql.type: :where takes a test, as (> (length) 50), not "1brs/B"',
      "1:9 seqql.type: search's P takes text, not 123",
      '1:24 seqql.type: :keep takes true or false, not "yes"',
      "1:15 seqql.type: first's n takes a whole number from 0, not 1.5",
      // capitals name something else (QE1a: Luna's chain): the argument's doc, not a did-you-mean
      '1:38 seqql.type: :side takes one of any, a, b, not "D" (the values are lower-case): which ends may be returned: a or b, the pairs\' columns (any: both); not a chain, which goes in Q',
      // a value the checker knows to be a guess (QE1a, QL2c, QE1b: models write it for the codons of a mapping)
      "1:40 seqql.type: :kind is a kind of sequence (unit, molecule, row, columns), and codon is not one: across a codon mapping leave :kind out, (onto Q gene) gives the codons",
      "1:25 seqql.type: :kind is a kind of sequence (unit, molecule, row, columns), and codon is not one: across a codon mapping leave :kind out, (onto Q gene) gives the codons",
    ]);
  });

  it("checks what a call gives against where it goes", () => {
    expect(
      [
        "(pick (items :where (> (length) 3)) :where (items))",
        "(items :where (= (all) 3))",
        "(positions :where (matches (index) 1))",
        "(onto 1brs/A (all))",
      ].map(first),
    ).toEqual([
      "1:44 seqql.type: :where takes a test, as (> (length) 50), not a selection (items)",
      "1:18 seqql.type: ='s a takes a value, not a selection (all)",
      "1:28 seqql.type: matches's s takes text, not a number (index)",
      "1:14 seqql.type: onto's R takes a ref, not a selection (all)",
    ]);
    // a value of unknown type (a field, a value) fits a number, a test or text
    expect(rendered("(items :where (and (> (value) 0.8) (matches (field tooltip) x)))")).toEqual([]);
  });

  it("uses the symbol's hint for a type error, except about a selection", () => {
    expect(first("(items :where (< (label) 3))")).toBe("1:18 seqql.type: < compares numbers: (< (length) 50)");
    expect(first("(slice 1brs/A 45A 50)")).toBe(
      "1:15 seqql.type: slice takes indices, whole numbers from 0: for labels like 45A use (range 1brs/A 45A)",
    );
    expect(first("(slice (> 1 2) 1 2)")).toBe("1:8 seqql.type: slice's Q takes a selection, not a test (> …)");
  });
});

describe("arity", () => {
  it("names a missing argument, by the hint when the symbol has one", () => {
    expect(
      ["(letters)", "(values :above 3)", "(group 1brs/A)", "(group)", "(and 1brs/A)", "(ids domains)", "(input)"].map(
        (t) => rendered(t),
      ),
    ).toEqual([
      ["1:1 seqql.arg.missing: letters needs the letters: (letters C)\n  (letters L :in Q)"],
      [
        "1:1 seqql.arg.missing: values needs :of: (values :of A :min x :max x :above x :below x)\n  (values :of A :min x :max x :above x :below x)",
      ],
      ["1:1 seqql.arg.missing: group needs :by: (group Q :by none)\n  (group Q :by G)"],
      [
        "1:1 seqql.arg.missing: group needs the selection: (group Q :by G)\n  (group Q :by G)",
        "1:1 seqql.arg.missing: group needs :by: (group Q :by none)\n  (group Q :by G)",
      ],
      ["1:1 seqql.arg.missing: and needs the selections they are intersected with: (and Q Q...)\n  (and Q Q...)"],
      ["1:1 seqql.arg.missing: ids needs at least one id: (ids (items :of domains) PF00069)\n  (ids Q id...)"],
      [
        "1:1 seqql.arg.missing: input needs the input's name: selection, focus, structure: (input name)\n  (input name)",
      ],
    ]);
  });

  it("takes slice's indices in start-stop pairs; equal ends select nothing; a start after its stop is the run's", () => {
    expect(rendered("(slice 1brs/A 9 20 30)")).toEqual([
      "1:20 seqql.arg.count: slice takes start-stop pairs: (slice 1brs/A 9 20)\n  (slice Q start stop...)",
    ]);
    expect(first("(slice 1brs/A 9 20 30 30)")).toBe(
      "1:20 seqql.slice.range: slice 30 30 selects nothing: the stop comes after the start, as (slice Q 9 20)",
    );
    // It wraps a cyclic system's origin, an error elsewhere: known when it runs (CM-10).
    expect(rendered("(slice 1brs/A 9 20 30 4)")).toEqual([]);
    expect(rendered("(slice 1brs/A 0 1 5 9)")).toEqual([]);
  });
});

describe("slots and (current)", () => {
  it("checks the slot a full name reads, which the binder accepts anywhere", () => {
    expect(rendered("(items :where (= (seq.system.label) x))")).toEqual([
      "1:18 seqql.slot: label reads a system: use it in (systems :where …) or (rows :where …)\n  (label)",
    ]);
    expect(rendered("(rows :where (= (seq.system.label) x))")).toEqual([]);
    expect(first("(seq.group.overlaps 1brs/A)")).toBe(
      "1:1 seqql.slot: overlaps reads a group: use it in (pick Q :where …)",
    );
    // a selection inside a test reads no slot: a test inside it needs its own binder
    expect(first("(pick (all) :where (overlaps (letters (seq.position.letter))))")).toBe(
      "1:39 seqql.slot: letter reads a position: use it in (positions :where …)",
    );
  });

  it("allows (current) only inside each's query", () => {
    expect(rendered("(each (items :of d) (and (current) (letters C)))")).toEqual([]);
    expect(rendered("(current)")).toEqual(["1:1 seqql.current: (current) is only meaningful inside each\n  (current)"]);
    expect(first("(each (expand (current) 5) (all))")).toBe(
      "1:15 seqql.current: (current) is only meaningful inside each",
    );
  });
});

describe("fields", () => {
  it("knows a system's own fields, and any other with :of or :ns", () => {
    expect(first("(systems :where (= (field auth) A))")).toBe(
      '1:27 seqql.field.unknown: a system has no field "auth"; a system\'s own fields: molecule, complex, alignment, row; others need :of or :ns',
    );
    expect(first("(systems :where (= (field molecle) p53))")).toBe(
      '1:27 seqql.field.unknown: a system has no field "molecle"; did you mean "molecule"?',
    );
    expect(rendered("(systems :where (= (field auth_asym_id :ns mvs) A))")).toEqual([]);
    expect(rendered("(rows :where (>= (field identity :of stats) 0.9))")).toEqual([]);
  });

  it("warns of an item field one edit from a property; an item's columns are the runtime's", () => {
    const r = checkText("(items :where (= (field lable) x))");
    expect(r.issues.map((i) => [i.severity, renderIssue(i)])).toEqual([
      ["warning", '1:25 seqql.field.near: an item has no field "lable"; did you mean "label"?\n  (field name)'],
    ]);
    expect(r.expanded).not.toBeNull();
    expect(Q.read("(items :where (= (field lable) x))").text).toBe("(items :where (= (field lable) x))");
  });
});

describe("stored JSON", () => {
  it("checks what no binder made: names, keys, the shape of rest arguments", () => {
    const bad = {
      head: "seq.filter.pick",
      args: { "0": { head: "seq.generator.item" }, "1": "x", by: "x", where: [{ head: "core.rel.eq" }] },
    } as unknown as Expression;
    expect(checkQuery(bad).issues.map((i) => [i.path, i.code, i.message])).toEqual([
      ["/args/0", "seqql.symbol.unknown", 'unknown symbol "seq.generator.item"; did you mean "seq.generator.items"?'],
      ["/args/1", "seqql.arg.count", "pick has no argument 1"],
      ["/args/by", "seqql.arg.unknown", "pick has no argument by"],
      ["/args/where", "seqql.type", ":where is one value, or a rest a list"],
    ]);
    expect(checkQuery(Q.seq.filter.pick("1brs/A", { where: true })).issues).toEqual([]);
  });

  it("gives the expansion for the runtime, and null on an error", () => {
    const e = Q.seq.modifier.onto(Q.seq.macro.sel("1brs/A[10-20]"), "msa");
    expect(checkQuery(e).expr).toEqual(Q.seq.modifier.onto(Q.seq.modifier.range("1brs/A", "10-20"), "msa"));
    expect(checkQuery(Q.seq.generator.current()).expr).toBeNull();
  });
});

describe("macros", () => {
  // Stories' chain: (chain A) is (systems :where (= (field auth_asym_id :ns mvs) A)); `(chain)` with a number
  // expands to a field compared with a selection here, to show where an expansion's issue goes.
  const chain = defineSymbol({
    name: "stories.chain",
    alias: "chain",
    kind: "macro",
    args: [{ name: "C", type: "str", positional: true, doc: "the chain" }],
    result: "query",
    signature: "(chain C)",
    doc: "units of chain C",
    examples: [],
    hints: [],
    expand: (call) => {
      const c = call.args?.["0"] as string;
      const value = c === "ALL" ? Q.seq.generator.all() : c;
      return Q.seq.generator.systems({
        where: Q.core.rel.eq(Q.seq.system.field("auth_asym_id", { ns: "mvs" }), value),
      });
    },
  });
  const { language } = createLanguage([{ symbols: [chain] }]);

  it("checks a macro's own arguments, then its expansion, placing the expansion's issues at the call", () => {
    expect(rendered("(and 1brs (chain 7))", language)).toEqual([
      "1:18 seqql.type: chain's C takes text, not 7\n  (chain C)",
    ]);
    expect(rendered("(and 1brs\n  (chain ALL))", language)).toEqual([
      "2:3 seqql.type: ='s b takes a value, not a selection (all)\n  (= a b)",
    ]);
    const r = checkText("(and 1brs (chain A))", language);
    expect([r.issues, printQuery(r.expanded as Expression, language)]).toEqual([
      [],
      "(and 1brs (systems :where (= (field auth_asym_id :ns mvs) A)))",
    ]);
  });
});

describe("Q.read", () => {
  it("reads, binds and checks, giving what a document stores", () => {
    expect(Q.read('(onto (sel "1brs/A[10-20]") msa)')).toEqual({
      seqql: "1.0",
      expr: Q.seq.modifier.onto(Q.seq.macro.sel("1brs/A[10-20]"), "msa"),
      text: '(onto (sel "1brs/A[10-20]") msa)',
    });
  });

  it("throws with the issues as rendered, from any stage", () => {
    const fail = (text: string) => {
      try {
        Q.read(text);
      } catch (e) {
        return [(e as QueryError).message, (e as QueryError).issues.length];
      }
      return null;
    };
    expect(fail("(items")).toEqual(['1:7 seqql.read.unclosed: unclosed "(" opened at 1:1: add ")"', 1]);
    expect(fail("(item)")).toEqual([
      '1:2 seqql.symbol.unknown: unknown symbol "item"; did you mean "items"?\n  (items :of A :where T)',
      1,
    ]);
    expect(fail("(and (current) (slice 1brs/A 1))")).toEqual([
      "1:6 seqql.current: (current) is only meaningful inside each\n  (current)\n1:30 seqql.arg.count: slice takes start-stop pairs: (slice 1brs/A 9 20)\n  (slice Q start stop...)",
      2,
    ]);
  });
});

describe("a text property compared with a number (CM-06)", () => {
  it("is seqql.type saying to quote it, for =, != and in; a quoted number checks clean", () => {
    const said = (t: string) => checkText(t).issues.map((i) => `${i.code}: ${i.message}`);
    expect(said("(positions :in p :where (= (label) 315))")).toEqual([
      'seqql.type: (label) is text: quote the number to compare, "315"',
    ]);
    expect(said("(positions :in p :where (in (label) 10 11))")).toEqual([
      'seqql.type: (label) is text: quote the number to compare, "10"',
      'seqql.type: (label) is text: quote the number to compare, "11"',
    ]);
    expect(said("(items :of sites :where (!= 7 (id)))")).toEqual([
      'seqql.type: (id) is text: quote the number to compare, "7"',
    ]);
    expect(said('(positions :in p :where (= (label) "315"))')).toEqual([]);
    expect(said("(positions :in p :where (= (index) 315))")).toEqual([]);
  });
});
