/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QL2b: the `sel` shorthand (seqspec:seqql#the-sel-shorthand): its grammar and expansion, the stored string normalized,
// and its errors at their column inside the string (escapes counted as written).
import { describe, expect, it } from "vitest";
import { checkText, type Expression, printQuery, Q, renderIssue } from "../../src/query/index.ts";
import { parseSel } from "../../src/query/sel.ts";

/** The expansion, printed. */
const expansion = (s: string) => {
  const r = parseSel(s);
  if (!("expr" in r)) throw new Error(r.problem.message);
  return printQuery(r.expr);
};
const rendered = (text: string) => checkText(text).issues.map(renderIssue);

describe("the grammar and its expansion", () => {
  it("expands every form of the table", () => {
    expect(
      [
        "1brs/A",
        "1brs/A[10-20,45,50-52]",
        "1brs/A[uniprot: 120-140]",
        "1brs/A@9:20",
        "1brs/A@9",
        "domains[#PF00069,#PF07714]",
        "$selection",
        "$selection[10-20]",
        "$selection@0:5",
        "1brs/A[10-20], 1brs/B@0:5, $focus",
      ].map(expansion),
    ).toEqual([
      "1brs/A",
      "(range 1brs/A 10-20 45 50-52)",
      "(range 1brs/A 120-140 :numbering uniprot)",
      "(slice 1brs/A 9 20)",
      "(slice 1brs/A 9 10)",
      "(ids domains PF00069 PF07714)",
      "(input selection)",
      "(range (input selection) 10-20)",
      "(slice (input selection) 0 5)",
      "(or (range 1brs/A 10-20) (slice 1brs/B 0 5) (input focus))",
    ]);
  });

  it("reads labels as the format writes them, numbers where the binder would read one", () => {
    const r = parseSel("p53[-5-10,45A,45,1.5,007,45.10,H1]");
    expect("expr" in r && r.expr).toEqual(Q.seq.modifier.range("p53", "-5-10", "45A", 45, 1.5, "007", "45.10", "H1"));
    expect(expansion("p53[-10--5]")).toBe("(range p53 -10--5)");
  });

  it("takes refs with colons, dots and slashes, and spaces after commas, inside brackets and at the ends", () => {
    expect(expansion("1brs/entity:1@0:5")).toBe("(slice 1brs/entity:1 0 5)");
    expect(expansion("ctx.site.hinge")).toBe("ctx.site.hinge");
    expect(expansion("  1brs/A[ uniprot :  10-20 ,  45 ] ,\n 1brs/B  ")).toBe(
      "(or (range 1brs/A 10-20 45 :numbering uniprot) 1brs/B)",
    );
  });

  it("is usable anywhere a selection is, and expands for checking", () => {
    const r = checkText('(pick (items :of variants) :where (overlaps (sel "domains[#PF00069]")))');
    expect(r.issues).toEqual([]);
    expect(printQuery(r.expanded as Expression)).toBe(
      "(pick (items :of variants) :where (overlaps (ids domains PF00069)))",
    );
    expect(printQuery(r.expr as Expression)).toBe(
      '(pick (items :of variants) :where (overlaps (sel "domains[#PF00069]")))',
    );
  });
});

describe("the stored string", () => {
  it("is normalized, spaces removed, by the binder and the builder", () => {
    const r = checkText('(onto (sel " 1brs/A[ 10-20, 45 ],  $selection ") msa)');
    expect(r.expr).toEqual(Q.seq.modifier.onto(Q.seq.macro.sel("1brs/A[10-20,45],$selection"), "msa"));
    expect(Q.seq.macro.sel("a[1, 2]")).toEqual({ head: "seq.macro.sel", args: { "0": "a[1,2]" } });
    expect(Q.read('(sel "a[1, 2]")').text).toBe('(sel "a[1,2]")');
  });

  it("is kept as written when it does not read, so its errors point into what was written", () => {
    expect(Q.seq.macro.sel("a [1]")).toEqual({ head: "seq.macro.sel", args: { "0": "a [1]" } });
    expect(rendered('(sel "a, b [1]")')).toEqual([
      '1:12 seqql.sel.syntax: unexpected "[": labels in brackets, indices after @: (sel "1brs/A[10-20]"), (sel "1brs/A@9:20")\n  (sel S)',
    ]);
  });
});

describe("errors inside the string", () => {
  it("names each problem at its column, and the unclosed bracket", () => {
    expect(
      [
        '(sel "1brs/A[10-20")',
        '(sel "1brs/A[10-]")',
        '(sel "1brs/A@9-20")',
        '(sel "1brs/A@:20")',
        '(sel "1brs/A@9:9")',
        '(sel "1brs/A[10,#PF1]")',
        '(sel "domains[#PF1,10]")',
        '(sel "domains[#]")',
        '(sel "1brs/A[]")',
        '(sel "$")',
        '(sel "a,")',
        '(sel "")',
        '(sel "p53-features[#DNA-binding domain]")',
        '(sel "domains[#PF1 #PF2]")',
      ].map((t) => rendered(t)[0]?.split("\n")[0]),
    ).toEqual([
      '1:19 seqql.sel.syntax: expected "]" to close "[" at 1:13',
      '1:17 seqql.sel.syntax: a range is two labels joined by "-": 10-20',
      '1:15 seqql.sel.syntax: unexpected "-": labels in brackets, indices after @: (sel "1brs/A[10-20]"), (sel "1brs/A@9:20")',
      "1:14 seqql.sel.syntax: a slice is start:stop, indices from 0: @9:20",
      "1:14 seqql.sel.syntax: a slice's ends are equal, which selects nothing: @9:20",
      "1:17 seqql.sel.mixed: a bracket holds labels or #ids, not both",
      "1:20 seqql.sel.mixed: a bracket holds labels or #ids, not both",
      "1:16 seqql.sel.syntax: an item id follows #: [#PF00069]",
      "1:14 seqql.sel.syntax: expected a label: 10-20, 45A",
      "1:8 seqql.sel.syntax: an input is $ and its name: $selection",
      '1:9 seqql.sel.syntax: the shorthand ends early: labels in brackets, indices after @: (sel "1brs/A[10-20]"), (sel "1brs/A@9:20")',
      '1:7 seqql.sel.syntax: the shorthand ends early: labels in brackets, indices after @: (sel "1brs/A[10-20]"), (sel "1brs/A@9:20")',
      "1:32 seqql.sel.syntax: an id is one word: #DNA-binding ends at the space; ids are the annotation's own, separated by commas; match labels with (pick (items :of A) :where (matches (label) …))",
      "1:19 seqql.sel.syntax: an id is one word: #PF1 ends at the space; ids are the annotation's own, separated by commas; match labels with (pick (items :of A) :where (matches (label) …))",
    ]);
  });

  it("counts escapes as written, on any line, with the issue at the call's path", () => {
    const r = checkText('(or 1brs/B\n  (sel "\\u0031brs/A[10-20"))');
    expect(r.issues.map((i) => [i.path, i.at, i.open, i.message])).toEqual([
      ["/args/1/0", { line: 2, column: 26 }, { line: 2, column: 20 }, 'expected "]" to close "[" at 2:20'],
    ]);
  });

  it("takes only a string: the checker says so before expanding", () => {
    expect(rendered("(sel 45)")).toEqual([
      '1:6 seqql.type: sel takes the shorthand as a string: (sel "1brs/A[10-20]")\n  (sel S)',
    ]);
    expect(rendered("(sel)")[0]).toMatch(/^1:1 seqql.arg.missing: sel takes the shorthand as a string/);
  });
});
