/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QL2a: the reader (seqspec:seqql#the-reader) and its errors (seqspec:seqql#errors), each rendered as `line:column
// code: message`.
import { describe, expect, it } from "vitest";
import { MAX_DEPTH, type QueryNode, readQuery, renderIssue } from "../../src/query/index.js";
import { MAX_TEXT } from "../../src/query/read.js";

/** A node as compact data: lists as arrays, keywords as ":name", atoms as values, `@line:column` on lists. */
function show(n: QueryNode): unknown {
  if (n.kind === "list") return [`@${n.at.line}:${n.at.column}`, ...n.items.map(show)];
  return n.kind === "keyword" ? `:${n.name}` : n.value;
}
const node = (text: string) => {
  const r = readQuery(text);
  expect(r.issues).toEqual([]);
  return r.node as QueryNode;
};
const error = (text: string) => {
  const r = readQuery(text);
  expect(r.node).toBeNull();
  expect(r.issues).toHaveLength(1);
  return r.issues[0] as NonNullable<(typeof r.issues)[0]>;
};
const rendered = (text: string) => renderIssue(error(text));

describe("reading", () => {
  it("reads lists, words, strings, numbers, booleans and keywords", () => {
    const n = node('(pick (items :of domains) :where (matches (label) "kinase domain" :case true))');
    expect(show(n)).toEqual([
      "@1:1",
      "pick",
      ["@1:7", "items", ":of", "domains"],
      ":where",
      ["@1:34", "matches", ["@1:43", "label"], "kinase domain", ":case", true],
    ]);
    expect(show(node("(f 45 -5 1.5e3 0.25 007 1. -5-10 45A false)"))).toEqual([
      "@1:1",
      "f",
      45,
      -5,
      1500,
      0.25,
      "007",
      "1.",
      "-5-10",
      "45A",
      false,
    ]);
  });

  it("reads refs as words: '/', ':' after the first character, '.', '-', '@'", () => {
    const refs = show(node("(f 1brs/A 1brs/entity:1 msa/3 ctx.site.hinge p53@0:50 uniprot:Domain)")) as unknown[];
    expect(refs.slice(2)).toEqual(["1brs/A", "1brs/entity:1", "msa/3", "ctx.site.hinge", "p53@0:50", "uniprot:Domain"]);
  });

  it("reads JSON strings with their escapes, keeping the raw text", () => {
    const n = node('(f "a\\"b\\\\c\\n\\u00e9" "123" "true")');
    expect(show(n)).toEqual(["@1:1", "f", 'a"b\\c\né', "123", "true"]);
    const s = n.kind === "list" ? n.items[1] : null;
    expect(s?.kind === "atom" && s.raw).toBe('"a\\"b\\\\c\\n\\u00e9"');
  });

  it("reads one top-level atom, and drops whitespace and comments", () => {
    expect(show(node("  1brs/A  "))).toBe("1brs/A");
    const n = node("; residues of A\n(and 1brs/A ; the unit\n\t(letters C))\r\n; done");
    expect(show(n)).toEqual(["@2:1", "and", "1brs/A", ["@3:2", "letters", "C"]]);
  });

  it("gives every node its position, lines and columns from 1", () => {
    const n = node("(a\n  (b :c\n     d))");
    if (n.kind !== "list") throw new Error("a list");
    const [, b] = n.items as [QueryNode, QueryNode & { kind: "list" }];
    expect(b.items.map((x) => [x.at.line, x.at.column, x.at.offset])).toEqual([
      [2, 4, 6],
      [2, 6, 8],
      [3, 6, 16],
    ]);
  });

  it("reads deep nesting up to the limit without recursion", () => {
    const deep = `${"(f ".repeat(MAX_DEPTH)}x${")".repeat(MAX_DEPTH)}`;
    expect(readQuery(deep).issues).toEqual([]);
  });
});

describe("reader errors, as rendered", () => {
  it("seqql.read.unclosed: at the end, naming the innermost open paren", () => {
    const text = "(pick (items :of domains) :where (matches (label) kinase)";
    expect(rendered(text)).toBe('1:58 seqql.read.unclosed: unclosed "(" opened at 1:1: add ")"');
    expect(error(text).open).toEqual({ line: 1, column: 1 });
    expect(rendered("(and 1brs/A\n  (letters C)\n  (not (sel x)")).toBe(
      '3:15 seqql.read.unclosed: unclosed "(" opened at 3:3: add ")"',
    );
  });

  it("seqql.read.extra_close", () => {
    expect(rendered("(and 1brs/A (letters C)))")).toBe('1:25 seqql.read.extra_close: unexpected ")": nothing is open');
    expect(rendered(")")).toBe('1:1 seqql.read.extra_close: unexpected ")": nothing is open');
  });

  it("seqql.read.trailing", () => {
    expect(rendered("(letters C) (letters H)")).toBe(
      "1:13 seqql.read.trailing: a query is one form; found more after it at 1:13",
    );
    expect(rendered("1brs/A 1brs/B")).toBe("1:8 seqql.read.trailing: a query is one form; found more after it at 1:8");
  });

  it("seqql.read.string: unterminated (also at a line end) and invalid escapes", () => {
    expect(rendered('(matches (label) "kinase)')).toBe("1:18 seqql.read.string: unterminated string started at 1:18");
    expect(rendered('(f "a\nb")')).toBe("1:4 seqql.read.string: unterminated string started at 1:4");
    expect(rendered('(f "a\\x")')).toBe(
      "1:4 seqql.read.string: invalid string started at 1:4: JSON's escapes only (\\\" \\\\ \\n \\t \\uXXXX)",
    );
  });

  it("seqql.read.reserved: brackets, commas and $ point at sel, spaces and input", () => {
    expect(rendered("(range 1brs/A[10-20])")).toBe(
      '1:14 seqql.read.reserved: "[" is not part of the language: write (sel "1brs/A[10-20]")',
    );
    expect(rendered("(f ])")).toBe(
      '1:4 seqql.read.reserved: "]" is not part of the language: write (sel "1brs/A[10-20]")',
    );
    expect(rendered("(letters C, H)")).toBe(
      '1:11 seqql.read.reserved: "," is not part of the language: separate arguments with spaces',
    );
    expect(rendered("(onto $selection msa)")).toBe(
      '1:7 seqql.read.reserved: "$selection" is not part of the language: write (input selection)',
    );
  });

  it("seqql.read.keyword: lower-case names only, inside a list", () => {
    expect(rendered("(items :Of domains)")).toBe('1:8 seqql.read.keyword: ":Of": keywords are lower-case, as :of');
    expect(rendered("(items :1)")).toBe(
      '1:8 seqql.read.keyword: ":1" is not a keyword: a keyword is ":" and a lower-case name, as :of',
    );
    expect(rendered("(f : x)")).toContain('":" is not a keyword');
    expect(rendered(":of")).toBe(`1:1 seqql.read.keyword: ":of": a keyword goes inside a list, after a symbol's name`);
  });

  it("seqql.read.head: a list starts with an unquoted word", () => {
    for (const [text, column] of [
      ["((items) x)", 2],
      ['("items")', 2],
      ["(45 x)", 2],
      ["(:of x)", 2],
      ["(true)", 2],
    ]) {
      expect(rendered(text as string)).toBe(
        `1:${column} seqql.read.head: a list starts with a symbol's name, as (items :of domains)`,
      );
    }
  });

  it("seqql.read.empty, .depth, .size", () => {
    expect(rendered("(and ())")).toBe(
      `1:6 seqql.read.empty: "()" is empty: a list starts with a symbol's name, as (all)`,
    );
    expect(rendered("  ; nothing\n")).toBe("seqql.read.empty: the query is empty: write one form, as (all)");
    expect(rendered(`${"(f ".repeat(MAX_DEPTH + 1)}x`)).toBe(
      `1:193 seqql.read.depth: lists nest deeper than ${MAX_DEPTH}`,
    );
    expect(rendered(`(f ${"x".repeat(MAX_TEXT)})`)).toBe(
      `seqql.read.size: the query is longer than ${MAX_TEXT} characters`,
    );
  });

  it("reports the path as the root and no signature", () => {
    const i = error("(f");
    expect([i.path, i.signature, i.severity]).toEqual(["", null, "error"]);
  });
});
