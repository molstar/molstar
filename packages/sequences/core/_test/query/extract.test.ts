/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Queries in replies (seqspec:seqql#queries-in-text): [SeqQL] blocks found with their place in the whole text, issues
// inside a block placed in the reply, the edge cases as instructive issues, ```seqql fences, and the whole text as the
// fallback.
import { describe, expect, it } from "vitest";
import { checkText, extractQueries, type QueryBlock, readQuery, renderIssue } from "../../src/query/index.js";

const rendered = (text: string) => extractQueries(text).issues.map(renderIssue);
/** Each block's issues as rendered, placed in the reply. */
const checked = (text: string) =>
  extractQueries(text).blocks.map((b) => checkText(b.text, undefined, b.at).issues.map(renderIssue));

describe("extractQueries", () => {
  it("finds a block in prose, with its offset, line and column in the whole text", () => {
    const reply = "There are three conservation tracks;\nI chose cons: [SeqQL](values :of cons :above 0.8)[/SeqQL].";
    const { blocks, issues } = extractQueries(reply);
    expect(issues).toEqual([]);
    expect(blocks).toEqual([
      { text: "(values :of cons :above 0.8)", at: { line: 2, column: 22, offset: 58 }, via: "tag" },
    ]);
    expect(reply.slice(blocks[0]?.at.offset)).toMatch(/^\(values/);
    expect(checked(reply)).toEqual([[]]);
  });

  it("returns several blocks in order, the tags matched in any case", () => {
    const reply = "Chain A's side: [SeqQL]1brs/A[/SeqQL]\nand D's: [seqql] (and 1brs/D (letters C)) [/SEQQL]";
    const { blocks, issues } = extractQueries(reply);
    expect(issues).toEqual([]);
    expect(blocks.map((b) => [b.text, b.at.line, b.at.column])).toEqual([
      ["1brs/A", 1, 24],
      [" (and 1brs/D (letters C)) ", 2, 17],
    ]);
  });

  it("places the reader's, binder's, checker's and sel's issues in the reply, messages' positions too", () => {
    const reply = [
      "Here:",
      "[SeqQL](items :of domains",
      '  :where (> (lenght) 50))[/SeqQL] and [SeqQL](onto (sel "1brs/A[10-20") msa :kind codon)[/SeqQL]',
      "[SeqQL]",
      "(and 1brs/A (values :of cons",
      "[/SeqQL]",
    ].join("\n");
    expect(checked(reply)).toEqual([
      ['3:14 seqql.symbol.unknown: unknown symbol "lenght"; did you mean "length"?\n  (length)'],
      [
        '3:70 seqql.sel.syntax: expected "]" to close "[" at 3:64\n  (sel S)',
        "3:83 seqql.type: :kind is a kind of sequence (unit, molecule, row, columns), and codon is not one: across a codon mapping leave :kind out, (onto Q gene) gives the codons\n  (onto Q R :kind K :keep b)",
      ],
      ['6:1 seqql.read.unclosed: unclosed "(" opened at 5:13: add ")"'],
    ]);
  });

  it("reads a block that is the whole text from 1:1, as checkText alone does", () => {
    const text = '(sel "1brs/A[10-20")';
    const [block] = extractQueries(`[SeqQL]${text}[/SeqQL]`).blocks as [QueryBlock];
    expect(checkText(text).issues.map((i) => i.at)).toEqual([{ line: 1, column: 19 }]);
    expect(checkText(block.text, undefined, block.at).issues.map((i) => i.at)).toEqual([{ line: 1, column: 26 }]);
    expect(readQuery("(all)", { line: 4, column: 3, offset: 40 }).node?.at).toEqual({ line: 4, column: 3, offset: 40 });
  });

  it("reports an unclosed tag at the end, with the tag it opened at", () => {
    const reply = 'First [SeqQL](all)[/SeqQL], then\n[SeqQL](sel "p53@0:50")';
    const { blocks, issues } = extractQueries(reply);
    expect(blocks.map((b) => b.text)).toEqual(["(all)"]);
    expect(issues.map(renderIssue)).toEqual([
      "2:24 seqql.reply.unclosed: unclosed [SeqQL] opened at 2:1: end the query with [/SeqQL]",
    ]);
    expect(issues[0]?.open).toEqual({ line: 2, column: 1 });
  });

  it("reports a nested tag, reading the inner block", () => {
    const { blocks, issues } = extractQueries("[SeqQL] the query: [SeqQL](all)[/SeqQL]");
    expect(blocks.map((b) => b.text)).toEqual(["(all)"]);
    expect(issues.map(renderIssue)).toEqual([
      "1:20 seqql.reply.nested: [SeqQL] inside the block opened at 1:1: close each block with [/SeqQL] before the next",
    ]);
  });

  it("reports an empty block and a close with no open", () => {
    expect(rendered("[SeqQL]  \n [/SeqQL] [/SeqQL]")).toEqual([
      "1:1 seqql.reply.empty: an empty block: put one query in it, as [SeqQL](all)[/SeqQL]",
      "2:11 seqql.reply.unopened: [/SeqQL] closes no block: start each query with [SeqQL]",
    ]);
  });

  it("falls back to the whole text only when it reads as one form", () => {
    expect(extractQueries("  (and p53 (letters C))\n").blocks).toEqual([
      { text: "  (and p53 (letters C))\n", at: { line: 1, column: 1, offset: 0 }, via: "text" },
    ]);
    // a query that reads but does not check is still the query: its issues are the checker's
    expect(extractQueries("(itemz :of domains)").blocks.map((b) => b.via)).toEqual(["text"]);
    expect(rendered("The cysteines: (and p53 (letters C))")).toEqual([
      "seqql.reply.none: no query found: wrap each query in [SeqQL]…[/SeqQL]",
    ]);
    expect(rendered("(and p53")).toEqual(["seqql.reply.none: no query found: wrap each query in [SeqQL]…[/SeqQL]"]);
    expect(rendered(" \n")).toEqual(["seqql.reply.none: the reply is empty: wrap each query in [SeqQL]…[/SeqQL]"]);
  });

  it("accepts ```seqql fences in a text without tags; tags win when there are both", () => {
    const fenced = "Two:\n```seqql\n(all)\n```\nand\n```SeqQL\n1brs/A\n```\n```lisp\n(empty)\n```";
    expect(extractQueries(fenced)).toEqual({
      blocks: [
        { text: "(all)\n", at: { line: 3, column: 1, offset: 14 }, via: "fence" },
        { text: "1brs/A\n", at: { line: 7, column: 1, offset: 37 }, via: "fence" },
      ],
      issues: [],
    });
    const both = "```seqql\n[SeqQL](all)[/SeqQL]\n```\n```seqql\n(empty)\n```";
    expect(extractQueries(both).blocks.map((b) => [b.text, b.via])).toEqual([["(all)", "tag"]]);
    expect(rendered("```seqql\n\n```")).toEqual([
      "1:1 seqql.reply.empty: an empty block: put one query in it, as [SeqQL](all)[/SeqQL]",
    ]);
  });
});
