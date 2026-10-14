/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The S-expression reader (seqspec:seqql#the-reader): text to a tree of lists, atoms and keywords with a position on
// every node; no table. Iterative, so no input can overflow its stack; the first error stops it (seqspec:seqql#errors).
// Never throws.
import { queryIssue } from "./issue.js";
import type { QueryIssue } from "./types.js";

/** 1-based line and column, and the 0-based offset into the text. */
export interface TextPosition {
  readonly line: number;
  readonly column: number;
  readonly offset: number;
}

export type QueryNode =
  /** `(head arg…)`: the head is an unquoted word (the reader checks it). */
  | { readonly kind: "list"; readonly at: TextPosition; readonly items: readonly QueryNode[] }
  /** A word, string, number or boolean; `raw` is the source text (a string's with its quotes and escapes). */
  | {
      readonly kind: "atom";
      readonly at: TextPosition;
      readonly value: string | number | boolean;
      readonly raw: string;
    }
  /** `:name`, without the colon. */
  | { readonly kind: "keyword"; readonly at: TextPosition; readonly name: string };

export interface QueryRead {
  readonly node: QueryNode | null;
  readonly issues: readonly QueryIssue[];
}

export const MAX_DEPTH = 64;
export const MAX_TEXT = 65536;
/** JSON's number grammar: a word that matches all of it is a number. */
export const NUMBER = /^-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?$/;
const DELIMITER = /[ \t\r\n()";[\],]/;
const KEYWORD = /^:[a-z][a-z0-9-]*$/;
const RESERVED: Record<string, string> = {
  "[": 'write (sel "1brs/A[10-20]")',
  "]": 'write (sel "1brs/A[10-20]")',
  ",": "separate arguments with spaces",
};

const where = (p: TextPosition) => `${p.line}:${p.column}`;

/** The start of a text in itself: line 1, column 1, offset 0. */
export const TEXT_START: TextPosition = { line: 1, column: 1, offset: 0 };

/** Reads one form. `start` is where `text` begins in a larger text (a block of a reply, seqspec:seqql#queries-in-text): every position, the
 * nodes' and the messages', is then in the larger text. */
export function readQuery(text: string, start: TextPosition = TEXT_START): QueryRead {
  const fail = (code: string, message: string, at: TextPosition | null, open: TextPosition | null = null) => ({
    node: null,
    issues: [queryIssue(code, message, { at, open })],
  });
  if (text.length > MAX_TEXT) return fail("seqql.read.size", `the query is longer than ${MAX_TEXT} characters`, null);
  const n = text.length;
  let i = 0;
  let line = start.line;
  let lineStart = 1 - start.column;
  const pos = (o: number): TextPosition => ({ line, column: o - lineStart + 1, offset: start.offset + o });
  const stack: { at: TextPosition; items: QueryNode[] }[] = [];
  let top: QueryNode | null = null;

  for (;;) {
    for (; i < n; i++) {
      const c = text[i];
      if (c === ";") while (i + 1 < n && text[i + 1] !== "\n") i++;
      else if (c === "\n") {
        line++;
        lineStart = i + 1;
      } else if (c !== " " && c !== "\t" && c !== "\r") break;
    }
    if (i >= n) break;
    const c = text[i] as string;
    const at = pos(i);
    // after the one form: a ")" closes nothing (below), anything else is more
    if (top !== null && c !== ")") {
      return fail("seqql.read.trailing", `a query is one form; found more after it at ${where(at)}`, at);
    }
    let node: QueryNode;
    if (c === "(") {
      if (stack.length >= MAX_DEPTH) return fail("seqql.read.depth", `lists nest deeper than ${MAX_DEPTH}`, at);
      stack.push({ at, items: [] });
      i++;
      continue;
    } else if (c === ")") {
      const open = stack.pop();
      if (!open) return fail("seqql.read.extra_close", 'unexpected ")": nothing is open', at);
      i++;
      const head = open.items[0];
      if (!head)
        return fail("seqql.read.empty", '"()" is empty: a list starts with a symbol\'s name, as (all)', open.at);
      if (head.kind !== "atom" || typeof head.value !== "string" || head.raw.startsWith('"')) {
        const message = "a list starts with a symbol's name, as (items :of domains)";
        return fail("seqql.read.head", message, head.at, open.at);
      }
      node = { kind: "list", at: open.at, items: open.items };
    } else if (c === '"') {
      let j = i + 1;
      while (j < n && text[j] !== '"' && text[j] !== "\n") j += text[j] === "\\" ? 2 : 1;
      if (j >= n || text[j] !== '"')
        return fail("seqql.read.string", `unterminated string started at ${where(at)}`, at);
      const raw = text.slice(i, j + 1);
      let value: string;
      try {
        value = JSON.parse(raw) as string;
      } catch {
        const message = `invalid string started at ${where(at)}: JSON's escapes only (\\" \\\\ \\n \\t \\uXXXX)`;
        return fail("seqql.read.string", message, at);
      }
      node = { kind: "atom", at, value, raw };
      i = j + 1;
    } else if (RESERVED[c]) {
      return fail("seqql.read.reserved", `"${c}" is not part of the language: ${RESERVED[c]}`, at);
    } else {
      let j = i;
      while (j < n && !DELIMITER.test(text[j] as string)) j++;
      const raw = text.slice(i, j);
      i = j;
      if (raw.startsWith("$")) {
        const message = `"${raw}" is not part of the language: write (input ${raw.slice(1) || "selection"})`;
        return fail("seqql.read.reserved", message, at);
      } else if (raw.startsWith(":")) {
        if (!KEYWORD.test(raw)) {
          const lower = raw.toLowerCase();
          const message = KEYWORD.test(lower)
            ? `"${raw}": keywords are lower-case, as ${lower}`
            : `"${raw}" is not a keyword: a keyword is ":" and a lower-case name, as :of`;
          return fail("seqql.read.keyword", message, at);
        }
        if (stack.length === 0) {
          return fail("seqql.read.keyword", `"${raw}": a keyword goes inside a list, after a symbol's name`, at);
        }
        node = { kind: "keyword", at, name: raw.slice(1) };
      } else {
        const num = NUMBER.test(raw) ? Number(raw) : Number.NaN;
        const value = raw === "true" ? true : raw === "false" ? false : Number.isFinite(num) ? num : raw;
        node = { kind: "atom", at, value, raw };
      }
    }
    if (stack.length > 0) stack[stack.length - 1]?.items.push(node);
    else top = node;
  }
  const open = stack[stack.length - 1];
  if (open) return fail("seqql.read.unclosed", `unclosed "(" opened at ${where(open.at)}: add ")"`, pos(n), open.at);
  if (top === null) return fail("seqql.read.empty", "the query is empty: write one form, as (all)", null);
  return { node: top, issues: [] };
}
