/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Query issues (seqspec:seqql#errors): making them, rendering them for a person or the agent, and the did-you-mean of
// unknown names.
import { seqqlSeverity } from "./codes.ts";
import type { QueryIssue } from "./types.ts";

export type QueryIssuePlace = Partial<Pick<QueryIssue, "path" | "at" | "open" | "signature" | "known">>;

/** An issue, its severity the code's registered one (`SeqqlIssueCodes`). */
export function queryIssue(code: string, message: string, place: QueryIssuePlace = {}): QueryIssue {
  const { path = "", at = null, open = null, signature = null, known } = place;
  const p = (x: { line: number; column: number } | null) => x && { line: x.line, column: x.column };
  const severity = seqqlSeverity(code);
  return { code, severity, path, message, at: p(at), open: p(open), signature, ...(known ? { known } : {}) };
}

/** One line, `line:column code: message`, then the signature, indented, when there is one. */
export function renderIssue(issue: QueryIssue): string {
  const at = issue.at ? `${issue.at.line}:${issue.at.column} ` : "";
  return `${at}${issue.code}: ${issue.message}${issue.signature ? `\n  ${issue.signature}` : ""}`;
}

/** Optimal-string-alignment distance: Damerau–Levenshtein with adjacent swaps. */
function distance(a: string, b: string): number {
  const d = Array.from({ length: a.length + 1 }, (_, i) => Array.from({ length: b.length + 1 }, (_, j) => i + j * +!i));
  const at = (i: number, j: number) => d[i]?.[j] as number;
  for (let i = 1; i <= a.length; i++) {
    for (let j = 1; j <= b.length; j++) {
      const swap = a[i - 1] === b[j - 2] && a[i - 2] === b[j - 1] ? at(i - 2, j - 2) + 1 : Infinity;
      (d[i] as number[])[j] = Math.min(
        at(i - 1, j) + 1,
        at(i, j - 1) + 1,
        at(i - 1, j - 1) + +(a[i - 1] !== b[j - 1]),
        swap,
      );
    }
  }
  return at(a.length, b.length);
}

/** The names within 2 edits of `word` (3 for words over 8 characters), nearest first, at most 8. */
export function suggest(word: string, names: Iterable<string>): string[] {
  const max = word.length > 8 ? 3 : 2;
  const near: [string, number][] = [];
  for (const name of new Set(names)) {
    const k = Math.abs(word.length - name.length) > max ? max + 1 : distance(word, name);
    if (k <= max && name !== word) near.push([name, k]);
  }
  return near
    .sort((x, y) => x[1] - y[1])
    .slice(0, 8)
    .map(([name]) => name);
}
