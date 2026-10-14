/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The TypeScript builder (seqspec:seqql#overview): the base table, callable, by full name, as MolQL's
// `MolScriptBuilder` re-exports its tables. `Q.seq.filter.pick(q, { where: t })` builds the JSON of
// seqspec:seqql#names. `Q.query` stores an expression with its canonical text; `Q.read` reads, binds and checks a text,
// or throws.
import { checkText } from './check.js';
import { renderIssue } from './issue.js';
import type { Language } from './language.js';
import { printQuery } from './print.js';
import { BASE_SYMBOLS, type BaseSymbolName } from './table.js';
import { type Expression, type QueryIssue, SEQQL_VERSION, type SeqQLQuery, type SymbolDef } from './types.js';

export type SymbolPath<N extends string, V> = N extends `${infer Head}.${infer Tail}`
  ? { readonly [K in Head]: SymbolPath<Tail, V> }
  : { readonly [K in N]: V };
export type UnionToIntersection<U> = (U extends unknown ? (u: U) => void : never) extends (i: infer I) => void
  ? I
  : never;

/** The symbols of a set of full names as a tree by their dotted parts. */
export type SymbolTree<N extends string> = UnionToIntersection<N extends unknown ? SymbolPath<N, SymbolDef> : never>;

/** A tree of symbols by full name: `tree.seq.item.label` is `seq.item.label`. */
export function symbolTree<N extends string>(symbols: Iterable<SymbolDef & { readonly name: N }>): SymbolTree<N> {
  const root: Record<string, unknown> = {};
  for (const s of symbols) {
    const parts = s.name.split('.');
    const leaf = parts.pop() as string;
    let node = root;
    for (const part of parts) {
      node[part] ??= {};
      node = node[part] as Record<string, unknown>;
    }
    node[leaf] = s;
  }
  return root as SymbolTree<N>;
}

/** What `Q.read` throws: its message the errors as rendered, one a line (seqspec:seqql#errors). */
export type QueryError = Error & { readonly issues: readonly QueryIssue[] };

export const Q: SymbolTree<BaseSymbolName> & {
  readonly query: (expr: Expression) => SeqQLQuery;
  readonly read: (text: string, language?: Language) => SeqQLQuery;
} = Object.assign(symbolTree(BASE_SYMBOLS), {
  /** What a document or story stores: the version, the JSON and its canonical text. */
  query: (expr: Expression): SeqQLQuery => ({ seqql: SEQQL_VERSION, expr, text: printQuery(expr) }),
  /** A text read, bound and checked, as stored (macros as written); throws a `QueryError` on an error. */
  read: (text: string, language?: Language): SeqQLQuery => {
    const r = checkText(text, language);
    if (r.expr === null || r.expanded === null) {
      const errors = r.issues.filter((i) => i.severity === 'error').map(renderIssue);
      throw Object.assign(new Error(errors.join('\n')), { issues: r.issues }) as QueryError;
    }
    return { seqql: SEQQL_VERSION, expr: r.expr, text: printQuery(r.expr, language) };
  },
});
