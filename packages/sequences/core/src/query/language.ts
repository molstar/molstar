/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A host's language: the base table plus its extensions' symbols, assembled per call with no module state
// (seqspec:seqql#extension), and the alias index the binder resolves short names from
// (seqspec:seqql#binding-and-checking).
import { queryIssue } from './issue.js';
import { BASE_SYMBOLS } from './table.js';
import type { QueryIssue, SymbolDef } from './types.js';

export interface Language {
  /** By full name. */
  readonly symbols: ReadonlyMap<string, SymbolDef>;
  /** By short name: the binder's candidates, in table order (the base, then the extensions in load order). */
  readonly aliases: ReadonlyMap<string, readonly SymbolDef[]>;
}

/** What `createLanguage` reads of an extension (its registry entry: `QueryExtension` in `@molstar/seq-core/model/query`). */
export interface LanguageExtension {
  readonly symbols?: readonly SymbolDef[];
}

/** Namespaces of the base symbols and the format's codes; an extension's full names use its own. */
export const RESERVED_NAMESPACES: ReadonlySet<string> = new Set(['core', 'seq', 'svs', 'seqv']);

const BASE_ALIASES: ReadonlySet<string> = new Set(BASE_SYMBOLS.flatMap((s) => (s.alias === null ? [] : [s.alias])));
const NUMBER = /^-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?$/;
// A reader word (seqspec:seqql#the-reader) with no dot, so an alias never reads as a full name.
const WORD = /^[^\s()";[\],:$.][^\s()";[\],.]*$/;

/** Why an extension's alias is refused, or null: a base short name, or a word the reader reads as something else. */
function aliasProblem(alias: string): string | null {
  if (BASE_ALIASES.has(alias)) return 'is a base short name';
  if (!WORD.test(alias) || NUMBER.test(alias) || alias === 'true' || alias === 'false') {
    return 'is not a word the reader reads as a name (no spaces, parentheses, brackets, quotes, commas, dots or a leading : or $; not a number, true or false)';
  }
  return null;
}

/** Whether the binder can tell two symbols with one short name apart: only properties reading different slots. */
const distinct = (a: SymbolDef, b: SymbolDef) => a.kind === 'property' && b.kind === 'property' && a.reads !== b.reads;

const issue = (code: string, message: string): QueryIssue => queryIssue(code, message);

/**
 * The language of the base table and the given extensions. Refused, each with an issue: a full name outside an
 * extension namespace (`seqql.symbol.namespace`) or defined twice (`seqql.symbol.duplicate`, the later entry left
 * out); an alias that is reserved (`seqql.alias.reserved`) or shared by symbols the binder cannot tell apart
 * (`seqql.alias.collision`: none of them gets it, whatever the load order). A refused alias leaves its symbol
 * reachable by full name.
 */
export function createLanguage(ext: readonly LanguageExtension[] = []): {
  language: Language;
  issues: readonly QueryIssue[];
} {
  const issues: QueryIssue[] = [];
  const symbols = new Map<string, SymbolDef>();
  for (const s of BASE_SYMBOLS) symbols.set(s.name, s);
  const added: SymbolDef[] = [];
  for (const entry of ext) {
    for (const s of entry.symbols ?? []) {
      const dot = s.name.indexOf('.');
      if (dot <= 0 || dot === s.name.length - 1 || RESERVED_NAMESPACES.has(s.name.slice(0, dot))) {
        const msg = `"${s.name}" is not an extension's full name: write <namespace>.<name>, outside core, seq, svs and seqv`;
        issues.push(issue('seqql.symbol.namespace', msg));
      } else if (symbols.has(s.name)) {
        issues.push(issue('seqql.symbol.duplicate', `"${s.name}" is defined twice; the later definition is left out`));
      } else {
        symbols.set(s.name, s);
        added.push(s);
      }
    }
  }

  const refused = new Set<SymbolDef>();
  const byAlias = new Map<string, SymbolDef[]>();
  for (const s of added) {
    if (s.alias === null) continue;
    const problem = aliasProblem(s.alias);
    if (problem) {
      issues.push(
        issue('seqql.alias.reserved', `the short name "${s.alias}" of ${s.name} ${problem}; use the full name`),
      );
      refused.add(s);
    } else byAlias.set(s.alias, [...(byAlias.get(s.alias) ?? []), s]);
  }
  for (const [alias, list] of byAlias) {
    const clash = list.filter((s) => list.some((t) => t !== s && !distinct(s, t)));
    if (clash.length === 0) continue;
    for (const s of clash) refused.add(s);
    const names = clash.map((s) => s.name).sort();
    const who = `${names.slice(0, -1).join(', ')} and ${names[names.length - 1]}`;
    const msg = `the short name "${alias}" is defined by ${who}, which the binder cannot tell apart; use the full names`;
    issues.push(issue('seqql.alias.collision', msg));
  }

  const aliases = new Map<string, SymbolDef[]>();
  for (const s of symbols.values()) {
    if (s.alias === null || refused.has(s)) continue;
    const list = aliases.get(s.alias);
    if (list) list.push(s);
    else aliases.set(s.alias, [s]);
  }
  return { language: { symbols, aliases }, issues };
}
