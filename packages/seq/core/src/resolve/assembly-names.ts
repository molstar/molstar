/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An assembly's names rule, order and include (_spec/README.md "Resolve", Assemblies): pure over a list of sequences,
// built once per coverage. A name is looked up exactly, then through the aliases (an alias table's, then each
// sequence's), then with `match: chr` with a leading `chr` added or dropped; an exact name always wins, an alias that
// names two sequences (or joins two of different lengths) is dropped and said, never guessed.

import { includePattern } from '@molstar/seqx-core/lang';

/** A sequence as the names rule reads it: its name, its length when known, its own aliases. */
export interface NamedSequence {
  readonly name: string;
  readonly length?: number;
  readonly aliases?: readonly string[];
}

export type NameMatch = 'chr' | 'exact';

/** The names of an assembly: each name and kept alias to its sequence's index; the aliases dropped. */
export interface AssemblyNames {
  readonly index: ReadonlyMap<string, number>;
  readonly match: NameMatch;
  /** The aliases dropped (naming two sequences, or two of different lengths), in the order met. */
  readonly dropped: readonly string[];
}

/**
 * The names index of `sequences`: their names first (the first of a name twice wins), then the alias table's
 * (`table`: a name and its aliases a row), then each sequence's own. Why aliases after names: an exact name always
 * wins, so `MT` stays the sequence named `MT` whatever alias says otherwise.
 */
export function assemblyNames(
  sequences: readonly NamedSequence[],
  match: NameMatch = 'chr',
  table: readonly NamedSequence[] = [],
): AssemblyNames {
  const index = new Map<string, number>();
  for (let i = 0; i < sequences.length; i++) if (!index.has(sequences[i].name)) index.set(sequences[i].name, i);
  const names = new Set(index.keys());
  // an alias's first claim, and whether a second claim made it ambiguous
  const claims = new Map<string, number>();
  const bad = new Set<string>();
  const claim = (alias: string, i: number) => {
    if (alias === sequences[i].name) return;
    const named = names.has(alias) ? index.get(alias)! : undefined;
    if (named !== undefined) {
      // a name of another sequence: the name wins; an alias joining two of different known lengths is said
      const a = sequences[named].length;
      const b = sequences[i].length;
      if (named !== i && a !== undefined && b !== undefined && a !== b) bad.add(alias);
      return;
    }
    const was = claims.get(alias);
    if (was === undefined) claims.set(alias, i);
    else if (was !== i) bad.add(alias);
  };
  for (const row of table) {
    const i = lookup(index, row.name, match);
    if (i === undefined) continue;
    for (const a of row.aliases ?? []) claim(a, i);
  }
  for (let i = 0; i < sequences.length; i++) for (const a of sequences[i].aliases ?? []) claim(a, i);
  for (const [alias, i] of claims) if (!bad.has(alias)) index.set(alias, i);
  return { index, match, dropped: [...bad] };
}

/** The index of the sequence a name is (exact, alias, then the `chr` prefix added or dropped), or undefined. */
export function lookupName(names: AssemblyNames, name: string): number | undefined {
  return lookup(names.index, name, names.match);
}

function lookup(index: ReadonlyMap<string, number>, name: string, match: NameMatch): number | undefined {
  const exact = index.get(name);
  if (exact !== undefined || match === 'exact') return exact;
  // `chr1` is `1`, `M` is `chrM`; `MT` is never `M` (only the prefix moves)
  const bare = /^chr/i.test(name) ? name.slice(3) : undefined;
  return bare ? index.get(bare) : index.get(`chr${name}`);
}

/**
 * The order of `natural` (_spec/README.md "Resolve", Assemblies): on the name less a leading `chr` (any case), whole
 * numbers ascending, then X, Y, W, Z, then M, MT, then the rest in the source's order. Stable.
 */
export function naturalOrder(names: readonly string[]): number[] {
  const rank = (name: string): [number, number] => {
    const n = name.replace(/^chr/i, '');
    if (/^\d+$/.test(n)) return [0, Number(n)];
    const sex = ['X', 'Y', 'W', 'Z'].indexOf(n.toUpperCase());
    if (sex >= 0) return [1, sex];
    const m = ['M', 'MT'].indexOf(n.toUpperCase());
    return m >= 0 ? [2, m] : [3, 0];
  };
  const ranks = names.map(rank);
  return names.map((_, i) => i).sort((a, b) => ranks[a][0] - ranks[b][0] || ranks[a][1] - ranks[b][1] || a - b);
}

/** An `include` as the language writes it: names and patterns (`^…`), or `{ names?, min_length? }`. */
export type IncludeValue = readonly string[] | { readonly names?: readonly string[]; readonly min_length?: number };

/**
 * Whether a sequence is kept by `include` (absent: every one): a name through the names rule, or a pattern (a string
 * starting with `^`, `includePattern`'s language: matched by a table, never a backtracking search; one it does not take
 * is the language check's `seq.lang.include.pattern`, read here as a name), and a length floor.
 */
export function includeFilter(
  include: IncludeValue | undefined,
  names: AssemblyNames,
): ((i: number, s: NamedSequence) => boolean) | undefined {
  if (include === undefined) return undefined;
  const list = Array.isArray(include) ? include : ((include as { names?: readonly string[] }).names ?? undefined);
  const floor = Array.isArray(include) ? undefined : (include as { min_length?: number }).min_length;
  const ids = new Set<number>();
  const patterns: ((name: string) => boolean)[] = [];
  for (const x of list ?? []) {
    const re = x.startsWith('^') ? includePattern(x) : undefined;
    if (re) patterns.push(re);
    else {
      const i = lookupName(names, x);
      if (i !== undefined) ids.add(i);
    }
  }
  return (i, s) =>
    (list === undefined || ids.has(i) || patterns.some((re) => re(s.name))) &&
    (floor === undefined || (s.length ?? 0) >= floor);
}
