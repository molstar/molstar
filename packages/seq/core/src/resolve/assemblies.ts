/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An assembly's entry (_spec/README.md "Resolve", Assemblies): its sequences derived, a system each, as a complex's
// units are: from its inline `sequences`, then its source's coverage (names, lengths, aliases), ordered, kept by
// `include`, a `sequence` child overriding the one its `name` is. Linear in the sequences, each once.

import { type Issue, issueOf, isColumn, type Node, printAddress, type Value } from '@molstar/core/repr/tree';
import type { SystemInfo } from '@molstar/seqx-core/data';
import {
  type AssemblyNames,
  assemblyNames,
  type IncludeValue,
  includeFilter,
  lookupName,
  type NameMatch,
  naturalOrder,
} from './assembly-names.js';
import { ResolveIssues } from './issues.js';
import { sourcePlan } from './model-entries.js';
import type { AssemblyEntry, AssemblySequence } from './genome.js';
import type { Placed } from './model-molecules.js';
import type { SystemEntry } from './model.js';

/** The address of an assembly's derived sequence: always the form, never an id step (_spec/README.md "Addressing"). */
export function assemblySequenceAddress(assembly: string, name: string): string {
  return printAddress([...assembly.split('/').filter(Boolean), { kind: 'sequence', params: { name } }]);
}

/** The system a name is in an assembly (its names, aliases, the `chr` prefix), or undefined. */
export function assemblySystem(entry: AssemblyEntry, name: string): string | undefined {
  const i = entry.find(name);
  return i === undefined ? undefined : entry.systemAt[i];
}

const p = (n: Node, key: string): Value | undefined => n.params?.[key];
const strs = (v: Value | undefined): string[] => (Array.isArray(v) ? v.filter((x) => typeof x === 'string') : []);

/**
 * An assembly's entry and its derived systems. `coverage` is its source's systems once open (undefined while it is
 * not, or with no source); `table` the alias source's. A sequence with no length makes it unavailable
 * (`seq.resolve.assembly.lengths`); an alias naming two is dropped (`seq.resolve.assembly.alias`).
 */
export function assemblyEntry(
  at: Placed,
  source: string | undefined,
  aliasSource: string | undefined,
  coverage: readonly SystemInfo[] | undefined,
  table: readonly SystemInfo[] | undefined,
): { entry: AssemblyEntry; systems: SystemEntry[]; issues: Issue[] } {
  const n = at.node;
  const issues: Issue[] = [];
  const say = (code: keyof typeof ResolveIssues.codes, message: string, key: string) =>
    issues.push(
      issueOf(ResolveIssues, code, message, {
        address: at.address,
        path: ['params', key],
        hint: ResolveIssues.codes[code].hint,
      }),
    );
  const match = (p(n, 'match') ?? 'chr') as NameMatch;
  const own = [p(n, 'name'), ...strs(p(n, 'aliases'))].filter((x): x is string => typeof x === 'string');
  const pending = source !== undefined && coverage === undefined;
  let sequences: AssemblySequence[] = [];
  const rows = p(n, 'sequences') as Record<string, Value> | undefined;
  if (rows) {
    const name = strs(rows.name);
    const length = isColumn(rows.length) ? (rows.length as ArrayLike<number>) : undefined;
    const aliases = Array.isArray(rows.aliases) ? rows.aliases : [];
    const cyclic = Array.isArray(rows.cyclic) ? rows.cyclic : [];
    for (let i = 0; i < name.length; i++)
      sequences.push({ name: name[i], length: length?.[i], aliases: strs(aliases[i]), cyclic: cyclic[i] === true });
  }
  // the source's sequences join the inline ones through the names rule, else follow them
  if (coverage?.length) {
    const inline = assemblyNames(sequences, match);
    for (const s of coverage) {
      const i = lookupName(inline, s.name);
      if (i === undefined)
        sequences.push({ name: s.name, length: s.length, aliases: s.aliases, sourceName: s.name, cyclic: s.cyclic });
      else sequences[i] = { ...sequences[i], length: sequences[i].length ?? s.length, sourceName: s.name };
    }
  }
  // one sequence a name: a later one of a name already given is dropped and said, never written over the first
  const seen = new Set<string>();
  const twice: string[] = [];
  sequences = sequences.filter((s) => (seen.has(s.name) ? (twice.push(s.name), false) : (seen.add(s.name), true)));
  if (twice.length)
    say(
      'seq.resolve.assembly.duplicate',
      `Sequences named twice: ${twice.slice(0, 3).join(', ')}`,
      rows ? 'sequences' : 'source',
    );
  const empty = (lookup: AssemblyNames, order: number[] = []) => ({
    entry: {
      address: at.address,
      names: own,
      ...(source ? { source } : {}),
      ...(aliasSource ? { aliasSource } : {}),
      pending,
      sequences,
      order,
      systems: [],
      systemAt: [],
      lookup,
      find: (name: string) => lookupName(lookup, name),
      include: (v: unknown) => includeFilter(v as IncludeValue | undefined, lookup),
    },
    systems: [],
    issues,
  });
  if (pending) return empty(assemblyNames([], match));
  // every sequence has a length: one without is never drawn with a guessed one
  const unsized = sequences.filter((s) => s.length === undefined);
  if (unsized.length) {
    const first = unsized.slice(0, 3).map((s) => s.name);
    say(
      'seq.resolve.assembly.lengths',
      `${unsized.length} sequences have no length (${first.join(', ')}${unsized.length > 3 ? ', …' : ''})`,
      rows ? 'sequences' : 'source',
    );
    return empty(assemblyNames([], match));
  }
  const lookup = assemblyNames(sequences, match, table);
  if (lookup.dropped.length)
    say('seq.resolve.assembly.alias', `Aliases dropped: ${lookup.dropped.slice(0, 3).join(', ')}`, 'sequences');
  const order = orderOf(p(n, 'order'), sequences, lookup);
  const keep = includeFilter(p(n, 'include') as IncludeValue | undefined, lookup);
  // the overrides: a `sequence` child by its `name` (its id when absent)
  const overrides = new Map<number, Node>();
  for (const c of n.children ?? []) {
    if (c.kind !== 'sequence') continue;
    const name = p(c, 'name') ?? c.id;
    const i = typeof name === 'string' ? lookupName(lookup, name) : undefined;
    if (i !== undefined) overrides.set(i, c);
  }
  const alphabet = (p(n, 'alphabet') ?? 'dna') as SystemEntry['alphabet'];
  const systems: SystemEntry[] = [];
  const systemAt: (string | undefined)[] = new Array(sequences.length);
  for (const i of order) {
    const s = sequences[i];
    if (keep && !keep(i, s)) continue;
    const o = overrides.get(i);
    const address = assemblySequenceAddress(at.address, s.name);
    const sourceName = s.sourceName ?? s.name;
    const label = o && p(o, 'label');
    systemAt[i] = address;
    systems.push({
      address,
      kind: 'sequence',
      alphabet: (o && (p(o, 'alphabet') as SystemEntry['alphabet'])) || alphabet,
      length: s.length,
      cyclic: o ? p(o, 'cyclic') === true : s.cyclic === true,
      name: sourceName,
      // matched against a source's coverage by either name (`isSystemInfo`)
      aliases: [...(sourceName === s.name ? [] : [s.name]), ...(s.aliases ?? [])],
      numberings: ['index', 'ordinal'],
      ...(source ? { source } : {}),
      ...(typeof label === 'string' ? { label } : {}),
      access: sourcePlan(source ?? '', sourceName, 'letters'),
    });
  }
  const entry: AssemblyEntry = {
    ...empty(lookup, order).entry,
    systems: systems.map((s) => s.address),
    systemAt,
  };
  return { entry, systems, issues };
}

/** An assembly's `order`: `natural` (absent), `length` (longest first), `source` (as listed), or names first. */
function orderOf(order: Value | undefined, sequences: readonly AssemblySequence[], lookup: AssemblyNames): number[] {
  const all = sequences.map((_, i) => i);
  if (order === 'source') return all;
  if (order === 'length') return all.sort((a, b) => (sequences[b].length ?? 0) - (sequences[a].length ?? 0) || a - b);
  const natural = naturalOrder(sequences.map((s) => s.name));
  if (!Array.isArray(order)) return natural;
  const first: number[] = [];
  const seen = new Set<number>();
  for (const name of strs(order)) {
    const i = lookupName(lookup, name);
    if (i !== undefined && !seen.has(i)) seen.add(i), first.push(i);
  }
  return [...first, ...natural.filter((i) => !seen.has(i))];
}
