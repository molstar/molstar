/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Data over an assembly, read by sequence (_spec/README.md "Resolve", Assemblies: data over an assembly), in the genome
// chunk: a `ranges`, `values`, `pairs`, `grid` or `matrix` node whose `target` (`a`, `x`) names an assembly gets a plan
// for each derived sequence, on the sequence's own system. Inline rows are split by their `target` field through the
// names rule (one inline source over every sequence named, its systems named by their addresses); a source is read
// by region on each sequence under the name it lists it by (its coverage's names mapped by the names rule both ways,
// `17` is `chr17`; before it opens, the sequence's printed name). A sequence `include` leaves out reads nothing.

import { exactAddress, type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import {
  inlineSeq,
  type SeqData,
  type SeqPair,
  SeqRuntime,
  type Shape,
  type ShapeInfo,
  type SystemInfo,
  seqToValue,
} from '@molstar/seqx-core/data';
import { channelAddress } from '@molstar/seqx-core/data/channels';
import type { AssemblyEntry } from './genome.js';
import { GenomeIssues } from './genome-issues.js';
import { pick } from './genome-model.js';
import { inlineRanges, sourcePlan } from './model-entries.js';
import type { Placed } from './model-molecules.js';
import { statesSystemOf, StatesSystemName } from './model-region.js';
import type { AccessPlan, AnnotationEntry, SystemEntry } from './model.js';

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/** What reading data over an assembly reads of the build. */
export interface AssemblyDataInput {
  readonly assemblies: ReadonlyMap<string, AssemblyEntry>;
  readonly systems: ReadonlyMap<string, SystemEntry>;
  /** The node's source's address, when it reads one. */
  readonly source?: string;
  /** The systems that source listed once open (`listedSystems`), else undefined. */
  readonly listed?: readonly SystemInfo[];
  /** The states a matrix source listed, once open. */
  readonly states?: number;
}

/**
 * A data node over an assembly as its annotation entry (`bySystem`: a plan per derived sequence), its states system for
 * a matrix, and what it says; undefined when its assembly is not one of the model's.
 */
export function assemblyAnnotation(
  at: Placed,
  input: AssemblyDataInput,
): { annotation: AnnotationEntry; states?: SystemEntry; issues: Issue[] } | undefined {
  const kind = at.node.kind as AnnotationEntry['shape'];
  const firstKey = kind === 'pairs' ? 'a' : kind === 'grid' ? 'x' : 'target';
  const secondKey = kind === 'pairs' ? 'b' : kind === 'grid' ? 'y' : undefined;
  const a = input.assemblies.get(exactAddress(p(at.node, firstKey)) ?? '');
  if (!a) return undefined;
  const second = secondKey === undefined ? undefined : exactAddress(p(at.node, secondKey));
  const b = second === undefined ? undefined : input.assemblies.get(second);
  const issues: Issue[] = [];
  const base = {
    address: at.address,
    shape: kind,
    target: a.address,
    ...(kind === 'pairs' || kind === 'grid' ? { ySystem: second ?? a.address } : {}),
  };
  if (input.source === undefined) {
    // inline: ranges split by their `target` field; the other kinds name no sequence a row (nothing read)
    if (kind !== 'ranges') return { annotation: { ...base, fields: [], bySystem: new Map() }, issues };
    return { annotation: inlineOverAssembly(at, a, base, input.systems), issues };
  }
  const field = p(at.node, 'field');
  const read = typeof field === 'string' && kind !== 'ranges' ? channelAddress(input.source, field) : input.source;
  const names = sourceNames(a, input.listed);
  const bySystem = new Map<string, AccessPlan>();
  const add = (e: AssemblyEntry, n: ReadonlyMap<number, string>) => {
    for (let i = 0; i < e.sequences.length; i++) {
      const system = e.systemAt[i];
      if (system !== undefined) bySystem.set(system, sourcePlan(read, n.get(i) ?? e.sequences[i].name, kind as Shape));
    }
  };
  add(a, names.byIndex);
  // a second assembly's sequences (synteny between two genomes) by its own names rule
  if (b && b !== a) add(b, sourceNames(b, input.listed).byIndex);
  if (names.unknown.length && !(b && b !== a))
    issues.push(
      issueOf(
        GenomeIssues,
        'seq.resolve.assembly.name',
        `The source names ${names.unknown.length} sequences the assembly lacks (${names.unknown.slice(0, 3).join(', ')})`,
        { address: at.address, hint: GenomeIssues.codes['seq.resolve.assembly.name'].hint },
      ),
    );
  const first = bySystem.values().next().value;
  const plan = first ?? sourcePlan(read, '', kind as Shape);
  const fields = kind === 'grid' ? ['value'] : [];
  if (kind !== 'matrix') return { annotation: { ...base, fields, access: plan, bySystem }, issues };
  // a matrix's states: a system of their own, as many as its source lists (one set over every sequence)
  const statesAddress = statesSystemOf(at.address);
  const states: SystemEntry = {
    address: statesAddress,
    kind: 'sequence',
    ...(input.states !== undefined ? { length: input.states } : {}),
    cyclic: false,
    name: StatesSystemName,
    numberings: ['index', 'ordinal'],
    access: sourcePlan(input.source, StatesSystemName, 'matrix'),
  };
  return { annotation: { ...base, fields, states: statesAddress, access: plan, bySystem }, states, issues };
}

/**
 * The names a source lists each sequence of `a` by (its index → the listed name), through the names rule both ways
 * (a listed name or alias found by `find`, the first listed wins), and the listed names no sequence is (in order).
 */
function sourceNames(
  a: AssemblyEntry,
  listed: readonly SystemInfo[] | undefined,
): { byIndex: ReadonlyMap<number, string>; unknown: readonly string[] } {
  const byIndex = new Map<number, string>();
  const unknown: string[] = [];
  for (const s of listed ?? []) {
    // a sided source's systems are its sides' (sided pairs): the names rule reads the name alone
    let i = a.find(s.name);
    for (const alias of s.aliases ?? []) if (i === undefined) i = a.find(alias);
    if (i === undefined) unknown.push(s.name);
    else if (!byIndex.has(i)) byIndex.set(i, s.name);
  }
  return { byIndex, unknown: [...new Set(unknown)] };
}

/**
 * Inline ranges over an assembly: the rows split by their `target` field through the names rule, each sequence's rows
 * a table on its system (named by its address), all of them one inline `seq` source; a row naming no sequence (or one
 * `include` leaves out) is read nowhere (counted by `checkAssemblyData`). A row's key is its row in the node's table.
 */
function inlineOverAssembly(
  at: Placed,
  a: AssemblyEntry,
  base: Pick<AnnotationEntry, 'address' | 'shape' | 'target'>,
  systems: ReadonlyMap<string, SystemEntry>,
): AnnotationEntry {
  const rows = (p(at.node, 'rows') ?? {}) as Record<string, Value>;
  const names = rows.target as ({ length: number; get?(i: number): string | null } & ArrayLike<string>) | undefined;
  const n = (rows.start as ArrayLike<number> | undefined)?.length ?? 0;
  const groups = new Map<string, number[]>();
  for (let i = 0; i < n && names; i++) {
    const name = names.get ? names.get(i) : names[i];
    const k = typeof name === 'string' ? a.find(name) : undefined;
    const system = k === undefined ? undefined : a.systemAt[k];
    if (system === undefined) continue;
    let g = groups.get(system);
    if (!g) groups.set(system, (g = []));
    g.push(i);
  }
  const listed: SystemInfo[] = [];
  const shapes: ShapeInfo[] = [];
  const pairs: SeqPair[] = [];
  const levels = new Set<number>([1]);
  const answers: { system: string; length: number; answer: SeqPair['answer'] }[] = [];
  let fields: readonly string[] = [];
  // in the assembly's order, so the inline source lists its systems as the blocks run
  for (const system of a.systems) {
    const g = groups.get(system);
    const s = systems.get(system);
    if (!g || !s) continue;
    const sub: Record<string, Value> = {};
    for (const [c, v] of Object.entries(rows)) if (c !== 'target') sub[c] = pick(v, g);
    const r = inlineRanges(sub, s);
    fields = r.fields;
    // a row's key: its row in the node's table, so two sequences' rows are never one item
    const table = r.table.row ? { ...r.table, row: Int32Array.from(r.table.row, (j) => g[j]) } : r.table;
    const one = inlineSeq({ name: system, length: s.length, cyclic: s.cyclic }, { ranges: table });
    listed.push(...one.coverage.systems);
    shapes.push(...one.coverage.shapes);
    for (const l of one.coverage.levels) levels.add(l);
    answers.push({ system, length: one.coverage.systems[0].length ?? 0, answer: one.pairs[0].answer });
  }
  // every sequence served at every level listed (its rows whole at each: the source's levels are one list)
  const sorted = [...levels].sort((x, y) => x - y);
  for (const level of sorted)
    for (const { system, length, answer } of answers)
      pairs.push({ request: { system, start: 0, end: length, level, shapes: ['ranges'] }, answer });
  const data: SeqData = { coverage: { systems: listed, shapes, levels: sorted }, pairs };
  const source = { address: at.address, format: 'seq', value: seqToValue(data) };
  const bySystem = new Map<string, AccessPlan>();
  for (const system of groups.keys())
    bySystem.set(system, { kind: 'inline', source, runtime: SeqRuntime, system, shape: 'ranges' });
  const first = bySystem.values().next().value;
  return { ...base, fields, ...(first ? { access: first } : {}), bySystem };
}
