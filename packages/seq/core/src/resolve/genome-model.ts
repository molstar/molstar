/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What resolve reads of data over assemblies, in the genome chunk (_spec/README.md "Resolve", Assemblies): a source
// naming another assembly, inline rows naming no sequence, and a map between assemblies as its pairs of sequences.

import { exactAddress, type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import { assemblySystem } from './assemblies.js';
import type { AssemblyEntry } from './genome.js';
import { ResolveIssues } from './issues.js';
import type { Placed } from './model-molecules.js';
import { overAssembly, type SystemEntry } from './model.js';

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/**
 * What resolve says of data over an assembly: a source naming another assembly (`seq.resolve.assembly.mismatch`), and
 * inline rows naming no sequence of it, counted once with the first (`seq.resolve.assembly.name`).
 */
export function checkAssemblyData(
  data: readonly Placed[],
  assemblies: ReadonlyMap<string, AssemblyEntry>,
  byAddress: ReadonlyMap<string, Placed>,
  issues: Issue[],
): void {
  const say = (code: 'seq.resolve.assembly.mismatch' | 'seq.resolve.assembly.name', message: string, at: Placed) =>
    issues.push(issueOf(ResolveIssues, code, message, { address: at.address, hint: ResolveIssues.codes[code].hint }));
  for (const at of data) {
    if (!overAssembly(at.node, assemblies)) continue;
    const a = assemblies.get(exactAddress(p(at.node, 'target') ?? p(at.node, 'a') ?? p(at.node, 'x')) ?? '');
    const src = byAddress.get(exactAddress(p(at.node, 'source')) ?? '')?.node;
    const says = src && p(src, 'assembly');
    if (a && typeof says === 'string' && a.names.length && !a.names.includes(says))
      say('seq.resolve.assembly.mismatch', `The source is on ${says}, not ${a.names[0]}`, at);
    const target = (p(at.node, 'rows') as Record<string, Value> | undefined)?.target;
    if (!a || a.pending || !a.lookup.index.size || !target) continue;
    const names = target as { length: number; get?(i: number): string | null } & ArrayLike<string>;
    const unknown = new Set<string>();
    for (let i = 0; i < names.length; i++) {
      const name = names.get ? names.get(i) : names[i];
      // a name of the assembly by the names rule (an `include` leaving its sequence out is no unknown name)
      if (typeof name === 'string' && a.find(name) === undefined) unknown.add(name);
    }
    if (unknown.size)
      say(
        'seq.resolve.assembly.name',
        `Rows name ${unknown.size} sequences the assembly lacks (${[...unknown].slice(0, 3).join(', ')})`,
        at,
      );
  }
}

/**
 * A map between assemblies (or an assembly and a sequence) as its pairs of sequences: its inline runs grouped by their
 * `from_name` and `to_name` (an end on a sequence needs none), each pair one entry over the two systems. A run naming
 * no sequence of its side is left out (the map's `seq.resolve.assembly.name`). A map read from a source is the store's.
 */
export function assemblyPairs(
  at: Placed,
  from: AssemblyEntry | undefined,
  to: AssemblyEntry | undefined,
  systems: ReadonlyMap<string, SystemEntry>,
  issues: Issue[],
): { key: string; from: string; to: string; at: Placed }[] {
  const table = p(at.node, 'segments') as Record<string, Value> | undefined;
  if (!table || p(at.node, 'method') !== 'segments') return [];
  const column = (v: Value | undefined) => v as { length: number; get?(i: number): string | null } | undefined;
  const fromNames = column(table.from_name);
  const toNames = column(table.to_name);
  const side = (e: AssemblyEntry | undefined, names: typeof fromNames, end: string, i: number) => {
    if (!e) return systems.has(end) ? end : undefined;
    const name = names?.get?.(i) ?? (names as unknown as string[] | undefined)?.[i];
    return typeof name === 'string' ? assemblySystem(e, name) : undefined;
  };
  const rows = (table.from as ArrayLike<number> | undefined)?.length ?? 0;
  const groups = new Map<string, { from: string; to: string; rows: number[] }>();
  let unknown = 0;
  for (let i = 0; i < rows; i++) {
    const f = side(from, fromNames, exactAddress(p(at.node, 'from')) ?? '', i);
    const t = side(to, toNames, exactAddress(p(at.node, 'to')) ?? '', i);
    if (f === undefined || t === undefined) {
      unknown++;
      continue;
    }
    const k = `${f}|${t}`;
    let g = groups.get(k);
    if (!g) groups.set(k, (g = { from: f, to: t, rows: [] }));
    g.rows.push(i);
  }
  if (unknown && !(from?.pending || to?.pending))
    issues.push(
      issueOf(ResolveIssues, 'seq.resolve.assembly.name', `${unknown} runs name no sequence of their assembly`, {
        address: at.address,
        hint: ResolveIssues.codes['seq.resolve.assembly.name'].hint,
      }),
    );
  return [...groups].map(([k, g]) => {
    const segments: Record<string, Value> = {};
    for (const [c, v] of Object.entries(table)) segments[c] = pick(v, g.rows);
    return {
      key: k,
      from: g.from,
      to: g.to,
      at: { address: at.address, node: { ...at.node, params: { ...at.node.params, segments } } },
    };
  });
}

/** A table column's rows `rows`, in its own form (a typed array, a category column read out, a list). */
function pick(v: Value, rows: readonly number[]): Value {
  const c = v as unknown as {
    get?(i: number): unknown;
    constructor: new (n: number) => ArrayLike<number> & Record<number, number>;
  };
  if (ArrayBuffer.isView(v)) {
    const out = new c.constructor(rows.length);
    rows.forEach((r, j) => (out[j] = (v as unknown as ArrayLike<number>)[r]));
    return out as unknown as Value;
  }
  if (typeof c.get === 'function') return rows.map((r) => c.get!(r)) as Value;
  return Array.isArray(v) ? (rows.map((r) => v[r]) as Value) : v;
}
