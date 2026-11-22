/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Inline pairs, grids and matrices over an assembly (_spec/README.md "Resolve", Assemblies: data over an assembly), in
// the genome chunk: each row's sequences by its name columns (`a_name`/`b_name`, `x_name`/`y_name`, `name`) through
// the names rule, a side on a sequence its own system; the rows split by sequence (pairs: each in the read of each
// system an end is on; grids: a sequence pair's values one grid; matrices: a sequence's items its positions in order)
// into one inline `seq` source over every system named, each system named by its address: the same access plan a
// source gets, so the store reads, bins (the pairs level, the grid's pyramid) and bounds them as it does a file's.

import {
  type CategoryColumn,
  categoryColumn,
  type Issue,
  issueOf,
  isColumn,
  type Node,
  type Value,
} from '@molstar/core/repr/tree';
import {
  binPairs,
  coarsenPairBins,
  DefaultSummarySettings,
  type MatrixItems,
  type PairsBins,
  type PairsList,
  type SeqData,
  type SeqPair,
  SeqRuntime,
  type Shape,
  type ShapeInfo,
  type SystemInfo,
  seqToValue,
} from '@molstar/seqx-core/data';
import type { AssemblyEntry } from './genome.js';
import { GenomeIssues } from './genome-issues.js';
import { ResolveIssues } from './issues.js';
import type { Placed } from './model-molecules.js';
import type { AccessPlan, SystemEntry } from './model.js';

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/** A side of inline data: an assembly (its rows by name) or one system (every row on it). */
export interface InlineSide {
  readonly assembly?: AssemblyEntry;
  /** The side's system when it is no assembly. */
  readonly system?: string;
}

type Names = { readonly length: number; get?(i: number): string | null } & ArrayLike<string | null>;

/** A row's system on a side: by its name through the names rule, else the side's own; undefined when none. */
function systemOn(side: InlineSide, names: Names | undefined, i: number, unknown: Set<string>): string | undefined {
  if (!side.assembly) return side.system;
  const name = names ? (names.get ? names.get(i) : names[i]) : null;
  if (typeof name !== 'string') return undefined;
  const k = side.assembly.find(name);
  // a name of the assembly by the names rule; one `include` leaves out reads nothing, and is no unknown name
  if (k === undefined) unknown.add(name);
  return k === undefined ? undefined : side.assembly.systemAt[k];
}

/** What inline data over an assembly resolves to: a plan a system named, its fields, what it says. */
export interface InlineOverAssembly {
  readonly bySystem: Map<string, AccessPlan>;
  readonly fields: readonly string[];
  readonly issues: Issue[];
}

/** The rows of an inline pairs, grid or matrix node over an assembly, by system, as one inline `seq` source. */
export function inlineOverAssembly(
  at: Placed,
  sides: readonly [InlineSide, InlineSide?],
  systems: ReadonlyMap<string, SystemEntry>,
  states?: number,
): InlineOverAssembly {
  const unknown = new Set<string>();
  const kind = at.node.kind;
  const length = (s: string) => systems.get(s)?.length ?? 0;
  // the systems listed in the assemblies' sequence order (a side on a sequence first), as a source's plans are, so
  // the blocks read and gather in one order whichever way the rows came
  const ranks = new Map<string, number>();
  for (const e of new Set(sides.map((x) => x?.assembly)))
    e?.systemAt.forEach((x, i) => x !== undefined && !ranks.has(x) && ranks.set(x, ranks.size + i));
  const rank = (s: string) => ranks.get(s) ?? -1;
  const order = (names: Iterable<string>) => [...names].sort((x, y) => rank(x) - rank(y));
  const built: { data: SeqData; sizes?: readonly string[] } =
    kind === 'pairs'
      ? pairsBySystem(at.node, sides, unknown, length, order)
      : kind === 'grid'
        ? gridBySystem(at.node, sides, unknown, length, order)
        : matrixBySystem(at.node, sides[0], unknown, length, order, states ?? 0);
  const shape = kind as Shape;
  const source = { address: at.address, format: 'seq', value: seqToValue(built.data) };
  const bySystem = new Map<string, AccessPlan>();
  for (const s of built.data.coverage.systems)
    bySystem.set(s.name, { kind: 'inline', source, runtime: SeqRuntime, system: s.name, shape });
  const issues: Issue[] = [];
  if (unknown.size)
    issues.push(
      issueOf(
        GenomeIssues,
        'seq.resolve.assembly.name',
        `Rows name ${unknown.size} sequences the assembly lacks (${[...unknown].slice(0, 3).join(', ')})`,
        { address: at.address, hint: GenomeIssues.codes['seq.resolve.assembly.name'].hint },
      ),
    );
  const sizes = built.sizes ?? [];
  if (sizes.length)
    issues.push(
      issueOf(ResolveIssues, 'seq.resolve.grid.size', `A sequence pair's grid has ${sizes[0]} values: left out`, {
        address: at.address,
        hint: ResolveIssues.codes['seq.resolve.grid.size'].hint,
      }),
    );
  return { bySystem, fields: kind === 'grid' ? ['value'] : [], issues };
}

/** One whole-system pair a system (its rows, level 1), in the systems' first-named order. */
function seqData(
  systems: readonly string[],
  length: (s: string) => number,
  shapes: ShapeInfo[],
  pairs: SeqPair[],
  levels: readonly number[] = [1],
): SeqData {
  const listed: SystemInfo[] = systems.map((name) => ({ name, length: length(name) }));
  return { coverage: { systems: listed, shapes, levels }, pairs };
}

/** A category column of strings, or undefined. */
function names(v: Value | undefined): Names | undefined {
  return v === undefined || v === null ? undefined : (v as unknown as Names);
}

const PairEnds = ['a_start', 'a_end', 'b_start', 'b_end'] as const;

/**
 * Pairs: each row on its two systems, kept in the read of each (a pair between two sequences is a chord between their
 * blocks, read from either end); its key its row in the node's table, so a selection names it at every level.
 */
function pairsBySystem(
  node: Node,
  sides: readonly [InlineSide, InlineSide?],
  unknown: Set<string>,
  length: (s: string) => number,
  order: (names: Iterable<string>) => string[],
) {
  const rows = (p(node, 'rows') ?? {}) as Record<string, Value>;
  const [as, ae, bs, be] = PairEnds.map((k) => (isColumn(rows[k]) ? (rows[k] as unknown as ArrayLike<number>) : []));
  const n = Math.min(as.length, ae.length, bs.length, be.length);
  const aNames = names(rows.a_name);
  const bNames = names(rows.b_name);
  const type = names(rows.type);
  const aStrand = names(rows.a_strand);
  const bStrand = names(rows.b_strand);
  const at = (c: Names | undefined, i: number) => (c ? (c.get ? c.get(i) : c[i]) : null);
  const strand = (c: Names | undefined, i: number) => (at(c, i) === '-' ? -1 : at(c, i) === '+' ? 1 : 0);
  const bySystem = new Map<string, number[]>();
  const aSys: string[] = [];
  const bSys: string[] = [];
  const add = (s: string, i: number) => {
    let g = bySystem.get(s);
    if (!g) bySystem.set(s, (g = []));
    if (g[g.length - 1] !== i) g.push(i);
  };
  for (let i = 0; i < n; i++) {
    const a = systemOn(sides[0], aNames, i, unknown);
    const b = systemOn(sides[1] ?? sides[0], bNames, i, unknown);
    aSys[i] = a ?? '';
    bSys[i] = b ?? '';
    // an end that is no interval, or on no sequence, leaves the row out (the checks said the one, the names the other)
    if (a === undefined || b === undefined || !(as[i] >= 0 && ae[i] > as[i] && bs[i] >= 0 && be[i] > bs[i])) continue;
    add(a, i);
    add(b, i);
  }
  const strands = aStrand !== undefined || bStrand !== undefined;
  const listOf = (g: readonly number[]): PairsList => ({
    form: 'list',
    aSystem: categoryColumn(g.map((i) => aSys[i])),
    bSystem: categoryColumn(g.map((i) => bSys[i])),
    aStart: Int32Array.from(g, (i) => as[i]),
    aEnd: Int32Array.from(g, (i) => ae[i]),
    bStart: Int32Array.from(g, (i) => bs[i]),
    bEnd: Int32Array.from(g, (i) => be[i]),
    ...(strands
      ? { aStrand: Int8Array.from(g, (i) => strand(aStrand, i)), bStrand: Int8Array.from(g, (i) => strand(bStrand, i)) }
      : {}),
    ...(type ? { type: categoryColumn(g.map((i) => at(type, i))) } : {}),
    row: Int32Array.from(g),
  });
  const systems = order(bySystem.keys());
  const shapes: ShapeInfo[] = [];
  const pairs: SeqPair[] = [];
  const request = (system: string, level: number) => ({
    system,
    start: 0,
    end: length(system),
    level,
    shapes: ['pairs' as Shape],
  });
  for (const system of systems) {
    // kinds kept by the bins, as a held table's are (`binPairs`' category)
    shapes.push({ system, shape: 'pairs', ...(type ? { binBy: 'type' } : {}) });
    pairs.push({ request: request(system, 1), answer: { ok: true, tables: { pairs: listOf(bySystem.get(system)!) } } });
  }
  // the pairs level (region-reads Decided 2, 7): bin-pair counts at levels `factor` apart up to the longest system,
  // each from the one below and every one listed, as a file's index lists its levels; a level-1 source alone would be
  // read at level 1 over a whole genome at the fit (past the store's tiles), its pyramid dropping sparse links' levels
  const all = [...new Set(systems.flatMap((s) => bySystem.get(s)!))].sort((x, y) => x - y);
  const longest = Math.max(0, ...systems.map(length));
  const f = DefaultSummarySettings.factor;
  const levels = [1];
  let bins: PairsBins | undefined;
  for (let level = f; level < longest; level *= f) {
    bins = bins ? coarsenPairBins(bins, f) : binPairs(listOf(all), level, { category: !!type, systems: true });
    levels.push(level);
    for (const system of systems)
      pairs.push({ request: request(system, level), answer: { ok: true, tables: { pairs: binsOn(bins, system) } } });
  }
  return { data: seqData(systems, length, shapes, pairs, levels) };
}

/** The bin pairs with an end on `system`. */
function binsOn(t: PairsBins, system: string): PairsBins {
  const rows: number[] = [];
  for (let i = 0; i < t.count.length; i++)
    if (t.aSystem?.get(i) === system || t.bSystem?.get(i) === system) rows.push(i);
  const take = <A extends Int32Array | Int8Array | Float32Array>(c: A): A =>
    (c.constructor as unknown as { from(a: number[]): A }).from(rows.map((i) => c[i]));
  const cat = (c: CategoryColumn) => categoryColumn(rows.map((i) => c.get(i)));
  return {
    form: 'bins',
    ...(t.aSystem ? { aSystem: cat(t.aSystem) } : {}),
    ...(t.bSystem ? { bSystem: cat(t.bSystem) } : {}),
    aBin: take(t.aBin),
    bBin: take(t.bBin),
    ...(t.strand ? { strand: take(t.strand) } : {}),
    ...(t.category ? { category: cat(t.category) } : {}),
    count: take(t.count),
    ...(t.sum ? { sum: take(t.sum) } : {}),
  };
}

/**
 * A grid: a sequence pair's values (in the node's order, x-major) one grid of `|x| × |y|`; a pair of the wrong count is
 * left out, said. Every pair of the systems named on x and y is served (one named nowhere empty), so no strip pair of
 * them is a gap.
 */
function gridBySystem(
  node: Node,
  sides: readonly [InlineSide, InlineSide?],
  unknown: Set<string>,
  length: (s: string) => number,
  order: (names: Iterable<string>) => string[],
) {
  const values = isColumn(p(node, 'values')) ? (p(node, 'values') as unknown as ArrayLike<number>) : [];
  const xNames = names(p(node, 'x_name'));
  const yNames = names(p(node, 'y_name'));
  const groups = new Map<string, { x: string; y: string; at: number[] }>();
  const xs = new Set<string>();
  const ys = new Set<string>();
  for (let i = 0; i < values.length; i++) {
    const x = systemOn(sides[0], xNames, i, unknown);
    const y = systemOn(sides[1] ?? sides[0], yNames, i, unknown);
    if (x === undefined || y === undefined) continue;
    const key = `${x}\u0000${y}`;
    let g = groups.get(key);
    if (!g) groups.set(key, (g = { x, y, at: [] }));
    g.at.push(i);
    xs.add(x);
    ys.add(y);
  }
  const pairs: SeqPair[] = [];
  const request = (x: string, y: string) => ({
    system: x,
    start: 0,
    end: length(x),
    level: 1,
    shapes: ['grid' as Shape],
    y: { system: y, start: 0, end: length(y) },
  });
  const sizes: string[] = [];
  const xo = order(xs);
  const yo = order(ys);
  for (const x of xo)
    for (const y of yo) {
      const g = groups.get(`${x}\u0000${y}`);
      const nx = length(x);
      const ny = length(y);
      const whole = g && g.at.length === nx * ny;
      if (g && !whole) sizes.push(`${g.at.length} of ${nx * ny}`);
      pairs.push({
        request: request(x, y),
        answer: {
          ok: true,
          tables: {
            grid: whole
              ? { form: 'dense', x: 0, y: 0, nx, ny, value: Float32Array.from(g.at, (i) => values[i]) }
              : sparseEmpty(),
          },
        },
      });
    }
  const shapes: ShapeInfo[] = xo.map((system) => ({ system, shape: 'grid', y: yo }));
  // the y systems are listed too: a section's y block reads its length from the coverage
  const listed = order(new Set([...xs, ...ys]));
  return { data: seqData(listed, length, shapes, pairs), sizes };
}

function sparseEmpty() {
  return {
    form: 'sparse' as const,
    x: 0,
    y: 0,
    xBin: new Int32Array(0),
    yBin: new Int32Array(0),
    value: new Float32Array(0),
    fill: Number.NaN,
  };
}

/**
 * A matrix: a sequence's items (in the node's order) its positions from 0, each with every state, as one matrix's
 * inline items; an item's key its index in the node, so two sequences' items are never one.
 */
function matrixBySystem(
  node: Node,
  side: InlineSide,
  unknown: Set<string>,
  length: (s: string) => number,
  order: (names: Iterable<string>) => string[],
  states: number,
) {
  const values = p(node, 'values');
  const codes = p(node, 'codes');
  const cells = isColumn(values) ? (values as Float32Array) : isColumn(codes) ? (codes as Uint8Array) : undefined;
  const items = cells && states ? Math.floor(cells.length / states) : 0;
  const itemNames = names(p(node, 'name'));
  const bySystem = new Map<string, number[]>();
  for (let i = 0; i < items; i++) {
    const s = systemOn(side, itemNames, i, unknown);
    if (s === undefined) continue;
    let g = bySystem.get(s);
    if (!g) bySystem.set(s, (g = []));
    // a sequence's items past its length are left out (a position each)
    if (g.length < length(s)) g.push(i);
  }
  const shapes: ShapeInfo[] = [];
  const pairs: SeqPair[] = [];
  for (const system of order(bySystem.keys())) {
    const g = bySystem.get(system)!;
    const take = <A extends Float32Array | Uint8Array>(c: A): A => {
      const out = new (c.constructor as new (n: number) => A)(g.length * states);
      g.forEach((i, j) => out.set(c.subarray(i * states, (i + 1) * states), j * states));
      return out;
    };
    const matrix: MatrixItems = {
      form: 'items',
      start: Int32Array.from(g, (_, j) => j),
      end: Int32Array.from(g, (_, j) => j + 1),
      states,
      ...(cells instanceof Float32Array ? { value: take(cells) } : {}),
      ...(cells instanceof Uint8Array ? { codes: take(cells) } : {}),
      row: Int32Array.from(g),
    };
    shapes.push({ system, shape: 'matrix' });
    pairs.push({
      request: { system, start: 0, end: length(system), level: 1, shapes: ['matrix'] },
      answer: { ok: true, tables: { matrix } },
    });
  }
  return { data: seqData(order(bySystem.keys()), length, shapes, pairs) };
}
