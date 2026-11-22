/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Inline pairs, grids and matrices over an assembly (_spec/README.md "Resolve", Assemblies: data over an assembly), in
// the genome chunk: each row's sequences by its name columns (`a_name`/`b_name`, `x_name`/`y_name`, `name`) through
// the names rule, a side on a sequence its own system. Pairs are held, as a complex's links are: each end on its
// sequence's unit (`PairsData.aUnit`, `bUnit` over the entry's `units`), every column kept, binned by the pairs level
// in the reader (region-reads "Held pairs get the same levels"). Grids and matrices are split by sequence (a sequence
// pair's values one grid, a sequence's items its positions in order) into one inline `seq` source over every system
// named, each system named by its address: the same access plan a source gets. Nothing is dropped silently: names the
// assembly lacks, rows naming none and ends past a sequence's length are counted and said once (EC-189).

import { type Issue, issueOf, isColumn, type Node, type Value } from '@molstar/core/repr/tree';
import {
  DefaultMaxBlockValues,
  type MatrixItems,
  type SeqData,
  type SeqPair,
  type Shape,
  type ShapeInfo,
  type SystemInfo,
  seqToValue,
} from '@molstar/seqx-core/data';
// the runtime by its own module (the inline sources' reader): no barrel between a lazy chunk and the eager code
import { SeqRuntime } from '@molstar/seqx-core/data/seq-runtime';
import type { AssemblyEntry } from './genome.js';
import { GenomeIssues } from './genome-issues.js';
import { type EndUnits, readPairs } from './held.js';
import { ResolveIssues } from './issues.js';
import type { Placed } from './model-molecules.js';
import type { AccessPlan, AnnotationEntry, SystemEntry } from './model.js';

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/** A side of inline data: an assembly (its rows by name) or one system (every row on it). */
export interface InlineSide {
  readonly assembly?: AssemblyEntry;
  /** The side's system when it is no assembly. */
  readonly system?: string;
}

type Names = { readonly length: number; get?(i: number): string | null } & ArrayLike<string | null>;

/** What a node's rows say, each counted once a node. */
interface Said {
  /** Names the assembly lacks, in order. */
  readonly unknown: Set<string>;
  /** Rows (items, values) naming no sequence: a null, a short column, no string. */
  readonly unnamed: Set<number>;
  /** Rows (items) past their sequence's length. */
  past: number;
}

/** The name at row `i` of a column (null where it has none). */
function nameAt(c: Names | undefined, i: number): unknown {
  return c ? (c.get ? c.get(i) : c[i]) : null;
}

/**
 * A row's system on a side: by its name through the names rule, else the side's own; undefined when none. A row
 * naming no string is counted (`unnamed`); a name the assembly lacks is said, never while the assembly is pending (its
 * names are not known yet); one `include` leaves out reads nothing and is no unknown name.
 */
function systemOn(side: InlineSide, name: unknown, row: number, said: Said): string | undefined {
  const a = side.assembly;
  if (!a) return side.system;
  if (a.pending || !a.lookup.index.size) return undefined;
  if (typeof name !== 'string') {
    said.unnamed.add(row);
    return undefined;
  }
  const k = a.find(name);
  if (k === undefined) said.unknown.add(name);
  return k === undefined ? undefined : a.systemAt[k];
}

/** What inline data over an assembly resolves to: the annotation's own params, its fields, what it says. */
export interface InlineOverAssembly {
  readonly entry: Partial<AnnotationEntry>;
  readonly fields: readonly string[];
  readonly issues: Issue[];
}

/** An inline pairs, grid or matrix node over an assembly: pairs held, grids and matrices one inline `seq` source. */
export function inlineOverAssembly(
  at: Placed,
  sides: readonly [InlineSide, InlineSide?],
  systems: ReadonlyMap<string, SystemEntry>,
  states?: number,
): InlineOverAssembly {
  const said: Said = { unknown: new Set(), unnamed: new Set(), past: 0 };
  const kind = at.node.kind;
  const length = (s: string) => systems.get(s)?.length ?? 0;
  // the systems in the assemblies' sequence order (a side on a sequence first), as a source's plans are, so the blocks
  // read and gather in one order whichever way the rows came
  const ranks = new Map<string, number>();
  for (const e of new Set(sides.map((x) => x?.assembly)))
    for (const x of e?.systemAt ?? []) if (x !== undefined && !ranks.has(x)) ranks.set(x, ranks.size);
  const rank = (s: string) => ranks.get(s) ?? -1;
  const order = (names: Iterable<string>) => [...names].sort((x, y) => rank(x) - rank(y));
  const issues: Issue[] = [];
  let entry: Partial<AnnotationEntry>;
  let fields: readonly string[] = [];
  if (kind === 'pairs') {
    const held = heldPairs(at.node, sides, said, length, ranks);
    entry = { pairs: held.pairs, units: held.units };
    fields = held.fields;
  } else {
    const built: { data: SeqData; sizes?: readonly string[] } =
      kind === 'grid'
        ? gridBySystem(at.node, sides, said, length, order)
        : matrixBySystem(at.node, sides[0], said, length, order, states ?? 0);
    const shape = kind as Shape;
    const source = { address: at.address, format: 'seq', value: seqToValue(built.data) };
    const bySystem = new Map<string, AccessPlan>();
    for (const s of built.data.coverage.systems)
      bySystem.set(s.name, { kind: 'inline', source, runtime: SeqRuntime, system: s.name, shape });
    const plan = bySystem.values().next().value;
    entry = { ...(plan ? { access: plan } : {}), bySystem };
    if (kind === 'grid') fields = ['value'];
    for (const size of built.sizes ?? [])
      issues.push(
        issueOf(ResolveIssues, 'seq.resolve.grid.size', size, {
          address: at.address,
          hint: ResolveIssues.codes['seq.resolve.grid.size'].hint,
        }),
      );
  }
  const say = (code: 'seq.resolve.assembly.name' | 'seq.resolve.assembly.bounds', message: string) =>
    issues.push(issueOf(GenomeIssues, code, message, { address: at.address, hint: GenomeIssues.codes[code].hint }));
  const what = kind === 'matrix' ? 'items' : kind === 'grid' ? 'values' : 'rows';
  if (said.unknown.size)
    say(
      'seq.resolve.assembly.name',
      `Rows name ${said.unknown.size} sequences the assembly lacks (${[...said.unknown].slice(0, 3).join(', ')})`,
    );
  if (said.unnamed.size) say('seq.resolve.assembly.name', `${said.unnamed.size} ${what} name no sequence: left out`);
  if (said.past) say('seq.resolve.assembly.bounds', `${said.past} ${what} are past their sequence’s length: left out`);
  return { entry, fields, issues };
}

/**
 * Pairs held: each end's unit its sequence's index in `units` (every system of the sides, in the assemblies' order),
 * read by `readPairs` as a complex's links are, so `id`, `label`, `type`, strands and the ride-along columns (a
 * `score` a theme reads) stay; an end past its sequence's length leaves its row out, counted.
 */
function heldPairs(
  node: Node,
  sides: readonly [InlineSide, InlineSide?],
  said: Said,
  length: (s: string) => number,
  ranks: ReadonlyMap<string, number>,
) {
  const units: string[] = [];
  const index = new Map<string, number>();
  const add = (s: string) => {
    if (!index.has(s)) index.set(s, units.push(s) - 1);
  };
  for (const side of sides) if (side?.system !== undefined) add(side.system);
  for (const s of ranks.keys()) add(s);
  // an end past its sequence's length leaves its row out, counted once a row (never drawn past its block)
  const rows = (p(node, 'rows') ?? {}) as Record<string, Value>;
  const endOf = (k: string) => (isColumn(rows[k]) ? (rows[k] as unknown as ArrayLike<number>) : []);
  const aEnd = endOf('a_end');
  const bEnd = endOf('b_end');
  const past = new Set<number>();
  const ends: EndUnits = {
    columns: ['a_name', 'b_name'],
    unit: (name, side, row) => {
      const s = systemOn(side === 0 ? sides[0] : (sides[1] ?? sides[0]), name, row, said);
      if (s === undefined) return -1;
      const l = length(s);
      if (l && (side ? bEnd : aEnd)[row] > l) {
        past.add(row);
        return -1;
      }
      return index.get(s) ?? -1;
    },
  };
  const pairs = readPairs(node, ends);
  said.past += past.size;
  const fields = ['type', ...(pairs.labels ? ['label'] : []), ...Object.keys(pairs.extra ?? {})];
  return { pairs, units, fields };
}

/** One whole-system pair a system (its rows, level 1), in the systems' order. */
function seqData(
  systems: readonly string[],
  length: (s: string) => number,
  shapes: ShapeInfo[],
  pairs: SeqPair[],
): SeqData {
  const listed: SystemInfo[] = systems.map((name) => ({ name, length: length(name) }));
  return { coverage: { systems: listed, shapes, levels: [1] }, pairs };
}

/** A category column of strings, or undefined. */
function names(v: Value | undefined): Names | undefined {
  return v === undefined || v === null ? undefined : (v as unknown as Names);
}

/**
 * A grid: a sequence pair's values (in the node's order, x-major) one grid of `|x| × |y|`, served only for the pairs
 * with values (the coverage lists each x's y systems); a pair past `DefaultMaxBlockValues` is refused before anything
 * is allocated, one of the wrong count left out, each said.
 */
function gridBySystem(
  node: Node,
  sides: readonly [InlineSide, InlineSide?],
  said: Said,
  length: (s: string) => number,
  order: (names: Iterable<string>) => string[],
) {
  const values = isColumn(p(node, 'values')) ? (p(node, 'values') as unknown as ArrayLike<number>) : [];
  const xNames = names(p(node, 'x_name'));
  const yNames = names(p(node, 'y_name'));
  const groups = new Map<string, Map<string, number[]>>();
  for (let i = 0; i < values.length; i++) {
    const x = systemOn(sides[0], nameAt(xNames, i), i, said);
    const y = systemOn(sides[1] ?? sides[0], nameAt(yNames, i), i, said);
    if (x === undefined || y === undefined) continue;
    let byY = groups.get(x);
    if (!byY) groups.set(x, (byY = new Map()));
    let g = byY.get(y);
    if (!g) byY.set(y, (g = []));
    g.push(i);
  }
  const pairs: SeqPair[] = [];
  const shapes: ShapeInfo[] = [];
  const sizes: string[] = [];
  const listed = new Set<string>();
  for (const x of order(groups.keys())) {
    const byY = groups.get(x)!;
    const ys: string[] = [];
    for (const y of order(byY.keys())) {
      const at = byY.get(y)!;
      const nx = length(x);
      const ny = length(y);
      if (nx * ny > DefaultMaxBlockValues) {
        sizes.push(`A sequence pair's grid of ${nx * ny} values is past ${DefaultMaxBlockValues}`);
        continue;
      }
      if (at.length !== nx * ny) {
        sizes.push(`A sequence pair's grid has ${at.length} of ${nx * ny} values: left out`);
        continue;
      }
      ys.push(y);
      listed.add(x).add(y);
      pairs.push({
        request: {
          system: x,
          start: 0,
          end: nx,
          level: 1,
          shapes: ['grid'],
          y: { system: y, start: 0, end: ny },
        },
        answer: {
          ok: true,
          tables: { grid: { form: 'dense', x: 0, y: 0, nx, ny, value: Float32Array.from(at, (i) => values[i]) } },
        },
      });
    }
    if (ys.length) shapes.push({ system: x, shape: 'grid', y: ys });
  }
  // the y systems are listed too: a section's y block reads its length from the coverage
  return { data: seqData(order(listed), length, shapes, pairs), sizes: sizes.slice(0, 1) };
}

/**
 * A matrix: a sequence's items (in the node's order) its positions from 0, each with every state, as one matrix's
 * inline items; an item's key its index in the node, so two sequences' items are never one. Items past a sequence's
 * length are left out, counted.
 */
function matrixBySystem(
  node: Node,
  side: InlineSide,
  said: Said,
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
    const s = systemOn(side, nameAt(itemNames, i), i, said);
    if (s === undefined) continue;
    let g = bySystem.get(s);
    if (!g) bySystem.set(s, (g = []));
    if (g.length < length(s)) g.push(i);
    else said.past++;
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
