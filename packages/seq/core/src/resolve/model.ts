/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The static model (_spec/README.md "Resolve"): tier 1, pure, over a checked spec. Every coordinate system the spec
// declares, its maps resolved to their encodings, paths composed on demand and kept, each source bound to its runtime
// and each annotation's access plan. It implements the query forms' model, so queries run on a real spec. Entries are
// keyed by address and cached by the hashes of the nodes they read: a resolver given an edited spec rebuilds only
// what the edit touched.

import type { ParamsFor } from '@molstar/core/repr/builder';
import {
  type Atom,
  categoryColumn,
  exactAddress,
  hashNode,
  type Issue,
  issueOf,
  isColumn,
  type Node,
  printAddress,
  type Value,
} from '@molstar/core/repr/tree';
import {
  type FormatRuntime,
  inlineSeq,
  type RangesField,
  type RangesList,
  SeqRuntime,
  type Shape,
  type SourceEntry,
  seqToValue,
  type Tables,
} from '@molstar/seqx-core/data';
import type { SeqKinds } from '@molstar/seqx-core/lang';
import {
  type ItemsRead,
  type Location,
  location,
  overlapsRanges,
  partOf,
  type QueryAnnotation,
  type QueryItem,
  type QueryModel,
  type QuerySystem,
  rangesLocation,
  unionLocations,
} from '@molstar/seqx-core/query';
import { codingLength, codonRuns, codonLinear } from './codon.js';
import { ResolveIssues } from './issues.js';
import { composeLinear, type Encoding, ImplicitLinear, isEmptyLinear, type Linear } from './linear.js';
import {
  DefaultPathBounds,
  type Edge,
  findPath,
  inversePath,
  type Path,
  type PathBounds,
  projectRanges,
  type Step,
} from './paths.js';

/** The length a system of unknown length (letters from a source, a streamed genome) answers the query forms with. */
export const OpenLength = 0x7fffffff;

/** A coordinate system the spec declares: a sequence, its length when its letters are inline. */
export interface SystemEntry {
  readonly address: string;
  readonly kind: 'sequence';
  readonly alphabet: 'protein' | 'dna' | 'rna';
  /** Undefined: open, known only once its source is opened. */
  readonly length?: number;
  readonly cyclic: boolean;
  readonly letters?: string;
  /** The source of its letters. */
  readonly source?: string;
  /** Its numberings: the built-ins (`index`, `ordinal`) until the language declares others. */
  readonly numberings: readonly string[];
  /** How its letters are read. */
  readonly access: AccessPlan;
}

/** A map, resolved: its kind, its encoding, and its linear part (fine places of the residue side for a codon map). */
export interface MapEntry {
  readonly address: string;
  readonly kind: 'offset' | 'codon';
  readonly from: string;
  readonly to: string;
  readonly encoding: Encoding;
  readonly step: Step;
  readonly linear: Linear;
}

/** The language's params of a source: one definition, the description's. */
export type SourceParams = ParamsFor<SeqKinds, 'source'>;

/** A source as the protocol reads it, with the rest of its params. */
export type BoundSourceEntry = SourceEntry & Omit<SourceParams, 'format' | 'uri' | 'value'>;

/** A source and the runtime bound to it by its format; unavailable without one. */
export interface SourceBinding {
  readonly address: string;
  readonly entry: BoundSourceEntry;
  readonly runtime?: FormatRuntime;
}

/**
 * How a node's data is read: inline data as a `seq` source of one tile (the whole system), or a source of the spec by
 * region. `system` is the name the source knows the system by.
 */
export type AccessPlan =
  | {
      readonly kind: 'inline';
      readonly source: SourceEntry;
      readonly runtime: FormatRuntime;
      readonly system: string;
      readonly shape: Shape;
    }
  | { readonly kind: 'source'; readonly source: string; readonly system: string; readonly shape: Shape };

/** An annotation: its shape, the system it is on, its fields and how it is read. */
export interface AnnotationEntry {
  readonly address: string;
  readonly shape: 'ranges' | 'values';
  readonly target: string;
  readonly fields: readonly string[];
  readonly access: AccessPlan;
  /** Inline ranges as items (a wrapped row is two rows of one item, EC-9). */
  readonly items?: readonly QueryItem[];
}

/** A path between two systems, or why there is none. */
export interface PathResult {
  readonly path?: Path;
  readonly issue?: Issue;
}

export interface Projection {
  readonly location: Location;
  readonly issues: readonly Issue[];
}

export interface StaticModel extends QueryModel {
  readonly spec: Node;
  /** The spec's hash: equal specs give the same model. */
  readonly hash: string;
  readonly systems: ReadonlyMap<string, SystemEntry>;
  readonly maps: ReadonlyMap<string, MapEntry>;
  readonly sources: ReadonlyMap<string, SourceBinding>;
  readonly annotations: ReadonlyMap<string, AnnotationEntry>;
  readonly issues: readonly Issue[];
  /** The shortest unambiguous path, composed once and kept. */
  path(from: string, to: string): PathResult;
  /**
   * A Location carried onto a system: each part along its path (rects are not carried yet); a part with no path is
   * left out with the path's issue. With `close`, the location is first closed (`close`).
   */
  project(loc: Location, to: string, options?: { close?: boolean }): Projection;
  /**
   * Closure, one step (EC-30): nucleotides gain their residues and lift to their whole codons, through the codon
   * maps between the systems named (every system when none are).
   */
  close(loc: Location, among?: readonly string[]): Location;
}

export interface ResolveOptions {
  /** The format runtimes the application composed; inline data always reads through the `seq` runtime. */
  readonly runtimes?: readonly FormatRuntime[];
  readonly bounds?: PathBounds;
}

export interface Resolver {
  resolve(spec: Node): StaticModel;
}

/** A resolve pass that keeps its entries between specs, by address and the hashes of the nodes each reads. */
export function createResolver(options: ResolveOptions = {}): Resolver {
  let entries = new Map<string, unknown>();
  let last: StaticModel | undefined;
  return {
    resolve(spec) {
      if (last && last.hash === hashNode(spec)) return last;
      const next = new Map<string, unknown>();
      const cached = <T>(key: string, make: () => T): T => {
        const v = (entries.has(key) ? entries.get(key) : make()) as T;
        next.set(key, v);
        return v;
      };
      last = build(spec, options, cached);
      entries = next;
      return last;
    },
  };
}

/** The static model of a checked spec, with no cache. */
export function resolveSpec(spec: Node, options: ResolveOptions = {}): StaticModel {
  return build(spec, options, (_, make) => make());
}

type Cached = <T>(key: string, make: () => T) => T;

interface Placed {
  readonly node: Node;
  readonly address: string;
}

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

function build(spec: Node, options: ResolveOptions, cached: Cached): StaticModel {
  const issues: Issue[] = [];
  const data: Placed[] = [];
  const byAddress = new Map<string, Placed>();
  for (const c of spec.children ?? []) {
    if (c.kind !== 'data' || c.id === undefined) continue;
    for (const n of c.children ?? []) {
      if (n.id === undefined) continue;
      const placed = { node: n, address: printAddress([c.id, n.id]) };
      data.push(placed);
      byAddress.set(placed.address, placed);
    }
  }
  const ref = (at: Placed, key: string, kinds: readonly string[]): Placed | undefined => {
    const a = exactAddress(p(at.node, key));
    const t = a === undefined ? undefined : byAddress.get(a);
    if (t && kinds.includes(t.node.kind)) return t;
    issues.push(
      issueOf(ResolveIssues, 'seq.resolve.ref.missing', 'The reference names no node of a kind it takes', {
        address: at.address,
        path: ['params', key],
        hint: ResolveIssues.codes['seq.resolve.ref.missing'].hint,
      }),
    );
    return undefined;
  };
  const key = (...nodes: Node[]) => nodes.map(hashNode).join(':');
  const runtimes = options.runtimes ?? [];

  // sources first: what reads by region names them
  const sources = new Map<string, SourceBinding>();
  for (const at of data) {
    if (at.node.kind !== 'source') continue;
    const binding = cached(`${at.address}|${key(at.node)}`, () => bindSource(at, runtimes));
    if (!binding.runtime)
      issues.push(
        issueOf(ResolveIssues, 'seq.resolve.source.runtime', 'No runtime reads the format of the source', {
          address: at.address,
          hint: ResolveIssues.codes['seq.resolve.source.runtime'].hint,
        }),
      );
    sources.set(at.address, binding);
  }

  const systems = new Map<string, SystemEntry>();
  for (const at of data) {
    if (at.node.kind !== 'sequence') continue;
    const src = p(at.node, 'source') === undefined ? undefined : ref(at, 'source', ['source']);
    systems.set(
      at.address,
      cached(`${at.address}|${key(at.node)}`, () => sequenceEntry(at, src?.address)),
    );
  }

  const maps = new Map<string, MapEntry>();
  for (const at of data) {
    if (at.node.kind !== 'map') continue;
    const from = ref(at, 'from', ['sequence']);
    const to = ref(at, 'to', ['sequence']);
    if (!from || !to) continue;
    const entry = cached(`${at.address}|${key(at.node, from.node, to.node)}`, () =>
      mapEntry(at, systems.get(from.address)!, systems.get(to.address)!),
    );
    maps.set(at.address, entry.map);
    issues.push(...entry.issues);
  }

  const annotations = new Map<string, AnnotationEntry>();
  for (const at of data) {
    if (at.node.kind !== 'ranges' && at.node.kind !== 'values') continue;
    const target = ref(at, 'target', ['sequence']);
    if (!target) continue;
    const src = p(at.node, 'source') === undefined ? undefined : ref(at, 'source', ['source']);
    annotations.set(
      at.address,
      cached(`${at.address}|${key(at.node, target.node)}`, () =>
        annotationEntry(at, systems.get(target.address)!, src?.address),
      ),
    );
  }

  return new Model(spec, systems, maps, sources, annotations, issues, options.bounds ?? DefaultPathBounds);
}

// ---- entries ----

function bindSource(at: Placed, runtimes: readonly FormatRuntime[]): SourceBinding {
  const params = { ...(at.node.params ?? {}) } as unknown as SourceParams;
  const entry: BoundSourceEntry = { ...params, address: at.address };
  const runtime = runtimes.find((r) => r.formats.includes(entry.format));
  return runtime ? { address: at.address, entry, runtime } : { address: at.address, entry };
}

/** Inline tables as a `seq` source of one tile over the whole system, named by the system's address. */
function inlinePlan(
  address: string,
  system: { readonly address: string; readonly length?: number; readonly cyclic: boolean },
  shape: Shape,
  tables: Tables,
): AccessPlan {
  const name = system.address;
  const value = seqToValue(inlineSeq({ name, length: system.length, cyclic: system.cyclic }, tables));
  return { kind: 'inline', source: { address, format: 'seq', value }, runtime: SeqRuntime, system: name, shape };
}

/** A source read by region; the system is named by the node's id (a source with one system per sequence). */
function sourcePlan(source: string, address: string, shape: Shape): AccessPlan {
  return { kind: 'source', source, system: address.slice(address.lastIndexOf('/') + 1), shape };
}

function sequenceEntry(at: Placed, source: string | undefined): SystemEntry {
  const letters = p(at.node, 'letters');
  const base = {
    address: at.address,
    kind: 'sequence' as const,
    alphabet: p(at.node, 'alphabet') as SystemEntry['alphabet'],
    cyclic: p(at.node, 'cyclic') === true,
    numberings: ['index', 'ordinal'],
  };
  if (typeof letters === 'string') {
    const tables: Tables = { letters: { form: 'letters', start: 0, letters } };
    const length = letters.length;
    return { ...base, length, letters, access: inlinePlan(at.address, { ...base, length }, 'letters', tables) };
  }
  return { ...base, ...(source ? { source } : {}), access: sourcePlan(source ?? '', at.address, 'letters') };
}

function mapEntry(at: Placed, from: SystemEntry, to: SystemEntry): { map: MapEntry; issues: Issue[] } {
  const issues: Issue[] = [];
  const report = (code: 'seq.resolve.map.empty' | 'seq.resolve.map.codon_length', message: string) =>
    issues.push(issueOf(ResolveIssues, code, message, { address: at.address, hint: ResolveIssues.codes[code].hint }));
  const method = p(at.node, 'method');
  let linear: Linear;
  let step: Step = 0;
  if (method === 'codon') {
    step = 1;
    const exons = p(at.node, 'exons') as unknown as { start: ArrayLike<number>; end: ArrayLike<number> };
    const strand = p(at.node, 'strand') === '-' ? '-' : '+';
    const coding = codingLength(codonRuns(exons.start, exons.end, strand));
    const n = from.length;
    if (coding % 3 !== 0 || (n !== undefined && coding !== 3 * n && coding !== 3 * n + 3))
      report('seq.resolve.map.codon_length', 'The exons are not three nucleotides a residue');
    linear = codonLinear(exons.start, exons.end, strand);
    // what lies past either sequence maps nowhere (EC-21)
    if (to.length !== undefined) linear = composeLinear(linear, new ImplicitLinear(0, 1, 0, to.length));
    if (n !== undefined) linear = composeLinear(new ImplicitLinear(0, 1, 0, 3 * n), linear);
  } else {
    const by = Number(p(at.node, 'by') ?? 0);
    const lf = from.length ?? Infinity;
    if (to.cyclic && to.length !== undefined && Number.isFinite(lf))
      linear = new ImplicitLinear(by, 1, 0, lf, to.length);
    else linear = new ImplicitLinear(by, 1, Math.max(0, -by), Math.min(lf, (to.length ?? Infinity) - by));
  }
  if (isEmptyLinear(linear)) report('seq.resolve.map.empty', 'The map pairs no position');
  const map: MapEntry = {
    address: at.address,
    kind: method === 'codon' ? 'codon' : 'offset',
    from: from.address,
    to: to.address,
    encoding: linear.encoding,
    step,
    linear,
  };
  return { map, issues };
}

const Strands: Readonly<Record<string, number>> = { '+': 1, '-': -1 };

function annotationEntry(at: Placed, target: SystemEntry, source: string | undefined): AnnotationEntry {
  const shape = at.node.kind as 'ranges' | 'values';
  const base = { address: at.address, shape, target: target.address };
  if (shape === 'values') {
    const values = p(at.node, 'values');
    if (!isColumn(values)) return { ...base, fields: [], access: sourcePlan(source ?? '', target.address, 'values') };
    const tables: Tables = { values: { form: 'positions', start: 0, value: values as Float32Array } };
    return { ...base, fields: [], access: inlinePlan(at.address, target, 'values', tables) };
  }
  const rows = p(at.node, 'rows') as Record<string, Value> | undefined;
  if (!rows) return { ...base, fields: [], access: sourcePlan(source ?? '', target.address, 'ranges') };
  const { items, table, fields } = inlineRanges(rows, target);
  return { ...base, fields, items, access: inlinePlan(at.address, target, 'ranges', { ranges: table }) };
}

/** Inline rows as items and a ranges table: a row that wraps a cyclic origin is two rows of one item (EC-9, EC-10). */
function inlineRanges(rows: Record<string, Value>, target: SystemEntry) {
  const start = rows.start as Int32Array;
  const end = rows.end as Int32Array;
  const id = rows.id as { get(i: number): string | null } | undefined;
  const strand = rows.strand as { get(i: number): string | null } | undefined;
  const label = Array.isArray(rows.label) ? (rows.label as (string | null)[]) : undefined;
  const category = rows.category as { get(i: number): string | null } | undefined;
  const fieldNames = [label && 'label', category && 'category', strand && 'strand'].filter((f): f is string => !!f);
  const items: QueryItem[] = [];
  const s: number[] = [];
  const e: number[] = [];
  const rowOf: number[] = [];
  const length = target.length ?? OpenLength;
  for (let i = 0; i < start.length; i++) {
    const fields: Record<string, Atom> = {};
    if (label) fields.label = label[i] ?? null;
    if (category) fields.category = category.get(i);
    if (strand) fields.strand = strand.get(i);
    const key = id?.get(i) ?? i;
    const pieces: [number, number][] =
      start[i] > end[i]
        ? [
            [start[i], length],
            [0, end[i]],
          ]
        : [[start[i], end[i]]];
    for (const [a, b] of pieces) {
      items.push({ key, start: a, end: b, ...(fieldNames.length ? { fields } : {}) });
      s.push(a);
      e.push(b);
      rowOf.push(i);
    }
  }
  const pick = <T>(f: (i: number) => T) => rowOf.map(f);
  const columns: RangesField[] = [];
  if (label) columns.push({ name: 'label', column: categoryColumn(pick((i) => label[i] ?? null)) });
  if (category) columns.push({ name: 'category', column: categoryColumn(pick((i) => category.get(i))) });
  const table: RangesList = {
    form: 'list',
    start: Int32Array.from(s),
    end: Int32Array.from(e),
    ...(strand ? { strand: Int8Array.from(pick((i) => Strands[strand.get(i) ?? ''] ?? 0)) } : {}),
    ...(id ? { id: categoryColumn(pick((i) => id.get(i))) } : { row: Int32Array.from(rowOf) }),
    ...(columns.length ? { fields: columns } : {}),
  };
  return { items, table, fields: fieldNames };
}

// ---- the model ----

class Model implements StaticModel {
  readonly hash: string;
  private readonly paths = new Map<string, PathResult>();
  private readonly edges = new Map<string, Edge[]>();
  private readonly querySystems = new Map<string, QuerySystem>();
  private readonly queryAnnotations = new Map<string, QueryAnnotation | undefined>();

  constructor(
    readonly spec: Node,
    readonly systems: ReadonlyMap<string, SystemEntry>,
    readonly maps: ReadonlyMap<string, MapEntry>,
    readonly sources: ReadonlyMap<string, SourceBinding>,
    readonly annotations: ReadonlyMap<string, AnnotationEntry>,
    readonly issues: readonly Issue[],
    private readonly bounds: PathBounds,
  ) {
    this.hash = hashNode(spec);
    for (const m of maps.values()) {
      const fwd: Edge = { map: m.address, from: m.from, to: m.to, backward: false, step: m.step, linear: m.linear };
      const back: Edge = {
        map: m.address,
        from: m.to,
        to: m.from,
        backward: true,
        step: -m.step as Step,
        linear: m.linear.inverse(),
      };
      for (const e of [fwd, back]) {
        let list = this.edges.get(e.from);
        if (!list) this.edges.set(e.from, (list = []));
        list.push(e);
      }
    }
  }

  path(from: string, to: string): PathResult {
    const k = `${from}\u0000${to}`;
    let r = this.paths.get(k);
    if (r) return r;
    const back = this.paths.get(`${to}\u0000${from}`);
    if (back) r = back.path ? { path: inversePath(back.path) } : back;
    else {
      const found = findPath({ edgesFrom: (s) => this.edges.get(s) ?? [] }, from, to, this.bounds);
      if (found.ok) r = { path: found.path };
      else {
        const code = found.why === 'none' ? 'seq.resolve.path.none' : 'seq.resolve.path.ambiguous';
        const message =
          found.why === 'none'
            ? 'No path of maps joins the systems'
            : found.why === 'limit'
              ? 'More shortest paths join the systems than are compared'
              : 'Several shortest paths join the systems and map them differently';
        r = { issue: issueOf(ResolveIssues, code, message, { address: to, hint: ResolveIssues.codes[code].hint }) };
      }
    }
    this.paths.set(k, r);
    return r;
  }

  project(loc: Location, to: string, options: { close?: boolean } = {}): Projection {
    const src = options.close ? this.close(loc) : loc;
    const issues: Issue[] = [];
    const parts: { system: string; ranges: number[] }[] = [];
    for (const part of src.parts) {
      const r = this.path(part.system, to);
      if (!r.path) {
        if (r.issue) issues.push(r.issue);
        continue;
      }
      parts.push({ system: to, ranges: projectRanges(r.path, part.ranges) });
    }
    return { location: location({ parts }), issues };
  }

  close(loc: Location, among?: readonly string[]): Location {
    const allowed = among ? new Set(among) : undefined;
    let out = loc;
    for (const m of this.maps.values()) {
      if (m.step === 0 || (allowed && (!allowed.has(m.from) || !allowed.has(m.to)))) continue;
      const nt = partOf(loc, m.to);
      if (!nt) continue;
      const residues = this.project(rangesLocation(m.to, nt.ranges), m.from).location;
      const codons = this.project(residues, m.to).location;
      out = unionLocations(out, unionLocations(residues, codons));
    }
    return out;
  }

  // ---- the query forms' model ----

  system(address: string): QuerySystem | undefined {
    let s = this.querySystems.get(address);
    if (s) return s;
    const e = this.systems.get(address);
    if (!e) return undefined;
    s = {
      address,
      length: e.length ?? OpenLength,
      cyclic: e.cyclic,
      numbering: () => undefined,
      numberingIds: () => [],
    };
    this.querySystems.set(address, s);
    return s;
  }

  annotation(address: string): QueryAnnotation | undefined {
    if (this.queryAnnotations.has(address)) return this.queryAnnotations.get(address);
    const e = this.annotations.get(address);
    const a = e?.shape === 'ranges' ? queryAnnotation(e, this.systems.get(e.target)) : undefined;
    this.queryAnnotations.set(address, a);
    return a;
  }

  /** Dynamic nodes are the overlay's: none in the static model. */
  node(): Node | undefined {
    return undefined;
  }

  /** Structures are not built yet: no chain, and the form reports the node missing. */
  chain(): Node | undefined {
    return undefined;
  }
}

/** Inline items answer at once; a source's are pending over what was asked (no data in tier 1). */
function queryAnnotation(e: AnnotationEntry, target: SystemEntry | undefined): QueryAnnotation {
  const length = target?.length ?? OpenLength;
  return {
    address: e.address,
    target: e.target,
    fields: e.fields,
    items(over?: Location): ItemsRead {
      if (!e.items) return { items: [], pending: over ?? rangesLocation(e.target, [0, length]) };
      const part = over && partOf(over, e.target);
      if (!over) return { items: e.items };
      if (!part) return { items: [] };
      return { items: e.items.filter((it) => overlapsRanges(part.ranges, it.start, it.end, length)) };
    },
  };
}
