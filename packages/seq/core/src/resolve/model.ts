/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The static model (_spec/README.md "Resolve"): tier 1, pure, over a checked spec. Every coordinate system the spec
// declares, its maps resolved to their encodings, paths composed on demand and kept, each source bound to its runtime
// and each annotation's access plan. It implements the query forms' model, so queries run on a real spec. Entries are
// keyed by address and cached by the hashes of the nodes they read, and paths by the map entries their search read: a
// resolver given an edited spec rebuilds only what the edit touched.

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
import {
  codingLength,
  codonLinear,
  codonRuns,
  DefaultGeneticCode,
  type GeneticCode,
  GeneticCodes,
  hasUnsplitExon,
} from './codon.js';
import { ResolveIssues } from './issues.js';
import { composeLinear, type Encoding, expand, ImplicitLinear, isEmptyLinear, type Linear, toRuns } from './linear.js';
import { composeEdges, type Edge, findPath, inversePath, type Path, projectRanges, type Step } from './paths.js';

/** The length a system of unknown length (letters from a source, a streamed genome) answers the query forms with. */
export const OpenLength = 0x7fffffff;

/** The choices of the resolve pass, each with its recommended default (`DefaultResolveConfig`). */
export interface ResolveConfig {
  /** Counted maps a path takes at most: 4. */
  readonly maxPathSteps: number;
  /** Shortest candidates compared at most before a path is refused as over the limit: 16. */
  readonly maxPathCandidates: number;
  /** How bounded maps are held: the cheapest encoding that answers (default), or all as runs, or all expanded. */
  readonly encoding: 'cheapest' | 'runs' | 'expanded';
  /**
   * The length the query forms read on a system whose length is unknown, so labels resolve there before its source is
   * opened: int32's max. Layout never reads it: an axis of unknown length is pending.
   */
  readonly openQueryLength: number;
}

export const DefaultResolveConfig: ResolveConfig = {
  maxPathSteps: 4,
  maxPathCandidates: 16,
  encoding: 'cheapest',
  openQueryLength: OpenLength,
};

/** A coordinate system the spec declares: a sequence, its length when its letters are inline. */
export interface SystemEntry {
  readonly address: string;
  readonly kind: 'sequence';
  readonly alphabet: 'protein' | 'dna' | 'rna';
  /** Undefined: unknown until its source is opened (pending, never a guess). */
  readonly length?: number;
  readonly cyclic: boolean;
  readonly letters?: string;
  /** The source of its letters. */
  readonly source?: string;
  /** The name its source knows it by: its `name`, else its id. */
  readonly name: string;
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
  /** A codon map's genetic code: its `table`, 1 when absent; undefined when that code is not built. */
  readonly code?: GeneticCode;
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
  readonly config: ResolveConfig;
  readonly systems: ReadonlyMap<string, SystemEntry>;
  readonly maps: ReadonlyMap<string, MapEntry>;
  readonly sources: ReadonlyMap<string, SourceBinding>;
  readonly annotations: ReadonlyMap<string, AnnotationEntry>;
  readonly issues: readonly Issue[];
  /** The shortest unambiguous path, composed once and kept (across resolves while the maps its search read hold). */
  path(from: string, to: string): PathResult;
  /**
   * A Location carried onto a system: each part along its path (rects are not carried yet); a part with no path is
   * left out with the path's issue. With `close`, the location is first closed (`close`).
   */
  project(loc: Location, to: string, options?: { close?: boolean }): Projection;
  /**
   * Closure, one step (EC-30): nucleotides gain their residues and lift to their whole codons, through each codon map
   * between the systems named (every system when none are), along the map itself.
   */
  close(loc: Location, among?: readonly string[]): Location;
}

export interface ResolveOptions {
  /** The format runtimes the application composed; inline data always reads through the `seq` runtime. */
  readonly runtimes?: readonly FormatRuntime[];
  readonly config?: Partial<ResolveConfig>;
}

export interface Resolver {
  resolve(spec: Node): StaticModel;
}

/** A path search's outcome and the map entries of the systems it read: valid while those are the same. */
interface PathEntry {
  readonly found: { readonly path: Path } | { readonly why: 'none' | 'ambiguous' | 'limit' };
  readonly deps: ReadonlyMap<string, readonly MapEntry[]>;
}

/** A resolve pass that keeps its entries between specs, by address and the hashes of the nodes each reads. */
export function createResolver(options: ResolveOptions = {}): Resolver {
  let entries = new Map<string, unknown>();
  let last: Model | undefined;
  return {
    resolve(spec) {
      if (last && last.hash === hashNode(spec)) return last;
      const next = new Map<string, unknown>();
      const cached = <T>(key: string, make: () => T): T => {
        const v = (entries.has(key) ? entries.get(key) : make()) as T;
        next.set(key, v);
        return v;
      };
      last = build(spec, options, cached, last?.pathEntries);
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

function build(spec: Node, options: ResolveOptions, cached: Cached, carried?: ReadonlyMap<string, PathEntry>): Model {
  const config: ResolveConfig = { ...DefaultResolveConfig, ...options.config };
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
    const bound = cached(`${at.address}|${key(at.node)}`, () => bindSource(at, runtimes));
    issues.push(...bound.issues);
    sources.set(at.address, bound.binding);
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
    const entry = cached(`${at.address}|${key(at.node, from.node, to.node)}|${config.encoding}`, () =>
      mapEntry(at, systems.get(from.address)!, systems.get(to.address)!, config),
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
    const entry = cached(`${at.address}|${key(at.node, target.node)}`, () =>
      annotationEntry(at, systems.get(target.address)!, src?.address),
    );
    annotations.set(at.address, entry.annotation);
    issues.push(...entry.issues);
  }

  return new Model(spec, config, systems, maps, sources, annotations, issues, carried);
}

// ---- entries ----

const SourceKeys = ['uri', 'value', 'index', 'fields', 'assembly'] as const;

/** A source's entry from its params, bound to the first runtime that reads its format; none is an issue. */
function bindSource(at: Placed, runtimes: readonly FormatRuntime[]): { binding: SourceBinding; issues: Issue[] } {
  const format = p(at.node, 'format');
  const entry: Record<string, unknown> = { address: at.address, format: typeof format === 'string' ? format : '' };
  for (const k of SourceKeys) if (p(at.node, k) !== undefined) entry[k] = p(at.node, k);
  const bound = entry as unknown as BoundSourceEntry;
  const place = (code: 'seq.resolve.source.format' | 'seq.resolve.source.runtime') => ({
    address: at.address,
    hint: ResolveIssues.codes[code].hint,
  });
  if (typeof format !== 'string' || !format)
    return {
      binding: { address: at.address, entry: bound },
      issues: [
        issueOf(
          ResolveIssues,
          'seq.resolve.source.format',
          'The source names no format',
          place('seq.resolve.source.format'),
        ),
      ],
    };
  const runtime = runtimes.find((r) => r.formats.includes(format));
  if (runtime) return { binding: { address: at.address, entry: bound, runtime }, issues: [] };
  return {
    binding: { address: at.address, entry: bound },
    issues: [
      issueOf(
        ResolveIssues,
        'seq.resolve.source.runtime',
        'No runtime reads the format of the source',
        place('seq.resolve.source.runtime'),
      ),
    ],
  };
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

/** A source read by region, the system by the name the source knows it by. */
function sourcePlan(source: string, system: string, shape: Shape): AccessPlan {
  return { kind: 'source', source, system, shape };
}

function sequenceEntry(at: Placed, source: string | undefined): SystemEntry {
  const letters = p(at.node, 'letters');
  const name = p(at.node, 'name');
  const base = {
    address: at.address,
    kind: 'sequence' as const,
    alphabet: p(at.node, 'alphabet') as SystemEntry['alphabet'],
    cyclic: p(at.node, 'cyclic') === true,
    name: typeof name === 'string' ? name : (at.node.id ?? ''),
    numberings: ['index', 'ordinal'],
  };
  if (typeof letters === 'string') {
    const tables: Tables = { letters: { form: 'letters', start: 0, letters } };
    const length = letters.length;
    return { ...base, length, letters, access: inlinePlan(at.address, { ...base, length }, 'letters', tables) };
  }
  return { ...base, ...(source ? { source } : {}), access: sourcePlan(source ?? '', base.name, 'letters') };
}

type MapCode =
  | 'seq.resolve.map.empty'
  | 'seq.resolve.map.codon_length'
  | 'seq.resolve.map.table'
  | 'seq.resolve.length.unknown';

function mapEntry(
  at: Placed,
  from: SystemEntry,
  to: SystemEntry,
  config: ResolveConfig,
): { map: MapEntry; issues: Issue[] } {
  const issues: Issue[] = [];
  const report = (code: MapCode, message: string) =>
    issues.push(issueOf(ResolveIssues, code, message, { address: at.address, hint: ResolveIssues.codes[code].hint }));
  const method = p(at.node, 'method');
  let linear: Linear;
  let step: Step = 0;
  let code: GeneticCode | undefined;
  if (method === 'codon') {
    step = 1;
    const exons = p(at.node, 'exons') as unknown as { start: ArrayLike<number>; end: ArrayLike<number> };
    const strand = p(at.node, 'strand') === '-' ? '-' : '+';
    const size = to.cyclic ? to.length : undefined;
    const unsplit = hasUnsplitExon(exons.start, exons.end, size);
    if (unsplit) report('seq.resolve.length.unknown', 'An exon crosses the origin of a gene of unknown length');
    const coding = codingLength(codonRuns(exons.start, exons.end, strand, size));
    const n = from.length;
    // a coding sequence missing a piece is not judged until it is whole
    if (!unsplit && (coding % 3 !== 0 || (n !== undefined && coding !== 3 * n && coding !== 3 * n + 3)))
      report('seq.resolve.map.codon_length', 'The exons are not three nucleotides a residue');
    const table = p(at.node, 'table');
    code = GeneticCodes[typeof table === 'number' ? table : DefaultGeneticCode];
    if (!code) report('seq.resolve.map.table', 'The genetic code of the map is not built');
    linear = codonLinear(exons.start, exons.end, strand, size);
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
  if (config.encoding !== 'cheapest' && linear.runs())
    linear = config.encoding === 'runs' ? toRuns(linear) : expand(linear, linear.lo, linear.hi);
  const map: MapEntry = {
    address: at.address,
    kind: method === 'codon' ? 'codon' : 'offset',
    from: from.address,
    to: to.address,
    encoding: linear.encoding,
    step,
    linear,
    ...(code ? { code } : {}),
  };
  return { map, issues };
}

const Strands: Readonly<Record<string, number>> = { '+': 1, '-': -1 };

function annotationEntry(
  at: Placed,
  target: SystemEntry,
  source: string | undefined,
): { annotation: AnnotationEntry; issues: Issue[] } {
  const shape = at.node.kind as 'ranges' | 'values';
  const base = { address: at.address, shape, target: target.address };
  const bySource = (s: Shape) => ({
    annotation: { ...base, fields: [], access: sourcePlan(source ?? '', target.name, s) },
    issues: [],
  });
  if (shape === 'values') {
    const values = p(at.node, 'values');
    if (!isColumn(values)) return bySource('values');
    const tables: Tables = { values: { form: 'positions', start: 0, value: values as Float32Array } };
    return {
      annotation: { ...base, fields: [], access: inlinePlan(at.address, target, 'values', tables) },
      issues: [],
    };
  }
  const rows = p(at.node, 'rows') as Record<string, Value> | undefined;
  if (!rows) return bySource('ranges');
  const { items, table, fields, wrappedUnknown } = inlineRanges(rows, target);
  const issues = wrappedUnknown
    ? [
        issueOf(ResolveIssues, 'seq.resolve.length.unknown', 'A row wraps the origin of a sequence of unknown length', {
          address: at.address,
          hint: ResolveIssues.codes['seq.resolve.length.unknown'].hint,
        }),
      ]
    : [];
  return {
    annotation: { ...base, fields, items, access: inlinePlan(at.address, target, 'ranges', { ranges: table }) },
    issues,
  };
}

/**
 * Inline rows as items and a ranges table: a row that wraps a cyclic origin is two rows of one item (EC-9, EC-10); on
 * a sequence of unknown length it is left out, and said.
 */
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
  const length = target.length;
  let wrappedUnknown = false;
  for (let i = 0; i < start.length; i++) {
    const fields: Record<string, Atom> = {};
    if (label) fields.label = label[i] ?? null;
    if (category) fields.category = category.get(i);
    if (strand) fields.strand = strand.get(i);
    const key = id?.get(i) ?? i;
    let pieces: [number, number][] = [[start[i], end[i]]];
    if (start[i] > end[i]) {
      if (length === undefined) {
        wrappedUnknown = true;
        continue;
      }
      pieces = [
        [start[i], length],
        [0, end[i]],
      ];
    }
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
  return { items, table, fields: fieldNames, wrappedUnknown };
}

// ---- the model ----

/** A codon map read both ways along itself: for closure. */
interface CodonPaths {
  readonly map: MapEntry;
  readonly down: Path;
  readonly up: Path;
}

class Model implements StaticModel {
  readonly hash: string;
  /** The path searches of this model, by `from`, `to`; the next resolve carries them. */
  readonly pathEntries = new Map<string, PathEntry>();
  private readonly results = new Map<string, PathResult>();
  private readonly edges = new Map<string, Edge[]>();
  private readonly edgeMaps = new Map<string, MapEntry[]>();
  private readonly codons: CodonPaths[] = [];
  private readonly querySystems = new Map<string, QuerySystem>();
  private readonly queryAnnotations = new Map<string, QueryAnnotation | undefined>();

  constructor(
    readonly spec: Node,
    readonly config: ResolveConfig,
    readonly systems: ReadonlyMap<string, SystemEntry>,
    readonly maps: ReadonlyMap<string, MapEntry>,
    readonly sources: ReadonlyMap<string, SourceBinding>,
    readonly annotations: ReadonlyMap<string, AnnotationEntry>,
    readonly issues: readonly Issue[],
    private readonly carried?: ReadonlyMap<string, PathEntry>,
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
        if (!list) {
          this.edges.set(e.from, (list = []));
          this.edgeMaps.set(e.from, []);
        }
        list.push(e);
        this.edgeMaps.get(e.from)!.push(m);
      }
      if (m.step !== 0)
        this.codons.push({ map: m, down: composeEdges(m.to, [back])!, up: composeEdges(m.from, [fwd])! });
    }
  }

  /** A carried search is valid while every system it read has the same map entries. */
  private valid(entry: PathEntry): boolean {
    for (const [system, maps] of entry.deps) {
      const now = this.edgeMaps.get(system) ?? [];
      if (now.length !== maps.length) return false;
      for (let i = 0; i < now.length; i++) if (now[i] !== maps[i]) return false;
    }
    return true;
  }

  private entry(from: string, to: string): PathEntry {
    const k = `${from}\u0000${to}`;
    let e = this.pathEntries.get(k);
    if (e) return e;
    const old = this.carried?.get(k);
    if (old && this.valid(old)) e = old;
    else {
      const deps = new Map<string, readonly MapEntry[]>();
      const r = findPath(
        {
          edgesFrom: (s) => {
            if (!deps.has(s)) deps.set(s, this.edgeMaps.get(s) ?? []);
            return this.edges.get(s) ?? [];
          },
        },
        from,
        to,
        { maxSteps: this.config.maxPathSteps, maxCandidates: this.config.maxPathCandidates },
      );
      e = { found: r.ok ? { path: r.path } : { why: r.why }, deps };
    }
    this.pathEntries.set(k, e);
    return e;
  }

  path(from: string, to: string): PathResult {
    const k = `${from}\u0000${to}`;
    let r = this.results.get(k);
    if (r) return r;
    // a path read the other way is the inverse of the one found
    const back = this.pathEntries.get(`${to}\u0000${from}`);
    const e = back && !this.pathEntries.has(k) ? back : this.entry(from, to);
    if ('path' in e.found) r = { path: e === back ? inversePath(e.found.path) : e.found.path };
    else r = { issue: pathIssue(e.found.why, from, to) };
    this.results.set(k, r);
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
    for (const { map, down, up } of this.codons) {
      if (allowed && (!allowed.has(map.from) || !allowed.has(map.to))) continue;
      const nt = partOf(loc, map.to);
      if (!nt) continue;
      const residues = projectRanges(down, nt.ranges);
      const codons = projectRanges(up, residues);
      out = unionLocations(
        out,
        location({
          parts: [
            { system: map.from, ranges: residues },
            { system: map.to, ranges: codons },
          ],
        }),
      );
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
      length: e.length ?? this.config.openQueryLength,
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
    const a = e?.shape === 'ranges' ? queryAnnotation(e, this.systems.get(e.target), this.config) : undefined;
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

/** A path's issue, at the system it starts from, with both systems as its places. */
function pathIssue(why: 'none' | 'ambiguous' | 'limit', from: string, to: string): Issue {
  const code =
    why === 'none'
      ? 'seq.resolve.path.none'
      : why === 'limit'
        ? 'seq.resolve.path.limit'
        : 'seq.resolve.path.ambiguous';
  const message =
    why === 'none'
      ? 'No path of maps joins the systems'
      : why === 'limit'
        ? 'More shortest paths join the systems than are compared'
        : 'Several shortest paths join the systems and map them differently';
  return {
    ...issueOf(ResolveIssues, code, message, { address: from, hint: ResolveIssues.codes[code].hint }),
    places: [{ address: from }, { address: to }],
  };
}

/** Inline items answer at once; a source's are pending over what was asked (no data in tier 1). */
function queryAnnotation(e: AnnotationEntry, target: SystemEntry | undefined, config: ResolveConfig): QueryAnnotation {
  const length = target?.length ?? config.openQueryLength;
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
